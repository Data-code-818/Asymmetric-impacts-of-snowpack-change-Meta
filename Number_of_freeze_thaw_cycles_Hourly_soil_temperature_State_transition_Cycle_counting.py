# -*- coding: utf-8 -*-
r"""
Freeze-thaw metrics from hourly soil temperature, January-July, 2000-2025.
Python 3.9+; numpy, pandas, xarray, netCDF4, rasterio, eccodes, h5netcdf,
h5py, and scipy. Complete download_FTC_boundary_52days.py before processing.
Version 3 rounds longitude before normalization to resolve the 180-degree boundary.
Outputs are saved under v3_fixed; earlier-version results are retained.

Local mean solar time = UTC + longitude/15 hours (no equation-of-time correction).
Daily freeze-thaw occurrence requires 24 valid hours and Tmin < 0 < Tmax.
Transitions require three consecutive hourly samples strictly below or above 0 C;
event time is the first sample. Zero interrupts candidates but retains confirmed
states; missing data reset states and incomplete cycles. A complete cycle requires
a thaw-to-freeze transition followed by a freeze-to-thaw transition, both starting
within the target local-date interval. Initially frozen conditions are not a start.

GRIB inputs take priority; daily NetCDF files fill gaps. Decodable, complete GRIB
messages in .part files are eligible. Alternative files fill failed reads and pixel
NaNs; each pixel-hour is counted once, without interpolation or zero imputation.
Input files are never modified. Pause downloads before running; input changes abort
processing. Read 24 hours beyond each period boundary where data are available.

Counts describe valid observations and may underestimate totals with incomplete
coverage. Quality rasters report valid hours, complete local days, and missing hours.
DFTC is NoData without complete days; transition metrics are NoData without valid
hours (NoData = 65535). File size and modification time govern index-cache reuse;
annual completion markers also validate input inventories and configuration.
Only one hemispheric hourly field is decoded at a time; state arrays require memory.
Run with --self-test for algorithm checks without reading input data.
"""
from __future__ import annotations
import os
import sys
import gc
import json
import csv
import hashlib
import logging
from pathlib import Path
from dataclasses import dataclass
from collections import defaultdict, OrderedDict
from typing import Dict, Tuple
import numpy as np
import pandas as pd

GRIB_ROOT = Path(r"F:\Shen\stl1_grib_7day")
NC_ROOT = Path(r"F:\Shen\stl1")
OUTPUT_DIR = Path(r"F:\Shen\FTC_hourly_soil_temperature")
START_YEAR, END_YEAR = 2000, 2025
MONTHS = tuple(range(1, 8))
HOLD_HOURS = 3
MIN_VALID_HOURS_PER_LOCAL_DAY = 24
NODATA_UINT16 = 65535
OVERWRITE_EXISTING = False
QUICK_TEST = False  # True: 2000 only, 30-35 N and 100-105 E; separate quick_test directory.
DAY_NS = 86_400_000_000_000
HOUR_NS = 3_600_000_000_000
ZERO_C_K = np.float32(273.15)
SCRIPT_VERSION = "local-ftc-3.0"
GRIB_INDEX_VERSION = "local-ftc-1.0"  # Reuse the existing GRIB message index.
NC_INDEX_VERSION = "unicode-filelike-v2"
GRID_STEP = 0.1
# Northern Hemisphere grid; exact coordinate matching without interpolation.
GRID_LAT_MIN, GRID_LAT_MAX = 0.0, 90.0
GRID_LON_MIN, GRID_LON_MAX = -180.0, 179.9


def period_bounds(year: int) -> Tuple[pd.Timestamp, pd.Timestamp]:
    """Return the local-date interval [start, end); MONTHS must be consecutive."""
    start = pd.Timestamp(year=year, month=MONTHS[0], day=1)
    if MONTHS[-1] == 12:
        end = pd.Timestamp(year=year + 1, month=1, day=1)
    else:
        end = pd.Timestamp(year=year, month=MONTHS[-1] + 1, day=1)
    return start, end


def timestamp_ns(ts: pd.Timestamp) -> int:
    return int(ts.to_datetime64().astype("datetime64[ns]").astype(np.int64))


def day_number(ts: pd.Timestamp) -> int:
    return int(ts.to_datetime64().astype("datetime64[D]").astype(np.int64))


def expected_output_paths(output_dir: Path, year: int) -> Dict[str, Path]:
    period_tag = f"M{MONTHS[0]:02d}-M{MONTHS[-1]:02d}"
    hold_tag = f"H{HOLD_HOURS}"
    return {
        "dftc": output_dir / "01_DFTC_days" / f"DFTC_days_{year}_{period_tag}.tif",
        "freeze": output_dir
        / "02_freeze_transitions"
        / f"freeze_transitions_{hold_tag}_{year}_{period_tag}.tif",
        "thaw": output_dir
        / "03_thaw_transitions"
        / f"thaw_transitions_{hold_tag}_{year}_{period_tag}.tif",
        "total": output_dir
        / "04_FT_transitions_total"
        / f"FT_transitions_total_{hold_tag}_{year}_{period_tag}.tif",
        "complete": output_dir
        / "05_complete_FTC"
        / f"complete_FTC_{hold_tag}_{year}_{period_tag}.tif",
    }


@dataclass
class FTCResults:
    dftc_days: np.ndarray
    freeze_transitions: np.ndarray
    thaw_transitions: np.ndarray
    total_transitions: np.ndarray
    complete_ftc: np.ndarray
    valid_any: np.ndarray


class FTCAccumulator:
    """Accumulate daily freeze-thaw occurrence, transitions, and complete cycles."""

    def __init__(
        self,
        latitudes: np.ndarray,
        longitudes: np.ndarray,
        target_start: pd.Timestamp,
        target_end_exclusive: pd.Timestamp,
        hold_hours: int,
    ) -> None:
        self.latitudes = np.asarray(latitudes, dtype=np.float64)
        self.longitudes = np.asarray(longitudes, dtype=np.float64)
        self.shape = (self.latitudes.size, self.longitudes.size)
        self.hold = int(hold_hours)

        self.target_day_start = day_number(target_start)
        self.target_day_end = day_number(target_end_exclusive)

        # Local solar time at hourly sampling resolution: UTC + longitude/15 hours.
        self.local_offset_ns = np.rint(
            self.longitudes / 15.0 * 3600.0 * 1_000_000_000.0
        ).astype(np.int64)

        # Confirmed states: -1 frozen, +1 thawed, 0 unconfirmed/no candidate.
        self.confirmed_state = np.zeros(self.shape, dtype=np.int8)
        self.candidate_state = np.zeros(self.shape, dtype=np.int8)
        self.candidate_run = np.zeros(self.shape, dtype=np.uint8)

        # Only freezing transitions within the target period initiate eligible cycles.
        self.open_freeze = np.zeros(self.shape, dtype=bool)

        self.freeze_count = np.zeros(self.shape, dtype=np.uint32)
        self.thaw_count = np.zeros(self.shape, dtype=np.uint32)
        self.complete_count = np.zeros(self.shape, dtype=np.uint32)

        # Local-day DFTC state.
        self.current_local_day = np.full(self.longitudes.size, -2_147_483_648, dtype=np.int32)
        self.day_has_frozen = np.zeros(self.shape, dtype=bool)
        self.day_has_thawed = np.zeros(self.shape, dtype=bool)
        self.day_valid_hours = np.zeros(self.shape, dtype=np.uint8)
        self.dftc_count = np.zeros(self.shape, dtype=np.uint32)

        self.valid_any = np.zeros(self.shape, dtype=bool)
        self.valid_hours = np.zeros(self.shape, dtype=np.uint32)
        self.valid_days = np.zeros(self.shape, dtype=np.uint32)
        self.previous_time_ns: int | None = None

    def _local_day_keys(self, utc_ns: int) -> np.ndarray:
        return ((utc_ns + self.local_offset_ns) // DAY_NS).astype(np.int32)

    def _event_allowed_columns(self, event_onset_ns: int) -> np.ndarray:
        keys = self._local_day_keys(event_onset_ns)
        return (keys >= self.target_day_start) & (keys < self.target_day_end)

    def _finalize_day_columns(self, columns: np.ndarray) -> None:
        if columns.size == 0:
            return
        old_days = self.current_local_day[columns]
        in_target = (old_days >= self.target_day_start) & (old_days < self.target_day_end)
        if np.any(in_target):
            good = (
                self.day_has_frozen[:, columns]
                & self.day_has_thawed[:, columns]
                & (self.day_valid_hours[:, columns] >= MIN_VALID_HOURS_PER_LOCAL_DAY)
                & in_target[np.newaxis, :]
            )
            self.dftc_count[:, columns] += good.astype(np.uint32)
            self.valid_days[:, columns] += (
                (self.day_valid_hours[:, columns] == 24) & in_target[np.newaxis, :]
            ).astype(np.uint32)

        self.day_has_frozen[:, columns] = False
        self.day_has_thawed[:, columns] = False
        self.day_valid_hours[:, columns] = 0

    def _update_daily_metric(self, utc_ns: int, valid: np.ndarray, sign: np.ndarray) -> None:
        local_days = self._local_day_keys(utc_ns)
        uninitialized = self.current_local_day == -2_147_483_648
        if np.any(uninitialized):
            self.current_local_day[uninitialized] = local_days[uninitialized]

        changed = local_days != self.current_local_day
        if np.any(changed):
            columns = np.flatnonzero(changed)
            self._finalize_day_columns(columns)
            self.current_local_day[columns] = local_days[columns]

        self.day_has_frozen |= valid & (sign == -1)
        self.day_has_thawed |= valid & (sign == 1)
        # Strictly increasing hourly timestamps; at most 24 samples per local day.
        np.add(
            self.day_valid_hours,
            valid.astype(np.uint8),
            out=self.day_valid_hours,
            casting="unsafe",
        )

    def _reset_transition_state_for_gap(self) -> None:
        self.confirmed_state.fill(0)
        self.candidate_state.fill(0)
        self.candidate_run.fill(0)
        self.open_freeze.fill(False)

    def process(self, timestamp: pd.Timestamp, temperature_k: np.ndarray) -> None:
        if temperature_k.shape != self.shape:
            raise ValueError(f"Hourly spatial shape mismatch: {temperature_k.shape} != {self.shape}")

        utc_ns = timestamp_ns(timestamp)
        if self.previous_time_ns is not None:
            delta = utc_ns - self.previous_time_ns
            if delta <= 0:
                raise RuntimeError(f"Input timestamps are not strictly increasing: {timestamp}")
            if delta != HOUR_NS:
                # Do not pair transitions or cycles across missing intervals.
                self._reset_transition_state_for_gap()
        self.previous_time_ns = utc_ns

        valid = np.isfinite(temperature_k)
        in_period = self._event_allowed_columns(utc_ns)[None, :]
        self.valid_any |= valid & in_period
        self.valid_hours += (valid & in_period).astype(np.uint32)

        sign = np.zeros(self.shape, dtype=np.int8)
        sign[valid & (temperature_k < ZERO_C_K)] = -1
        sign[valid & (temperature_k > ZERO_C_K)] = 1

        self._update_daily_metric(utc_ns, valid, sign)

        # NaN resets confirmed states and incomplete cycles; exact zero resets candidates only.
        invalid = ~valid
        if np.any(invalid):
            self.confirmed_state[invalid] = 0
            self.candidate_state[invalid] = 0
            self.candidate_run[invalid] = 0
            self.open_freeze[invalid] = False

        neutral = valid & (sign == 0)
        if np.any(neutral):
            self.candidate_state[neutral] = 0
            self.candidate_run[neutral] = 0

        nonzero = valid & (sign != 0)
        same_candidate = nonzero & (self.candidate_state == sign)
        new_candidate = nonzero & ~same_candidate

        if np.any(same_candidate):
            incremented = self.candidate_run[same_candidate].astype(np.uint16) + 1
            self.candidate_run[same_candidate] = np.minimum(incremented, self.hold).astype(
                np.uint8
            )
        if np.any(new_candidate):
            self.candidate_state[new_candidate] = sign[new_candidate]
            self.candidate_run[new_candidate] = 1

        ready = nonzero & (self.candidate_run >= self.hold)
        changed_state = (
            ready
            & (self.confirmed_state != 0)
            & (self.candidate_state != self.confirmed_state)
        )

        if np.any(changed_state):
            # Assign events to the first hour of the H-hour confirmation sequence.
            onset_ns = utc_ns - (self.hold - 1) * HOUR_NS
            allowed_columns = self._event_allowed_columns(onset_ns)
            allowed = allowed_columns[np.newaxis, :]

            freeze_event = changed_state & (self.candidate_state == -1)
            thaw_event = changed_state & (self.candidate_state == 1)
            freeze_in_period = freeze_event & allowed
            thaw_in_period = thaw_event & allowed

            self.freeze_count += freeze_in_period.astype(np.uint32)
            self.thaw_count += thaw_in_period.astype(np.uint32)

            # Close cycles initiated by earlier eligible freezing events, then clear pending state.
            completed = thaw_in_period & self.open_freeze
            self.complete_count += completed.astype(np.uint32)
            self.open_freeze[thaw_event] = False

            # Freezing outside the target period cannot initiate an eligible complete cycle.
            self.open_freeze[freeze_event] = False
            self.open_freeze[freeze_in_period] = True

        # Initial states establish a baseline; count subsequent changes before updating state.
        state_can_update = ready & (self.candidate_state != self.confirmed_state)
        self.confirmed_state[state_can_update] = self.candidate_state[state_can_update]

    def finish(self) -> FTCResults:
        initialized_columns = np.flatnonzero(self.current_local_day != -2_147_483_648)
        self._finalize_day_columns(initialized_columns)
        total = self.freeze_count + self.thaw_count
        return FTCResults(
            dftc_days=self.dftc_count,
            freeze_transitions=self.freeze_count,
            thaw_transitions=self.thaw_count,
            total_transitions=total,
            complete_ftc=self.complete_count,
            valid_any=self.valid_any,
        )


# =============================================================================
# GeoTIFF output
# =============================================================================

def write_count_geotiff(
    path: Path,
    array: np.ndarray,
    valid_mask: np.ndarray,
    latitudes: np.ndarray,
    longitudes: np.ndarray,
    tags: Dict[str, str],
) -> None:
    # Defer rasterio import to allow independent algorithm self-tests.
    import rasterio
    from rasterio.transform import from_origin

    path.parent.mkdir(parents=True, exist_ok=True)
    if array.shape != valid_mask.shape:
        raise ValueError("Output array and validity mask shapes differ.")
    if int(np.nanmax(array)) >= NODATA_UINT16:
        raise OverflowError(f"{path.name} exceeds the uint16 range.")

    lats = np.asarray(latitudes, dtype=np.float64)
    lons = np.asarray(longitudes, dtype=np.float64)
    xres = float(np.median(np.diff(lons)))
    yres = float(abs(np.median(np.diff(lats))))
    if not np.allclose(np.diff(lons), xres, atol=1e-7):
        raise RuntimeError("Irregular longitude grid; cannot write GeoTIFF directly.")
    if not np.allclose(np.abs(np.diff(lats)), yres, atol=1e-7):
        raise RuntimeError("Irregular latitude grid; cannot write GeoTIFF directly.")

    # The first GeoTIFF row corresponds to the northernmost latitude.
    if lats[0] < lats[-1]:
        data_to_write = np.flipud(array)
        mask_to_write = np.flipud(valid_mask)
    else:
        data_to_write = array
        mask_to_write = valid_mask

    out = np.full(array.shape, NODATA_UINT16, dtype=np.uint16)
    converted = data_to_write.astype(np.uint16)
    out[mask_to_write] = converted[mask_to_write]

    west = float(lons[0] - xres / 2.0)
    north = float(np.max(lats) + yres / 2.0)
    transform = from_origin(west, north, xres, yres)

    temp_path = path.with_suffix(path.suffix + ".part")
    profile = {
        "driver": "GTiff",
        "height": out.shape[0],
        "width": out.shape[1],
        "count": 1,
        "dtype": "uint16",
        "crs": "EPSG:4326",
        "transform": transform,
        "nodata": NODATA_UINT16,
        "compress": "DEFLATE",
        "predictor": 2,
        "zlevel": 6,
        "tiled": True,
        "blockxsize": 256,
        "blockysize": 256,
        "BIGTIFF": "IF_SAFER",
    }
    with rasterio.open(temp_path, "w", **profile) as dst:
        dst.write(out, 1)
        dst.update_tags(**tags)

    os.replace(temp_path, path)


# -----------------------------------------------------------------------------
# Index internal valid times; filenames are used only for discovery.
# -----------------------------------------------------------------------------
def normalize_longitudes(values):
    """Round the 0.1-degree grid before mapping +180 to -180."""
    rounded = np.round(np.asarray(values, dtype=np.float64), 6)
    normalized = np.round((rounded + 180.0) % 360.0 - 180.0, 6)
    return np.where(normalized >= 180.0, -180.0, normalized)


def hour_key(value):
    ts = pd.Timestamp(value)
    if ts.value % HOUR_NS:
        raise ValueError(f"Timestamp is not an exact UTC hour: {ts}")
    return int(ts.value // HOUR_NS)


def signature(path):
    st = Path(path).stat()
    return [st.st_size, st.st_mtime_ns]


def open_nc_local(path):
    """Open Unicode paths in Python and pass binary handles to the backend."""
    import xarray as xr
    raw = open(path, 'rb')
    try:
        magic = raw.read(8)
        raw.seek(0)
        if magic[:4] == b'CDF\x01' or magic[:4] == b'CDF\x02':
            ds = xr.open_dataset(raw, engine='scipy', decode_times=True,
                                 mask_and_scale=True, cache=False)
        elif magic == b'\x89HDF\r\n\x1a\n':
            ds = xr.open_dataset(raw, engine='h5netcdf', decode_times=True,
                                 mask_and_scale=True, cache=False)
        else:
            raise ValueError(f'Unsupported NetCDF3/NetCDF4 file; header={magic.hex()}')
        backend_close = ds._close
        def close_both():
            try:
                if backend_close:
                    backend_close()
            finally:
                raw.close()
        ds.set_close(close_both)
        return ds
    except BaseException:
        raw.close()
        raise


def nc_variable(ds):
    for name in ('stl1', 'soil_temperature_level_1', 'sktl1', 'var139'):
        if name in ds:
            return ds[name]
    names = [n for n, v in ds.data_vars.items()
             if str(v.attrs.get('GRIB_paramId', '')) == '139']
    if len(names) == 1:
        return ds[names[0]]
    raise ValueError(f"stl1 variable not identified; available variables: {list(ds.data_vars)}")


def nc_time(da):
    # Require a one-dimensional time or valid_time coordinate.
    for name in ('valid_time', 'time'):
        if name in da.coords and da[name].ndim == 1:
            times = np.asarray(da[name].values)
            return da[name].dims[0], [hour_key(t) for t in times]
    raise ValueError('A one-dimensional hourly time/valid_time coordinate is required.')


def grib_hour(gid):
    import eccodes as ec
    date = int(ec.codes_get(gid, 'validityDate'))
    clock = int(ec.codes_get(gid, 'validityTime'))
    return hour_key(pd.to_datetime(f'{date:08d}{clock:04d}', format='%Y%m%d%H%M'))


def scan_one(path, kind):
    rows = []
    warning = ''
    if kind == 'nc':
        import xarray as xr
        with open_nc_local(path) as ds:
            _, times = nc_time(nc_variable(ds))
            rows = [[h, i] for i, h in enumerate(times)]
    else:
        import eccodes as ec
        with open(path, 'rb') as stream:
            while True:
                gid = None
                try:
                    gid = ec.codes_grib_new_from_file(stream)
                    if gid is None:
                        break
                    if int(ec.codes_get(gid, 'paramId')) != 139:
                        continue
                    offset = int(ec.codes_get_message_offset(gid))
                    length = int(ec.codes_get_message_size(gid))
                    # Verify the message terminator in addition to ecCodes decoding.
                    with open(path, 'rb') as check:
                        check.seek(offset + length - 4)
                        if check.read(4) != b'7777':
                            raise ValueError('GRIB message lacks the 7777 terminator.')
                    rows.append([grib_hour(gid), offset])
                except Exception as exc:
                    warning = f'Scan stopped at damaged/truncated data; retained {len(rows)} records: {exc}'
                    break
                finally:
                    if gid is not None:
                        ec.codes_release(gid)
    return rows, warning


def build_index(out, logger):
    for root in (GRIB_ROOT, NC_ROOT):
        if not root.exists():
            logger.warning('Directory does not exist: %s', root)
    discovered = []
    for root, kind in ((GRIB_ROOT, 'grib'), (NC_ROOT, 'nc')):
        if root.exists():
            for p in root.rglob('*'):
                name = p.name.lower()
                if p.is_file() and ((kind == 'grib' and name.endswith(
                        ('.grib', '.grb', '.grib2', '.grib.part', '.grb.part', '.grib2.part')))
                        or (kind == 'nc' and name.endswith('.nc'))):
                    discovered.append((str(p.resolve()), kind))
    if not discovered:
        raise FileNotFoundError('No files with supported extensions found in either input directory.')
    # Priority: GRIB, GRIB.part, NC; sort equal-priority paths and deduplicate hours.
    discovered.sort(key=lambda x: (2 if x[1] == 'nc' else int(x[0].endswith('.part')), x[0]))
    cache_path = out / 'file_index_cache.json'
    try:
        candidates = [cache_path, OUTPUT_DIR/'v2_fixed'/'file_index_cache.json', OUTPUT_DIR/'file_index_cache.json']
        source_cache = next((p for p in candidates if p.exists()), cache_path)
        old = json.loads(source_cache.read_text(encoding='utf-8'))
    except (OSError, ValueError):
        old = {}
    current, timeline, report = {}, defaultdict(list), []
    for i, (path, kind) in enumerate(discovered, 1):
        sig = signature(path)
        item = old.get(path)
        index_version = GRIB_INDEX_VERSION if kind == 'grib' else NC_INDEX_VERSION
        if not item or item.get('signature') != sig or item.get('version') != index_version or (kind == 'nc' and not item.get('rows')):
            try:
                rows, warning = scan_one(path, kind)
            except Exception as exc:
                rows, warning = [], f'{type(exc).__name__}: {exc}'
            if signature(path) != sig:
                raise RuntimeError(f'File changed during indexing; pause downloads and rerun: {path}')
            item = dict(signature=sig, rows=rows, warning=warning, kind=kind, version=index_version)
        current[path] = item
        if item['warning']:
            logger.warning('%s | %s', path, item['warning'])
        report.append(dict(path=path, kind=kind, indexed_hours=len(item['rows']),
                           warning=item['warning']))
        for h, offset in item['rows']:
            timeline[h].append((path, kind, offset))
        if i % 50 == 0 or i == len(discovered):
            logger.info('Indexed/cached %d/%d files; %d UTC hours available', i, len(discovered), len(timeline))
    temp = cache_path.with_suffix('.tmp')
    temp.write_text(json.dumps(current, ensure_ascii=False), encoding='utf-8')
    os.replace(temp, cache_path)
    pd.DataFrame(report).to_csv(out / 'File_scan_inventory.csv', index=False, encoding='utf-8-sig')
    nc_items = [v for v in current.values() if v['kind'] == 'nc']
    if nc_items and not any(v['rows'] for v in nc_items):
        raise RuntimeError('No NetCDF files could be read; processing stopped. See the file inventory for errors.')
    digest = hashlib.sha256(json.dumps(current, sort_keys=True).encode()).hexdigest()
    return timeline, current, digest


class HourReader:
    """Keep at most two file handles; decode single hours rather than seven-day arrays."""
    def __init__(self, inventory):
        self.inventory = inventory
        self.handles = OrderedDict()
        self.grid = (np.round(np.arange(GRID_LAT_MAX, GRID_LAT_MIN-GRID_STEP/2, -GRID_STEP), 6),
                     np.round(np.arange(GRID_LON_MIN, GRID_LON_MAX+GRID_STEP/2, GRID_STEP), 6))
        if QUICK_TEST:
            a, b = self.grid
            self.grid = (a[(a>=30)&(a<=35)], b[(b>=100)&(b<=105)])
        self.layout_records = {}
        self.mapping = None
        self.mapping_key = None
        self.failed_files = set()

    def close(self):
        for obj in self.handles.values():
            obj.close()
        self.handles.clear()

    def handle(self, path, kind):
        if signature(path) != self.inventory[path]['signature']:
            raise RuntimeError(f'INPUT_CHANGED: Input changed; pause downloads and rerun: {path}')
        if path in self.handles:
            self.handles.move_to_end(path)
            return self.handles[path]
        if len(self.handles) >= 2:
            _, obj = self.handles.popitem(last=False)
            obj.close()
        if kind == 'grib':
            obj = open(path, 'rb')
        else:
            import xarray as xr
            obj = open_nc_local(path)
        self.handles[path] = obj
        return obj

    def align(self, data, lat, lon):
        """Align same-resolution regional subsets exactly; leave exterior cells as NaN."""
        lat = np.round(np.asarray(lat, dtype=float), 6)
        lon = normalize_longitudes(lon)
        if data.shape != (len(lat), len(lon)):
            raise ValueError('GRID_UNSUPPORTED: Temperature shape does not match coordinate lengths.')
        if not (np.isfinite(lat).all() and np.isfinite(lon).all()):
            raise ValueError('GRID_UNSUPPORTED: Nonfinite longitude or latitude coordinates.')
        # Reject duplicate coordinates before assignment.
        if len(np.unique(lat)) != len(lat) or len(np.unique(lon)) != len(lon):
            raise ValueError('GRID_UNSUPPORTED: Duplicate coordinates (possibly 0/360 or -180/180).')
        for name, coord in (('latitude',lat),('longitude',lon)):
            steps = np.diff(np.sort(coord))
            # Allow one longitude-seam gap; remaining spacing must be 0.1 degrees.
            regular = np.isclose(steps, GRID_STEP, atol=1e-5)
            seam = (name == 'longitude' and np.count_nonzero(~regular) <= 1 and regular.any()
                    and np.all(steps[~regular] > GRID_STEP))
            if len(steps) and not (regular.all() or seam):
                raise ValueError(f'GRID_UNSUPPORTED: {name} spacing range {steps.min():.6f}—{steps.max():.6f}，'
                                 f'target spacing {GRID_STEP}; no interpolation applied')
        target_lat, target_lon = self.grid
        rows = np.rint((target_lat[0]-lat)/GRID_STEP).astype(int)
        cols = np.rint((lon-target_lon[0])/GRID_STEP).astype(int)
        iy = np.flatnonzero((rows>=0)&(rows<len(target_lat)))
        ix = np.flatnonzero((cols>=0)&(cols<len(target_lon)))
        if not len(iy) or not len(ix):
            raise ValueError('GRID_UNSUPPORTED: Input and target domains do not overlap.')
        if not np.allclose(lat[iy],target_lat[rows[iy]],atol=1e-5,rtol=0) or not np.allclose(
                lon[ix],target_lon[cols[ix]],atol=1e-5,rtol=0):
            raise ValueError('GRID_UNSUPPORTED: Shifted grid origin prevents exact alignment; no interpolation applied')
        path = self.active_path
        if path not in self.layout_records:
            rec = dict(path=path, rows=len(lat), cols=len(lon),
                       lat_min=float(lat.min()),lat_max=float(lat.max()),
                       lon_min=float(lon.min()),lon_max=float(lon.max()),
                       target_rows=len(target_lat),target_cols=len(target_lon),
                       matched_cells=len(iy)*len(ix),
                       full_target=len(iy)*len(ix)==len(target_lat)*len(target_lon))
            self.layout_records[path] = rec
            if not rec['full_target']:
                logging.getLogger('FTC_LOCAL').warning('Regional grid aligned by coordinates (alternatives may fill gaps): %s | %s',path,rec)
        out = np.full((len(target_lat),len(target_lon)),np.nan,dtype=np.float32)
        out[np.ix_(rows[iy],cols[ix])] = data[np.ix_(iy,ix)]
        return out

    def read(self, source, expected_h):
        path, kind, pos = source
        self.active_path = path
        obj = self.handle(path, kind)
        if kind == 'grib':
            import eccodes as ec
            obj.seek(pos)
            gid = ec.codes_grib_new_from_file(obj)
            if gid is None:
                raise ValueError('No complete GRIB message at this offset.')
            try:
                if grib_hour(gid) != expected_h:
                    raise ValueError('Index and actual valid time differ.')
                if str(ec.codes_get(gid, 'gridType')) != 'regular_ll':
                    raise ValueError('Only regular_ll latitude-longitude grids are supported.')
                key = (path, str(ec.codes_get(gid, 'md5GridSection')))
                if key != self.mapping_key:
                    la = np.round(ec.codes_get_array(gid, 'latitudes'), 6)
                    lo = normalize_longitudes(ec.codes_get_array(gid, 'longitudes'))
                    lat, lon = np.unique(la), np.unique(lo)
                    flat = np.searchsorted(lat, la) * len(lon) + np.searchsorted(lon, lo)
                    if len(flat) != len(lat) * len(lon) or len(np.unique(flat)) != len(flat):
                        raise ValueError('GRIB grid is not a complete rectangle.')
                    self.mapping = (lat, lon, np.argsort(flat))
                    self.mapping_key = key
                lat, lon, order = self.mapping
                values = np.asarray(ec.codes_get_values(gid), dtype=np.float32)
                if int(ec.codes_get(gid, 'bitmapPresent')):
                    bitmap = ec.codes_get_array(gid, 'bitmap')
                    values[bitmap == 0] = np.nan
                unit = str(ec.codes_get(gid, 'units'))
                data = values[order].reshape(len(lat), len(lon))
            finally:
                ec.codes_release(gid)
        else:
            da = nc_variable(obj)
            dim, times = nc_time(da)
            if times[pos] != expected_h:
                raise ValueError('NCIndex and actual valid time differ.')
            da = da.isel({dim: pos})
            # Merge nonmissing ERA5 expver values, prioritizing version 1 over version 5.
            if 'expver' in da.dims:
                versions = sorted(da.expver.values.tolist())
                combined = da.sel(expver=versions[0], drop=True)
                for version in versions[1:]:
                    combined = combined.fillna(da.sel(expver=version, drop=True))
                da = combined
            latname = next((n for n in ('latitude', 'lat') if n in da.coords), None)
            lonname = next((n for n in ('longitude', 'lon') if n in da.coords), None)
            if latname is None or lonname is None:
                raise ValueError('NetCDF lacks longitude/latitude coordinates.')
            if da[latname].ndim != 1 or da[lonname].ndim != 1:
                raise ValueError('NetCDF requires one-dimensional regular longitude/latitude coordinates.')
            ydim, xdim = da[latname].dims[0], da[lonname].dims[0]
            for d in list(da.dims):
                if d not in (ydim, xdim):
                    if da.sizes[d] != 1:
                        raise ValueError(f'NetCDF has an unsupported nonsingleton dimension: {d}')
                    da = da.isel({d: 0})
            unit = str(da.attrs.get('units', ''))
            data = np.asarray(da.transpose(ydim, xdim).values, dtype=np.float32)
            lat, lon = da[latname].values, da[lonname].values
        unit = unit.lower().strip().replace(' ', '').replace('°', 'deg')
        if unit in ('c', 'degc', 'celsius', 'degree_celsius', 'degrees_celsius'):
            data = data + ZERO_C_K
        elif unit not in ('k', 'kelvin', 'degree_kelvin', 'degrees_kelvin'):
            raise ValueError(f'Cannot determine temperature units: {unit!r}; units are not inferred from value magnitude')
        # Exclude unmasked sentinel values using physically impossible absolute temperatures.
        data = np.where((data > 0) & (data < 1000), data, np.nan)
        return self.align(data, lat, lon)


def self_test():
    np.testing.assert_allclose(normalize_longitudes(
        [179.9999999, 180.0, -180.0, 359.9, 0.0]), [-180., -180., -180., -.1, 0.])

    def run(values, lon=0, start='2000-01-01', end='2000-01-03', origin=None, gaps=None):
        a = FTCAccumulator(np.array([0.]), np.array([lon]), pd.Timestamp(start), pd.Timestamp(end), 3)
        t = pd.Timestamp(origin or start)
        for i, v in enumerate(values):
            if gaps and i in gaps:
                continue
            a.process(t + pd.Timedelta(hours=i), np.array([[v + float(ZERO_C_K)]], dtype=np.float32))
        return a, a.finish()
    _, r = run([1]*3 + [-1]*3 + [1]*3)
    assert (r.freeze_transitions.item(), r.thaw_transitions.item(), r.complete_ftc.item()) == (1,1,1)
    _, r = run([-1]*3 + [1]*3)
    assert (r.freeze_transitions.item(), r.thaw_transitions.item(), r.complete_ftc.item()) == (0,1,0)
    for middle in ([-1]*2, [-1,0,-1], [-1,float('nan'),-1]):
        _, r = run([1]*3 + middle + [1]*3)
        assert r.total_transitions.item() == 0
    _, r = run([1]*3 + [-1]*3 + [float('nan')] + [1]*3)
    assert r.complete_ftc.item() == 0
    _, r = run([1]*3 + [-1]*3 + [1]*4, gaps={6})
    assert r.complete_ftc.item() == 0
    a,r = run([-1]*12+[1]*12)
    assert r.dftc_days.item() == 1 and a.valid_days.item() == 1
    a,r = run([-1]*12+[1]*11)
    assert r.dftc_days.item() == 0 and a.valid_days.item() == 0
    a,r = run([-1]*12+[1]*12, lon=90, origin='1999-12-31 18:00', end='2000-01-02')
    assert r.dftc_days.item() == 1 and a.valid_hours.item() == 24
    # Pre-period freezing cannot initiate an eligible cycle even if confirmed in-period.
    _,r = run([1]*3+[-1]*3+[1]*3, origin='1999-12-31 20:00', end='2000-01-02')
    assert r.freeze_transitions.item() == 0 and r.complete_ftc.item() == 0
    # Assign freezing starting before period end to that period even if confirmed later.
    _,r = run([1]*3+[-1]*3, origin='2000-01-01 20:00', end='2000-01-02')
    assert r.freeze_transitions.item() == 1
    print('Self-tests passed: three-hour confirmation, brief fluctuations, zero, missing data, initial freezing, solar days, and period boundaries.')


def main():
    if '--self-test' in sys.argv:
        self_test()
        return
    if not (1 <= MONTHS[0] <= MONTHS[-1] <= 12) or MONTHS != tuple(range(MONTHS[0], MONTHS[-1]+1)):
        raise ValueError('MONTHS must be consecutive ascending months between 1 and 12.')
    if START_YEAR > END_YEAR or not 1 <= HOLD_HOURS <= 24:
        raise ValueError('Invalid year range or HOLD_HOURS.')
    import eccodes
    import xarray
    import rasterio
    import h5netcdf
    import h5py
    import scipy
    out = OUTPUT_DIR / ('v3_quick_test' if QUICK_TEST else 'v3_fixed')
    out.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(level=logging.INFO, format='%(asctime)s | %(levelname)s | %(message)s',
                        handlers=[logging.StreamHandler(sys.stdout),
                                  logging.FileHandler(out/'FTC_run_log.txt', encoding='utf-8')])
    logger = logging.getLogger('FTC_LOCAL')
    logger.info('Counts cover valid observations; assess completeness using 06_quality. Months: %s', MONTHS)
    timeline, inventory, digest = build_index(out, logger)
    configuration = dict(version=SCRIPT_VERSION, inventory=digest, months=MONTHS,
                         hold=HOLD_HOURS, daily_hours=MIN_VALID_HOURS_PER_LOCAL_DAY,
                         quick_test=QUICK_TEST,grid=[GRID_STEP, GRID_LAT_MIN, GRID_LAT_MAX, GRID_LON_MIN, GRID_LON_MAX])
    fingerprint = hashlib.sha256(json.dumps(configuration, sort_keys=True).encode()).hexdigest()
    years = [START_YEAR] if QUICK_TEST else range(START_YEAR, END_YEAR+1)
    for year in years:
        marker = out / f'{year}_completed.json'
        if not OVERWRITE_EXISTING and marker.exists():
            try:
                done = json.loads(marker.read_text(encoding='utf-8'))
                if done['fingerprint'] == fingerprint and all(Path(p).exists() for p in done['outputs']):
                    logger.info('%d results match the input version; skipping', year)
                    continue
            except (ValueError, KeyError):
                pass
        # Remove completion markers before recalculation to prevent stale-result reuse.
        marker.unlink(missing_ok=True)
        start, end = period_bounds(year)
        h0, h1 = hour_key(start-pd.Timedelta(hours=24)), hour_key(end+pd.Timedelta(hours=24))
        reader = HourReader(inventory)
        acc = None
        source_counts = defaultdict(int)
        reportpath = out / f'{year}_Hourly_sources_and_missing_data.csv'
        try:
            with open(reportpath, 'w', newline='', encoding='utf-8-sig') as report:
                writer = csv.writer(report)
                writer.writerow(['UTC_time','status','sources_used','finite_pixels','filled_from_backups','errors'])
                for h in range(h0, h1):
                    stamp = pd.Timestamp(h*HOUR_NS)
                    merged = None
                    used, errors = [], []
                    filled = 0
                    for source in timeline.get(h, []):
                        try:
                            data = reader.read(source, h)
                        except Exception as exc:
                            if 'INPUT_CHANGED:' in str(exc):
                                raise
                            errors.append(f'{source[0]}: {type(exc).__name__}: {exc}')
                            if source[0] not in reader.failed_files:
                                logger.warning('Read failed; trying an alternative: %s | %s', source[0], exc)
                                reader.failed_files.add(source[0])
                            continue
                        if merged is None:
                            merged = data
                            used.append(source[0])
                        else:
                            fill = ~np.isfinite(merged) & np.isfinite(data)
                            nfill = int(np.count_nonzero(fill))
                            if nfill:
                                merged[fill] = data[fill]
                                filled += nfill
                                used.append(source[0])
                        if np.isfinite(merged).all():
                            break
                    if merged is not None:
                        if acc is None:
                            acc = FTCAccumulator(*reader.grid, start, end, HOLD_HOURS)
                        acc.process(stamp, merged)
                        finite = int(np.count_nonzero(np.isfinite(merged)))
                        status = 'read' if finite else 'all_nan'
                    else:
                        finite, status = 0, 'missing'
                        # Skip missing hours; the next valid timestamp triggers gap detection and reset.
                    source_counts[status] += 1
                    writer.writerow([str(stamp),status,' | '.join(used),finite,filled,' | '.join(errors)])
                    if (h-h0) % 24 == 0:
                        logger.info('%d | %s | %s | valid pixels %d', year, stamp, status, finite)
            if acc is None:
                logger.error('%d has no decodable data; saved the missing-data report only, without results or a completion marker', year)
                continue
            results = acc.finish()
            paths = expected_output_paths(out, year)
            arrays = [results.dftc_days, results.freeze_transitions, results.thaw_transitions,
                      results.total_transitions, results.complete_ftc]
            tags = dict(source='local ERA5-Land stl1; GRIB then daily NetCDF fallback',
                        year=str(year), months=str(MONTHS), hold_hours=str(HOLD_HOURS),
                        event_time='first hour of confirmed run',
                        local_day='UTC + longitude/15 hours; mean solar day',
                        coverage='observed counts only; see 06_quality; gaps reset states',
                        complete_cycle='freeze transition then thaw transition; both onsets in target period',
                        daily_requirement='24 finite hourly samples', fingerprint=fingerprint)
            for i, (path, array) in enumerate(zip(paths.values(), arrays)):
                mask = acc.valid_days > 0 if i == 0 else acc.valid_any
                write_count_geotiff(path, array, mask, *reader.grid, tags)
            expected_hours = int((end-start).total_seconds()//3600)
            quality = {'valid_hours': acc.valid_hours, 'valid_local_days': acc.valid_days,
                       'missing_hours': expected_hours-acc.valid_hours}
            output_paths = [str(p) for p in paths.values()]
            for name, array in quality.items():
                p = out/'06_quality'/f'{name}_{year}_M{MONTHS[0]:02d}-M{MONTHS[-1]:02d}.tif'
                # Retain NoData for unobserved pixels to avoid classifying oceans as missing land data.
                write_count_geotiff(p, array, acc.valid_any, *reader.grid, tags)
                output_paths.append(str(p))
            land = acc.valid_any
            summary = dict(year=year, months=MONTHS, expected_hours=expected_hours,
                           expected_local_days=expected_hours//24,
                           observed_pixels=int(land.sum()),
                           fully_covered_pixels=int(((acc.valid_hours == expected_hours)&land).sum()),
                           hourly_read_status_including_padding=dict(source_counts),
                           warning='All five metrics count valid observations; missing data may cause underestimation. Coverage is recorded in 06_quality.')
            summary_path = out/f'{year}_Coverage_summary.json'
            summary_path.write_text(json.dumps(summary,ensure_ascii=False,indent=2), encoding='utf-8')
            output_paths += [str(reportpath), str(summary_path)]
            # Recheck annual input signatures before marking completion.
            used_paths = {s[0] for h in range(h0,h1) for s in timeline.get(h,[])}
            for p in used_paths:
                if signature(p) != inventory[p]['signature']:
                    raise RuntimeError(f'Input changed during processing; rerun: {p}')
            marker.write_text(json.dumps(dict(fingerprint=fingerprint,outputs=output_paths),
                                         ensure_ascii=False,indent=2), encoding='utf-8')
            logger.info('%d complete: %d observed pixels; %d pixels covering all target hours',
                        year, summary['observed_pixels'], summary['fully_covered_pixels'])
        finally:
            reader.close()
            if reader.layout_records:
                pd.DataFrame(list(reader.layout_records.values())).to_csv(
                    out/f'{year}_Grid_diagnostics.csv', index=False, encoding='utf-8-sig')
            del acc
            gc.collect()
    logger.info('Processing complete. Output directory: %s', out)


if __name__ == '__main__':
    main()
