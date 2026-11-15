"""Soil temperature: Northern Hemisphere snow-depth analysis.

1. Align inputs and apply the common snow and glacier mask.
2. Estimate pixelwise partial correlations with AR(1) adjustment.
3. Fit the equal-amplitude panel model and report predictor VIFs.
"""
from __future__ import annotations
import os
for _key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ.setdefault(_key, '1')
import argparse
import calendar
import gc
import json
import math
import re
from datetime import datetime
from pathlib import Path
from contextlib import contextmanager, ExitStack
import numpy as np
import pandas as pd
from scipy import stats
import rasterio
from rasterio.transform import Affine, from_bounds
from rasterio.warp import reproject, transform as transform_points
from rasterio.enums import Resampling
from rasterio.windows import Window
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Circle

# 1. Indicator and path settings
INDICATOR = 'ST'
TARGET_VAR = 'stl1'
CONTROLS = ['t2m','ssr','v10']
RESPONSE_LABEL = 'Soil temperature response (°C)'

WINTER = {'previous':[11,12], 'current':[1,2,3]}
PERIODS = {
    'Freeze_NovMar': {'previous':[11,12], 'current':[1,2,3]},
    'FreezeThaw_AprJun': {'previous':[], 'current':[4,5,6]},
    'Growing_JulOct': {'previous':[], 'current':[7,8,9,10]},
    'Annual_NovOct': {'previous':[11,12], 'current':list(range(1,11))},
}

PARTIAL_YEARS = list(range(2001,2026))
PANEL_YEARS = list(range(2001,2023))
REQUIRE_COMPLETE_PARTIAL = True
MIN_PARTIAL_YEARS = len(PARTIAL_YEARS)
MIN_PANEL_OBSERVATIONS = 2
ADJUSTED_RESPONSE_DEFINITION = 'model_residual_plus_original_snow_contribution'

PROJECT_ROOT = Path(__file__).resolve().parent
INPUT_ROOT = Path(os.environ.get('SNOWPACK_INPUT_ROOT', PROJECT_ROOT / 'data')).expanduser().resolve()
OUTPUT_ROOT = Path(os.environ.get('SNOWPACK_OUTPUT_ROOT', PROJECT_ROOT / 'results')).expanduser().resolve()

ERA_ROOT = INPUT_ROOT / 'era5_land'
FREQUENCY_CANDIDATES = [
    INPUT_ROOT / 'snow_cover_frequency' / name
    for name in ('snowc_frequency_2000_2025.tif', 'snowc_frequency_2000_2025_common_grid.tif',
                 'snowc_frequency.tif')
]
WORLD_SHP = INPUT_ROOT / 'boundaries' / 'global_all_country.shp'
WORLD_GDF = None
SNOW_STACK_DIR = INPUT_ROOT / 'snow_off' / 'annual_stack'
FTC_ROOT = INPUT_ROOT / 'freeze_thaw_cycles'
SDE_ROOT = INPUT_ROOT / 'daily_snow_depth'
LEGACY_OUT_ROOT = PROJECT_ROOT / 'previous_results' / INDICATOR
OUT_ROOT = OUTPUT_ROOT / INDICATOR
GLACIER_MASK_PATH = INPUT_ROOT / 'glacier_mask' / 'cicecap_v015_1279_4_regMIR.nc'
GLACIER_MASK_URL = ('https://confluence.ecmwf.int/download/attachments/140385202/'
                    'cicecap_v015_1279_4_regMIR.nc?api=v2&modificationDate=1711562368310&version=1')
GLACIER_FRACTION_THRESHOLD = 0.50


ANALYSIS_MIN_LAT = 0.0
MAP_MIN_LAT = 24.0
SNOW_FREQUENCY_THRESHOLD = 0.30
FREQUENCY_YEARS = list(range(2000, 2026))
BLOCK_ROWS = 8
PARTIAL_PIXEL_CHUNK = 1000
MIN_ADJACENT_PAIRS = 3
CAP_NEFF_AT_N = False
ALPHA = 0.05
MIN_BIN_N = 10
MIN_MODEL_N = 100
MIN_MODEL_CLUSTERS = 2
AMP_CENTERS_CM = np.array([1, 5, 10, 15, 20, 25, 30, 40, 50, 60], float)
# Bin intervals are (0, 3), left-closed interior bins, and [55, 60] cm.
AMP_EDGES_CM = np.array([0, 3, 7.5, 12.5, 17.5, 22.5, 27.5, 35, 45, 55, 60], float)
MAX_AMPLITUDE_CM = 60.0
# CR1 counts explicit regressors, excluding absorbed fixed effects.

CR1_COUNT_ABSORBED_FE = False
NODATA = -9999.0
NEGATIVE_COLOR = '#5A769E'
POSITIVE_COLOR = '#DF7E21'
FIGSIZE = (9.0, 9.0)
FIG_DPI = 600
MAP_AX_RECT = [0.10, 0.10, 0.80, 0.80]
MAP_PARALLELS = [60, 80]
MAP_MERIDIANS = [-135, -90, -45, 0, 45, 90, 135, 180]
MAP_BACKGROUND_COLOR = '#FFFFFF'
MAP_GRID_COLOR = '#D3D3D3'
MAP_COUNTRY_COLOR = MAP_GRID_COLOR
MAP_BORDER_COLOR = '#000000'
MAP_RASTER_PIXELS = 2400
POLAR_CRS = '+proj=stere +lat_0=90 +lat_ts=90 +lon_0=0 +datum=WGS84 +units=m +no_defs'
CORR_CMAP = LinearSegmentedColormap.from_list('snow_blue_white_orange',
    [NEGATIVE_COLOR, '#FFFFFF', POSITIVE_COLOR], N=256)
CORR_CMAP.set_bad((1, 1, 1, 0))
_available_fonts = {f.name for f in matplotlib.font_manager.fontManager.ttflist}
_plot_font = next((name for name in ['Arial', 'Liberation Sans', 'Noto Sans', 'DejaVu Sans']
                   if name in _available_fonts), 'sans-serif')
plt.rcParams.update({'font.family': _plot_font,
                     'axes.unicode_minus': False, 'font.size': 12})
RUN_DIR = None
SOURCE_RECORDS = {}


def log(message):
    """Log progress to the console and current run directory."""
    text = f'[{datetime.now():%H:%M:%S}] {message}'
    print(text, flush=True)
    if RUN_DIR is not None:
        with (RUN_DIR / 'run_log.txt').open('a', encoding='utf-8') as fh:
            fh.write(text + '\n')


def record_source(path):
    p = Path(path).resolve()
    if str(p) not in SOURCE_RECORDS:
        s = p.stat()
        SOURCE_RECORDS[str(p)] = {'path': str(p), 'size_bytes': s.st_size, 'mtime_ns': s.st_mtime_ns}


def monthly_path(var, year, month):
    folder = ERA_ROOT / var
    path = folder / f'{year}_{month:02d}.tif'
    if path.is_file():
        return path
    date = f'{year}{month:02d}01'
    matches = sorted(folder.glob(f'*_{var}_{date}.tif'))
    if len(matches) == 1:
        return matches[0]
    if len(matches) > 1:
        raise ValueError(f'Multiple monthly rasters for {var} ({year}-{month:02d}): {matches}')
    raise FileNotFoundError(f'Missing monthly raster: {path}')


def window_months(year, spec):
    return [(year-1, m) for m in spec['previous']] + [(year, m) for m in spec['current']]


def template_from_path(path):
    with rasterio.open(path) as src:
        if src.crs is None or src.crs.to_epsg() != 4326:
            raise ValueError(f'The template must use EPSG:4326: {path}')
        tr = src.transform
        if tr.a <= 0 or tr.e >= 0 or abs(tr.b) > 1e-10 or abs(tr.d) > 1e-10:
            raise ValueError(f'The template must be a north-up grid with increasing longitude: {path}')
        lon = tr.c + (np.arange(src.width) + .5) * tr.a
        lat = tr.f + (np.arange(src.height) + .5) * tr.e
        profile = {'driver': 'GTiff', 'width': src.width, 'height': src.height,
                   'count': 1, 'dtype': 'float32', 'crs': src.crs, 'transform': tr,
                   'nodata': NODATA, 'compress': 'lzw', 'BIGTIFF': 'IF_SAFER'}
    return {'profile': profile, 'height': len(lat), 'width': len(lon), 'lon': lon, 'lat': lat,
            'transform': tr, 'crs': profile['crs']}


def wrap_geographic(a, tr, crs):
    """Wrap a global 0-360 degree raster to the -180 to 180 degree range."""
    width = a.shape[-1]
    if crs.to_epsg() == 4326 and tr.c > -90 and tr.c + tr.a * width > 181:
        if abs(tr.a*width-360) > abs(tr.a)*1.1:
            raise ValueError('Longitude wrapping requires a global grid.')
        centers = tr.c + (np.arange(width) + .5)*tr.a
        split = int(np.searchsorted(centers, 180))
        a = np.concatenate([a[..., split:], a[..., :split]], axis=-1)
        tr = Affine(tr.a, 0, tr.c + split*tr.a - 360, 0, tr.e, tr.f)
    return a, tr


def coordinate_indices(source, target, cyclic=False):
    """Match coincident grid centers for count and daily data."""
    source = np.asarray(source, float)
    target = np.asarray(target, float)
    if cyclic:
        source = np.mod(np.round(source, 7), 360)
        target = np.mod(np.round(target, 7), 360)
    order = np.argsort(source)
    s = source[order]
    if np.any(np.diff(s) < 1e-7):
        raise ValueError('Duplicate grid coordinates. Check the longitude seam.')
    j = np.clip(np.searchsorted(s, target), 0, len(s)-1)
    jm = np.maximum(j-1, 0)
    j = np.where(np.abs(s[jm]-target) < np.abs(s[j]-target), jm, j)
    good = np.abs(s[j]-target) < 1e-5
    interior = (target >= s[0]-1e-5) & (target <= s[-1]+1e-5)
    if np.any(interior & ~good):
        raise ValueError('Count or daily snow-depth grid centers do not match the template.')
    return order[j], good


def read_aligned(path, template, band=1, categorical=False, exact=False):
    """Align a raster band using bilinear or nearest-neighbor resampling as appropriate."""
    record_source(path)
    with rasterio.open(path) as src:
        if band > src.count or src.crs is None:
            raise ValueError(f'Missing band or CRS: {path}, band {band}')
        a = src.read(band, masked=True).astype('float32').filled(np.nan)
        a[~np.isfinite(a)] = np.nan
        tr, crs = src.transform, src.crs
    if exact:
        if crs.to_epsg() != 4326 or tr.b or tr.d:
            raise ValueError(f'Exact matching requires a regular geographic grid: {path}')
        lon = tr.c + (np.arange(a.shape[1])+.5)*tr.a
        lat = tr.f + (np.arange(a.shape[0])+.5)*tr.e
        ri, rv = coordinate_indices(lat, template['lat'])
        ci, cv = coordinate_indices(lon, template['lon'], True)
        out = a[np.ix_(ri, ci)].copy()
        out[~rv] = np.nan
        out[:, ~cv] = np.nan
        return out
    a, tr = wrap_geographic(a, tr, crs)
    # Restore the template longitude order after reprojection.
    tx = (template['lon'] + 180) % 360 - 180
    order = np.argsort(tx)
    target_tr = Affine(template['transform'].a, 0, tx[order[0]]-template['transform'].a*.5,
                       0, template['transform'].e, template['transform'].f)
    if a.shape == (template['height'], template['width']) and tr.almost_equals(target_tr) and crs == template['crs']:
        out = a
    else:
        out = np.full((template['height'], template['width']), np.nan, np.float32)
        reproject(source=a, destination=out, src_transform=tr, src_crs=crs,
                  dst_transform=target_tr, dst_crs=template['crs'], src_nodata=np.nan,
                  dst_nodata=np.nan, resampling=Resampling.nearest if categorical else Resampling.bilinear)
    inv = np.argsort(order)
    return out[:, inv]


def write_tif(path, array, template):
    with rasterio.open(path, 'w', **template['profile']) as dst:
        dst.write(np.where(np.isfinite(array), array, NODATA).astype('float32'), 1)


def monthly_average(var, year, spec, template, require_complete=False):
    """Average available months with equal weights; optionally require every month."""
    shape = (template['height'], template['width'])
    total = np.zeros(shape, float)
    count = np.zeros(shape, np.uint16)
    months = window_months(year, spec)
    for yy, mm in months:
        a = read_aligned(monthly_path(var, yy, mm), template, categorical=(var=='snowc'))
        ok = np.isfinite(a)
        total += np.where(ok, a, 0)
        count += ok
    keep = count == len(months) if require_complete else count > 0
    return np.divide(total, count, out=np.full(shape, np.nan), where=keep).astype('float32')


def snow_frequency(template):
    """Read snow frequency or calculate it from the monthly snow-cover rasters."""
    for path in FREQUENCY_CANDIDATES:
        if Path(path).is_file():
            log(f'Reading snow-cover frequency: {path}')
            a = read_aligned(path, template, categorical=True)
            good = np.isfinite(a)
            if np.any((a[good]<0) | (a[good]>1)):
                raise ValueError('Snow-cover frequency must be between 0 and 1.')
            return a
    log('Calculating snow-cover frequency on the common reference grid.')
    reference = template_from_path(monthly_path('snowc', FREQUENCY_YEARS[0], 1))
    shape = (reference['height'], reference['width'])
    n = np.zeros(shape, np.uint16)
    snow = np.zeros(shape, np.uint16)
    for year in FREQUENCY_YEARS:
        for month in range(1,13):
            a = read_aligned(monthly_path('snowc', year, month), reference, categorical=True)
            good = np.isfinite(a)
            n += good
            snow += good & (a>0)
        log(f'Snow-cover frequency: processed {year}')
    freq = np.divide(snow, n, out=np.full(shape, np.nan, np.float32), where=n>0)
    return align_mask_grid(freq, reference['transform'], template)


# 2. Annual inputs
def discover_ftc_files():
    ftcs = list(FTC_ROOT.rglob('*.tif'))
    ncs = list(SDE_ROOT.rglob('*.nc'))
    counts, depths = {}, {}
    for year in PARTIAL_YEARS:
        hits = [p for p in ftcs if re.search(r'complete_FTC_H3_'+str(year)+r'_M01-M07(?:\D|$)', p.stem, re.I)]
        if len(hits) != 1:
            raise ValueError(f'Expected one FTC raster for {year}; found {len(hits)}: {hits}')
        counts[year] = hits[0]
        for yy, mm in [(year-1,12)]+[(year,m) for m in range(1,5)]:
            if (yy,mm) not in depths:
                hits = [p for p in ncs if re.search(r'sde_'+str(yy)+f'_{mm:02d}_00UTC', p.name, re.I)]
                if len(hits) != 1:
                    raise ValueError(f'Expected one daily snow-depth file for {yy}-{mm:02d}; found {len(hits)}: {hits}')
                depths[yy,mm] = hits[0]
    return counts, depths


@contextmanager
def open_daily_netcdf(path):
    import xarray as xr
    with Path(path).open('rb') as fh:
        magic = fh.read(8)
        fh.seek(0)
        if magic == b'\x89HDF\r\n\x1a\n':
            engine, opts = 'h5netcdf', {}
        elif magic[:4] in (b'CDF\x01', b'CDF\x02'):
            engine, opts = 'scipy', {'mmap': False}
        else:
            raise ValueError(f'Unsupported NetCDF format: {path}')
        with xr.open_dataset(fh, engine=engine, backend_kwargs=opts, decode_times=True) as ds:
            yield ds


def daily_max_month(path, year, month, template):
    """Calculate the monthly maximum of complete daily 00 UTC snow depth."""
    record_source(path)
    with open_daily_netcdf(path) as ds:
        names = [n for n in ds.data_vars if n=='sde' or ds[n].attrs.get('GRIB_shortName')=='sde']
        if len(names) != 1:
            raise ValueError(f'Could not identify a unique sde variable: {path}')
        da = ds[names[0]]
        unit = str(da.attrs.get('units','')).strip().lower()
        factor = {'m':1., 'metre':1., 'meter':1., 'cm':.01}.get(unit)
        if factor is None:
            raise ValueError(f'Unsupported snow-depth unit {unit!r}: {path}')
        latn = next((n for n in ['latitude','lat'] if n in da.dims), None)
        lonn = next((n for n in ['longitude','lon'] if n in da.dims), None)
        tn = next((n for n in ['valid_time','time'] if n in da.dims), None)
        if None in [latn,lonn,tn]:
            raise ValueError(f'Unsupported NetCDF dimensions: {da.dims}')
        for dim in list(da.dims):
            if dim not in [latn,lonn,tn]:
                if da.sizes[dim] != 1:
                    raise ValueError(f'Unexpected non-singleton dimension: {dim}')
                da = da.isel({dim:0}, drop=True)
        times = pd.DatetimeIndex(da[tn].values)
        expected = pd.date_range(f'{year}-{month:02d}-01', periods=calendar.monthrange(year,month)[1], freq='D')
        if times.has_duplicates or not times.sort_values().equals(expected):
            raise ValueError(f'The daily series has missing dates or non-midnight UTC timestamps: {path}')
        ri, rv = coordinate_indices(da[latn].values, template['lat'])
        ci, cv = coordinate_indices(da[lonn].values, template['lon'], True)
        out = np.full((template['height'], template['width']), np.nan, np.float32)
        for r0 in range(0, template['height'], BLOCK_ROWS):
            r1 = min(r0+BLOCK_ROWS, template['height'])
            if not np.any(template['lat'][r0:r1]>=ANALYSIS_MIN_LAT):
                continue
            a = da.isel({latn:ri[r0:r1]}).transpose(tn,latn,lonn).values.astype(float)[:,:,ci]*factor
            a[:,~rv[r0:r1],:] = np.nan
            a[:,:,~cv] = np.nan
            a[(a<0)&(a>=-1e-12)] = 0
            a[a < -1e-12] = np.nan
            complete = np.all(np.isfinite(a), axis=0)
            mx = np.max(np.where(np.isfinite(a),a,-np.inf), axis=0)
            out[r0:r1] = np.where(complete, mx, np.nan)
    return out


def snowoff_path(year):
    return SNOW_STACK_DIR / f'SnowOff_PreSnowDepth_{year}_NH.tif'


def prepare_annual(template, domain, inputs):
    """Build annual raster stacks for the current run."""
    folder = RUN_DIR/'00_annual_inputs'
    folder.mkdir()
    nh = domain
    np.save(folder/'domain.npy', nh)
    qa = []
    for year in PARTIAL_YEARS:
        if INDICATOR in ('ST','SWC'):
            snow = monthly_average('sde', year, WINTER, template)
        elif INDICATOR == 'SnowOff':
            # Read snow-off date, PreMax and QC only.
            snow = read_aligned(snowoff_path(year), template, band=4, exact=True)
        else:
            snow = np.full(nh.shape, -np.inf, np.float32)
            complete = nh.copy()
            for yy,mm in [(year-1,12)]+[(year,m) for m in range(1,5)]:
                a = daily_max_month(inputs['sde'][yy,mm], yy, mm, template)
                complete &= np.isfinite(a)
                snow = np.fmax(snow,a)
            snow[~complete] = np.nan
        snow[(snow<0)&(snow>=-1e-12)] = 0
        snow[(snow<0)|~nh] = np.nan
        for period,spec in PERIODS.items():
            dest = folder/period
            dest.mkdir(exist_ok=True)
            covs = [monthly_average(v,year,spec,template,require_complete=INDICATOR in ('FTC','SnowOff'))
                    for v in CONTROLS]
            if INDICATOR in ('ST','SWC'):
                y = monthly_average(TARGET_VAR,year,spec,template)
            elif INDICATOR == 'SnowOff':
                y = read_aligned(snowoff_path(year), template, band=1, exact=True)
                qc = read_aligned(snowoff_path(year), template, band=7, exact=True)
                good = np.isfinite(qc) & (qc==1) & (y>=1) & (y<=366)
                y[~good] = np.nan
                snow[~good] = np.nan
            else:
                y = read_aligned(inputs['ftc'][year],template,exact=True)
                finite = np.isfinite(y)
                if np.any(y[finite]<0) or np.any(np.abs(y[finite]-np.rint(y[finite]))>1e-5):
                    raise ValueError(f'FTC counts for {year} must be nonnegative integers.')
            y[~nh] = np.nan
            fields = [snow, y] + covs
            dest_file = dest/f'{year}.npy'
            arr = np.lib.format.open_memmap(dest_file, mode='w+', dtype='float32',
                                           shape=(len(fields),)+nh.shape)
            for j,a in enumerate(fields):
                arr[j] = np.where(nh,a,np.nan)
            arr.flush()
            del arr
            valid = nh & np.all(np.isfinite(np.stack(fields)),axis=0)
            qa.append({'period':period,'year':year,'valid_pixels':int(valid.sum())})
            del covs,fields
        log(f'Prepared annual inputs for {year}')
        gc.collect()
    pd.DataFrame(qa).to_csv(RUN_DIR/'annual_input_QA.csv',index=False,encoding='utf-8-sig')
    return folder


def load_map_boundaries():
    """Load national boundaries for mapping."""
    global WORLD_GDF
    if not WORLD_SHP.is_file():
        raise FileNotFoundError(f'Country boundary shapefile not found: {WORLD_SHP}')
    try:
        import geopandas as gpd
    except ImportError as exc:
        raise ImportError('Install geopandas to draw country boundaries.') from exc
    WORLD_GDF = gpd.read_file(WORLD_SHP)
    if WORLD_GDF.empty:
        raise ValueError(f'Country boundary shapefile is empty: {WORLD_SHP}')
    if WORLD_GDF.crs is None:
        raise ValueError(f'Country boundary shapefile has no CRS: {WORLD_SHP}')
    WORLD_GDF = WORLD_GDF.to_crs(POLAR_CRS)


def preflight():
    """Check source files required to build annual raster stacks."""
    paths, missing, inputs = set(), [], {}
    if INDICATOR == 'FTC':
        inputs['ftc'], inputs['sde'] = discover_ftc_files()
        paths.update(inputs['ftc'].values())
        paths.update(inputs['sde'].values())
        base = inputs['ftc'][PARTIAL_YEARS[0]]
    elif INDICATOR == 'SnowOff':
        paths.update(snowoff_path(y) for y in PARTIAL_YEARS)
        base = snowoff_path(PARTIAL_YEARS[0])
    else:
        base = monthly_path(TARGET_VAR,PARTIAL_YEARS[0],1)
    for year in PARTIAL_YEARS:
        wanted = set()
        if INDICATOR in ('ST','SWC'):
            wanted.update(('sde',yy,mm) for yy,mm in window_months(year,WINTER))
        for spec in PERIODS.values():
            for var in CONTROLS + ([TARGET_VAR] if INDICATOR in ('ST','SWC') else []):
                wanted.update((var,yy,mm) for yy,mm in window_months(year,spec))
        for var,yy,mm in wanted:
            try:
                paths.add(monthly_path(var,yy,mm))
            except FileNotFoundError as exc:
                missing.append(str(exc))
    if not any(Path(p).is_file() for p in FREQUENCY_CANDIDATES):
        for year in FREQUENCY_YEARS:
            for month in range(1,13):
                try:
                    paths.add(monthly_path('snowc',year,month))
                except FileNotFoundError as exc:
                    missing.append(str(exc))
    missing += [str(p) for p in paths if not p.is_file()]
    if missing:
        raise FileNotFoundError(f'Missing {len(missing)} inputs (first 30):\n'+'\n'.join(missing[:30]))
    template = template_from_path(base)
    if INDICATOR=='SnowOff':
        for p in (snowoff_path(y) for y in PARTIAL_YEARS):
            with rasterio.open(p) as src:
                if src.count<7:
                    raise ValueError(f'Snow-off raster has fewer than seven bands: {p}')
    for path in paths:
        record_source(path)
    log(f'Input checks passed for {len(paths)} files.')
    return template, inputs


# 3. Partial correlation and AR(1) effective sample size
def lag1_columns(residuals, years):
    """Estimate lag-one residual correlations from consecutive year pairs."""
    adjacent = np.diff(np.asarray(years)) == 1
    out = np.full(residuals.shape[0],np.nan)
    if adjacent.sum() < MIN_ADJACENT_PAIRS:
        return out
    a = residuals[:,:-1][:,adjacent].copy()
    b = residuals[:,1:][:,adjacent].copy()
    a -= a.mean(axis=1,keepdims=True)
    b -= b.mean(axis=1,keepdims=True)
    den = np.sqrt(np.sum(a*a,axis=1)*np.sum(b*b,axis=1))
    out = np.divide(np.sum(a*b,axis=1),den,out=out,where=den>0)
    # Bound lag-one correlations away from singular endpoints.
    return np.clip(out,-.99,.99)


def effective_n_product(n, phi_x, phi_y):
    """Apply the two-series AR(1) approximation of Bretherton et al. (1999), Eq. 31."""
    product = np.asarray(phi_x)*np.asarray(phi_y)
    raw = float(n)*(1.0-product)/(1.0+product)
    used = np.minimum(float(n),raw) if CAP_NEFF_AT_N else raw.copy()
    return raw, used


PARTIAL_FIELDS = ['r','p','n_valid','n_eff','n_eff_raw','phi_x','phi_y','df','control_rank']


def partial_block(cube, years, domain):
    """Estimate partial correlations by grouping pixels with the same missing-year pattern."""
    t,d,p = cube.shape
    k = d-2
    valid = np.all(np.isfinite(cube),axis=1) & domain[None,:]
    nvalid = valid.sum(axis=0)
    result = {key:np.full(p,np.nan) for key in PARTIAL_FIELDS}
    result['n_valid'][domain] = nvalid[domain]
    minimum = t if REQUIRE_COMPLETE_PARTIAL else MIN_PARTIAL_YEARS
    eligible = np.flatnonzero(domain & (nvalid>=max(minimum,k+3)))
    if eligible.size==0:
        return result
    patterns, group = np.unique(valid[:,eligible].T,axis=0,return_inverse=True)
    for gid,pattern in enumerate(patterns):
        ti = np.flatnonzero(pattern)
        pix = eligible[group==gid]
        for start in range(0,len(pix),PARTIAL_PIXEL_CHUNK):
            ix = pix[start:start+PARTIAL_PIXEL_CHUNK]
            a = cube[ti][:,:,ix].transpose(2,0,1).astype(float)
            x,y,z = a[:,:,0],a[:,:,1],a[:,:,2:]
            z -= z.mean(axis=1,keepdims=True)
            scale = np.sqrt(np.mean(z*z,axis=1,keepdims=True))
            z /= np.where(scale>0,scale,1)
            design = np.concatenate([np.ones((len(ix),len(ti),1)),z],axis=2)
            u,s,vh = np.linalg.svd(design,full_matrices=False)
            good_s = s > (max(design.shape[1:])*np.finfo(float).eps*s[:,0:1])
            rank = good_s.sum(axis=1)
            ux = np.einsum('ntk,nt->nk',u,x)*good_s
            uy = np.einsum('ntk,nt->nk',u,y)*good_s
            rx = x-np.einsum('ntk,nk->nt',u,ux)
            ry = y-np.einsum('ntk,nk->nt',u,uy)
            rx -= rx.mean(axis=1,keepdims=True)
            ry -= ry.mean(axis=1,keepdims=True)
            sx,sy = np.sum(rx*rx,axis=1),np.sum(ry*ry,axis=1)
            # Require residual variation relative to the original series.
            ssx = np.sum((x-x.mean(axis=1,keepdims=True))**2,axis=1)
            ssy = np.sum((y-y.mean(axis=1,keepdims=True))**2,axis=1)
            estimable = (sx>np.maximum(ssx,1e-30)*1e-12)&(sy>np.maximum(ssy,1e-30)*1e-12)
            r = np.divide(np.sum(rx*ry,axis=1),np.sqrt(sx*sy),out=np.full(len(ix),np.nan),where=estimable)
            r = np.clip(r,-1,1)
            px,py = lag1_columns(rx,np.asarray(years)[ti]),lag1_columns(ry,np.asarray(years)[ti])
            raw,neff = effective_n_product(len(ti),px,py)
            df = neff-(rank-1)-2
            ok = np.isfinite(r)&np.isfinite(df)&(df>0)
            pv = np.full(len(ix),np.nan)
            tt = np.abs(r[ok])*np.sqrt(df[ok]/np.maximum(1-r[ok]**2,1e-15))
            pv[ok] = 2*stats.t.sf(tt,df[ok])
            for key,values in [('r',r),('p',pv),('n_eff',neff),('n_eff_raw',raw),
                               ('phi_x',px),('phi_y',py),('df',df),('control_rank',rank-1)]:
                result[key][ix] = values
    return result


# 4. Panel model and spatially clustered inference
def pixel_demean(a, valid):
    """Demean each variable across valid years within each pixel."""
    n = valid.sum(axis=1)
    mean = np.sum(np.where(valid[:,:,None],a,0),axis=1)/np.maximum(n[:,None],1)
    return np.where(valid[:,:,None],a-mean[:,None,:],0)


def panel_chunk(cube, years, row0, width, template):
    """Construct valid pixel-year observations and directional snow terms."""
    ids = [years.index(y) for y in PANEL_YEARS]
    if INDICATOR=='SnowOff':
        prev = [years.index(y-1) for y in PANEL_YEARS]
        snow = (cube[ids,0]-cube[prev,0]).T*100.
        response = (cube[ids,1]-cube[prev,1]).T
        controls = cube[ids,2:].transpose(2,0,1)
    else:
        raw = cube[ids].transpose(2,0,1).astype(float)
        good = np.all(np.isfinite(raw),axis=2)
        anomalies = pixel_demean(raw,good)
        snow = anomalies[:,:,0]*100.
        response = raw[:,:,1]
        controls = raw[:,:,2:]
        response = np.where(good,response,np.nan)
    valid = np.isfinite(snow)&np.isfinite(response)&np.all(np.isfinite(controls),axis=2)
    valid &= (np.abs(snow)>0)&(np.abs(snow)<=MAX_AMPLITUDE_CM)
    valid &= (valid.sum(axis=1)>=MIN_PANEL_OBSERVATIONS)[:,None]
    keep = np.flatnonzero(valid.any(axis=1))
    if not len(keep):
        return None
    valid,snow,response,controls = valid[keep],snow[keep],response[keep],controls[keep]
    # Use year differences for Snow-off and pixel anomalies for the other indicators.
    if INDICATOR=='SnowOff':
        raw_y = response.copy()
    else:
        raw_y = anomalies[keep,:,1]
    inc,dec = np.maximum(snow,0),np.maximum(-snow,0)
    fields = np.concatenate([response[:,:,None],inc[:,:,None],dec[:,:,None],controls],axis=2)
    w = pixel_demean(fields,valid)
    row = row0+keep//width
    col = keep%width
    lon = (template['lon'][col]+180)%360-180
    lat = template['lat'][row]
    cluster = np.floor(lat).astype(int)*360+np.floor(lon+180).astype(int)
    return {'valid':valid,'w':w,'snow':snow,'raw_y':raw_y,'cluster':cluster,
            'pixel_id':row.astype(np.int64)*width+col}


def year_schur(paths, t, d):
    """Absorb pixel and year fixed effects in blocks using a Schur complement."""
    lap = np.zeros((t,t))
    rhs = np.zeros((t,d))
    nobs = npix = 0
    for path in paths:
        with np.load(path) as a:
            v,w = a['valid'],a['w']
            vf = v.astype(float)
            n = vf.sum(axis=1)
            lap += np.diag(vf.sum(axis=0))-vf.T@(vf/n[:,None])
            rhs += w.sum(axis=0)
            nobs += int(v.sum())
            npix += len(v)

    gamma = np.linalg.pinv(lap,hermitian=True)@rhs
    rank = int(np.linalg.matrix_rank(lap))
    return gamma,nobs,npix,rank


def within_from_chunk(a, gamma):
    v,w = a['valid'],a['w']
    avg = v.astype(float)@gamma/v.sum(axis=1)[:,None]
    return np.where(v[:,:,None],w-gamma[None,:,:]+avg[:,None,:],0)


def fit_panel(paths, t, period, out):
    k = 2+len(CONTROLS)
    gamma,nobs,npix,year_rank = year_schur(paths,t,k+1)
    gram = np.zeros((k,k)); cross = np.zeros(k)
    clusters = set()
    for path in paths:
        with np.load(path) as a:
            w = within_from_chunk(a,gamma)
            q = w[a['valid']]
            gram += q[:,1:].T@q[:,1:]
            cross += q[:,1:].T@q[:,0]
            clusters.update(a['cluster'].tolist())
    scales = np.sqrt(np.maximum(np.diag(gram)/nobs,0))
    export_vif(gram, scales, nobs, npix, period, out)
    if nobs < MIN_MODEL_N:
        raise ValueError(f'Too few valid panel observations: {nobs}<{MIN_MODEL_N}')
    if np.any(scales<=1e-14):
        raise ValueError('At least one predictor is constant after absorbing fixed effects.')
    gscaled = gram/np.outer(scales,scales)
    if np.linalg.matrix_rank(gscaled)<k:
        raise ValueError('The fixed-effect-transformed predictor matrix is rank deficient.')
    bread = np.linalg.inv(gscaled)
    beta_scaled = bread@(cross/scales)
    beta = beta_scaled/scales
    all_clusters = np.array(sorted(clusters))
    g = len(all_clusters)
    if g<MIN_MODEL_CLUSTERS:
        raise ValueError(f'Too few spatial clusters: {g}')
    scores = np.zeros((g,k)); sse = 0.
    max_year_mean = np.zeros((t,k+1)); max_pixel_mean = 0.
    for path in paths:
        with np.load(path) as a:
            w = within_from_chunk(a,gamma)
            v = a['valid']
            e = w[:,:,0]-w[:,:,1:]@beta
            e[~v] = 0
            sse += float(np.sum(e*e))
            score = np.einsum('ptk,pt->pk',w[:,:,1:]/scales,e)
            np.add.at(scores,np.searchsorted(all_clusters,a['cluster']),score)
            max_year_mean += w.sum(axis=0)
            max_pixel_mean = max(max_pixel_mean,float(np.max(np.abs(w.sum(axis=1)))))
    absorbed_rank = npix+year_rank
    k_correction = k+absorbed_rank if CR1_COUNT_ABSORBED_FE else k
    if nobs<=k+absorbed_rank:
        raise ValueError('The fixed-effects model has no residual degrees of freedom.')
    cr1 = g/(g-1)*(nobs-1)/(nobs-k_correction)
    cov_scaled = bread@(scores.T@scores)@bread*cr1
    cov = cov_scaled/np.outer(scales,scales)
    se = np.sqrt(np.maximum(np.diag(cov),0))
    df = g-1
    p = 2*stats.t.sf(np.abs(np.divide(beta,se,out=np.full(k,np.nan),where=se>0)),df)
    return {'beta':beta,'cov':cov,'se':se,'p':p,'df':df,'gamma':gamma,
            'n_obs':nobs,'n_pixels':npix,'n_clusters':g,'absorbed_fe_rank':absorbed_rank,
            'cr1_factor':cr1,'condition_number':float(np.linalg.cond(gscaled)),
            'max_abs_pixel_within_sum':max_pixel_mean,
            'max_abs_year_within_sum':float(np.max(np.abs(max_year_mean))),
            'rmse':math.sqrt(sse/(nobs-k-absorbed_rank))}


def model_contrasts(fit):
    b,cov,df = fit['beta'],fit['cov'],fit['df']
    contrast = b[0]+b[1]
    se = math.sqrt(max(cov[0,0]+cov[1,1]+2*cov[0,1],0))
    p = 2*stats.t.sf(abs(contrast/se),df) if se>0 else np.nan
    rho = rho_se = low = high = np.nan
    if abs(b[0])>1e-12:
        rho = abs(b[1])/abs(b[0])
        if abs(b[1])>1e-12:
            grad = np.zeros(len(b))
            grad[0] = -abs(b[1])*np.sign(b[0])/abs(b[0])**2
            grad[1] = np.sign(b[1])/abs(b[0])
            rho_se = math.sqrt(max(float(grad@cov@grad),0))
            crit = stats.t.ppf(.975,df)
            low,high = max(0,rho-crit*rho_se),rho+crit*rho_se
    return {'beta_inc':b[0],'beta_dec':b[1],'beta_inc_plus_beta_dec':contrast,
            'se_inc_plus_dec':se,'p_asymmetry_cluster':p,'rho_snow':rho,
            'rho_se_delta':rho_se,'rho_ci_low_delta':low,'rho_ci_high_delta':high}


def new_bin_accumulator():
    return {(method,sign,b):[0,0.,0.,0.,set()] for method in ['raw','adjusted']
            for sign in ['increase','decrease'] for b in range(len(AMP_CENTERS_CM))}


def update_bins(acc, method, snow, y, valid, clusters):
    amp = np.abs(snow)
    bins = np.searchsorted(AMP_EDGES_CM,amp,side='right')-1
    bins[amp==MAX_AMPLITUDE_CM] = len(AMP_CENTERS_CM)-1
    for sign,sgn in [('increase',snow>0),('decrease',snow<0)]:
        for b in range(len(AMP_CENTERS_CM)):
            mask = valid&sgn&(bins==b)&np.isfinite(y)
            n = int(mask.sum())
            if n:
                vals = y[mask]
                q = acc[method,sign,b]
                q[0] += n; q[1] += float(vals.sum()); q[2] += float(np.sum(vals*vals))
                q[3] += float(amp[mask].sum())
                q[4].update(clusters[mask.any(axis=1)].tolist())


def final_bins(acc):
    rows=[]
    for (method,sign,b),(n,total,ss,amplitude,clusters) in acc.items():
        mean = total/n if n else np.nan
        # Report descriptive standard errors within bins.
        sd = math.sqrt(max((ss-total*total/n)/(n-1),0)) if n>1 else np.nan
        se = sd/math.sqrt(n) if n else np.nan
        ci95 = 1.96*se if np.isfinite(se) else np.nan
        pvalue = (2*stats.t.sf(abs(mean/se), n-1)
                  if n >= 2 and np.isfinite(se) and se > 0 else np.nan)
        rows.append({'method':method,'direction':sign,'center_cm':AMP_CENTERS_CM[b],
                     'lower_cm':AMP_EDGES_CM[b],'upper_cm':AMP_EDGES_CM[b+1],
                     'n_pixel_years':n,'n_spatial_clusters':len(clusters),
                     'actual_mean_amplitude_cm':amplitude/n if n else np.nan,
                     'mean_response':mean,'SD':sd,'SE_descriptive':se,
                     'CI_low_descriptive':mean-ci95 if np.isfinite(ci95) else np.nan,
                     'CI_high_descriptive':mean+ci95 if np.isfinite(ci95) else np.nan,
                     'ci95':ci95,'pvalue':pvalue,'meets_min_n':n>=MIN_BIN_N})
    return pd.DataFrame(rows)


def bin_ratios(bins):
    rows=[]
    for method in ['raw','adjusted']:
        for amp in AMP_CENTERS_CM:
            sub=bins[(bins.method==method)&(bins.center_cm==amp)].set_index('direction')
            a,b=sub.loc['increase'],sub.loc['decrease']
            ok=a.meets_min_n and b.meets_min_n and abs(a.mean_response)>1e-12
            rows.append({'method':method,'center_cm':amp,'n_increase':int(a.n_pixel_years),
                         'n_decrease':int(b.n_pixel_years),'mean_increase':a.mean_response,
                         'mean_decrease':b.mean_response,
                         'actual_amplitude_increase_cm':a.actual_mean_amplitude_cm,
                         'actual_amplitude_decrease_cm':b.actual_mean_amplitude_cm,
                         'ratio_abs_decrease_over_increase':abs(b.mean_response)/abs(a.mean_response) if ok else np.nan,
                         'ratio_valid':bool(ok)})
    return pd.DataFrame(rows)


def analyze_panel(paths, period, out):
    fit=fit_panel(paths,len(PANEL_YEARS),period,out)
    asym=model_contrasts(fit)
    names=['snow_increase_cm','snow_decrease_cm']+CONTROLS
    crit=stats.t.ppf(.975,fit['df'])
    coef=pd.DataFrame({'variable':names,'beta':fit['beta'],'SE_cluster':fit['se'],
                       'p_cluster':fit['p'],'CI_low_cluster':fit['beta']-crit*fit['se'],
                       'CI_high_cluster':fit['beta']+crit*fit['se']})
    coef.to_csv(out/'FE_coefficients.csv',index=False,encoding='utf-8-sig')
    pd.DataFrame(fit['cov'],index=names,columns=names).to_csv(out/'FE_covariance_cluster.csv',encoding='utf-8-sig')
    # Table S9 compatibility output: preserve the established field names and order.
    summary = {
        'Trait': INDICATOR,
        'Season': period,
        'beta_inc': fit['beta'][0],
        'beta_dec': fit['beta'][1],
        'se_inc_cluster': fit['se'][0],
        'se_dec_cluster': fit['se'][1],
        'p_inc_cluster': fit['p'][0],
        'p_dec_cluster': fit['p'][1],
        'contrast_beta_inc_plus_beta_dec': asym['beta_inc_plus_beta_dec'],
        'p_asymmetry_cluster': asym['p_asymmetry_cluster'],
        'rho_snow': asym['rho_snow'],
        'rho_se': asym['rho_se_delta'],
        'rho_ci_low': asym['rho_ci_low_delta'],
        'rho_ci_high': asym['rho_ci_high_delta'],
    }
    table_s9_columns = ['Trait','Season','beta_inc','beta_dec','se_inc_cluster','se_dec_cluster',
                        'p_inc_cluster','p_dec_cluster','contrast_beta_inc_plus_beta_dec',
                        'p_asymmetry_cluster','rho_snow','rho_se','rho_ci_low','rho_ci_high']
    pd.DataFrame([summary], columns=table_s9_columns).to_csv(
        out/'FE_model_summary.csv', index=False, encoding='utf-8-sig')

    # Keep model diagnostics in a separate file so the established summary schema is unchanged.
    diagnostics = {key: fit[key] for key in ['n_obs','n_pixels','n_clusters','absorbed_fe_rank','df',
                  'cr1_factor','condition_number','max_abs_pixel_within_sum',
                  'max_abs_year_within_sum','rmse']}
    diagnostics.update({'indicator':INDICATOR,'period':period,'controls':'+'.join(CONTROLS),
                        'first_year':PANEL_YEARS[0],'last_year':PANEL_YEARS[-1],
                        'adjustment_definition':ADJUSTED_RESPONSE_DEFINITION,
                        'beta_inc_plus_beta_dec':asym['beta_inc_plus_beta_dec'],
                        'se_inc_plus_dec':asym['se_inc_plus_dec'],
                        'p_asymmetry_cluster':asym['p_asymmetry_cluster'],
                        'rho_snow':asym['rho_snow'],'rho_se_delta':asym['rho_se_delta'],
                        'rho_ci_low_delta':asym['rho_ci_low_delta'],
                        'rho_ci_high_delta':asym['rho_ci_high_delta']})
    pd.DataFrame([diagnostics]).to_csv(
        out/'FE_model_diagnostics.csv', index=False, encoding='utf-8-sig')
    acc=new_bin_accumulator()
    for path in paths:
        with np.load(path) as a:
            w=within_from_chunk(a,fit['gamma'])
            error=w[:,:,0]-w[:,:,1:]@fit['beta']
            # Add the original snow-term contribution back to the full-model residual.
            snow=a['snow']
            adjusted=error+np.maximum(snow,0)*fit['beta'][0]+np.maximum(-snow,0)*fit['beta'][1]
            update_bins(acc,'raw',a['snow'],a['raw_y'],a['valid'],a['cluster'])
            update_bins(acc,'adjusted',a['snow'],adjusted,a['valid'],a['cluster'])
    bins=final_bins(acc)
    ratios=bin_ratios(bins)
    # Keep the established binned-response fields and include method (raw/adjusted)
    # and Season so all periods can be combined into one indicator-level table.
    binned_export=bins.rename(columns={'direction':'sign_group','center_cm':'amp_bin_cm',
        'n_pixel_years':'n','mean_response':'mean','SE_descriptive':'se'}).copy()
    binned_export['Season']=period
    binned_export=binned_export[['method','sign_group','amp_bin_cm','n','mean','se','ci95','pvalue','Season']]
    binned_export.to_csv(out/'equal_amplitude_binned_response.csv',index=False,encoding='utf-8-sig')
    ratios.insert(0,'Season',period)
    ratios.to_csv(out/'equal_amplitude_binned_ratio.csv',index=False,encoding='utf-8-sig')
    predicted=[]
    for amp in range(1,61):
        predicted.append({'amplitude_cm':amp,'predicted_increase_response':amp*fit['beta'][0],
                          'predicted_decrease_response':amp*fit['beta'][1],
                          'increase_CI_low':amp*coef.CI_low_cluster.iloc[0],
                          'increase_CI_high':amp*coef.CI_high_cluster.iloc[0],
                          'decrease_CI_low':amp*coef.CI_low_cluster.iloc[1],
                          'decrease_CI_high':amp*coef.CI_high_cluster.iloc[1],**asym})
    predictions=pd.DataFrame(predicted)
    predictions.insert(0,'Season',period)
    predictions.to_csv(out/'equal_amplitude_model_prediction_1_60cm.csv',index=False,encoding='utf-8-sig')
    draw_response_figures(bins,ratios,predictions,out)
    log(f'{period}: fitted the equal-amplitude model with {fit["n_pixels"]} pixels, {fit["n_obs"]} observations, and {fit["n_clusters"]} spatial clusters.')
    return summary


# 5. Maps and response plots
def save_figure(fig, base):
    """Save figures at the configured size."""
    for ext in ['png','pdf','tif']:
        opts={'pil_kwargs':{'compression':'tiff_lzw'}} if ext=='tif' else {}
        fig.savefig(str(base)+'.'+ext,dpi=FIG_DPI,facecolor='white',**opts)
    plt.close(fig)


def polar_map(r_path, template, base):
    """Render the partial correlation map in a North Polar Stereographic projection."""
    with rasterio.open(r_path) as src:
        r=src.read(1,masked=True).filled(np.nan)
    # Apply the latitude cutoff for display only.
    r[template['lat']<MAP_MIN_LAT]=np.nan
    r,tr=wrap_geographic(r,template['transform'],template['crs'])
    xx,yy=transform_points('EPSG:4326',POLAR_CRS,[0],[MAP_MIN_LAT])
    radius=math.hypot(xx[0],yy[0])
    target=np.full((MAP_RASTER_PIXELS,MAP_RASTER_PIXELS),np.nan,np.float32)
    reproject(r,target,src_transform=tr,src_crs=template['crs'],src_nodata=np.nan,
              dst_transform=from_bounds(-radius,-radius,radius,radius,MAP_RASTER_PIXELS,MAP_RASTER_PIXELS),
              dst_crs=POLAR_CRS,dst_nodata=np.nan,resampling=Resampling.nearest)

    fig=plt.figure(figsize=FIGSIZE,facecolor=MAP_BACKGROUND_COLOR)
    ax=fig.add_axes(MAP_AX_RECT,facecolor=MAP_BACKGROUND_COLOR)
    ax.set_aspect('equal',adjustable='box')
    ax.set_xlim(-radius,radius);ax.set_ylim(-radius,radius)
    clip_circle=Circle((0,0),radius,facecolor=MAP_BACKGROUND_COLOR,edgecolor='none',linewidth=0,zorder=0)
    ax.add_patch(clip_circle)

    image=ax.imshow(target,extent=[-radius,radius,-radius,radius],origin='upper',
                    cmap=CORR_CMAP,vmin=-1,vmax=1,zorder=2,interpolation='nearest')
    image.set_clip_path(clip_circle)


    for lat in MAP_PARALLELS:
        lon=np.linspace(-180,180,721)
        gx,gy=transform_points('EPSG:4326',POLAR_CRS,lon.tolist(),[lat]*len(lon))
        ax.plot(gx,gy,color=MAP_GRID_COLOR,lw=.45,zorder=3)[0].set_clip_path(clip_circle)
    for lon in MAP_MERIDIANS:
        lats=np.linspace(MAP_MIN_LAT,90,201)
        gx,gy=transform_points('EPSG:4326',POLAR_CRS,[lon]*len(lats),lats.tolist())
        ax.plot(gx,gy,color=MAP_GRID_COLOR,lw=.45,zorder=3)[0].set_clip_path(clip_circle)

    if WORLD_GDF is None or WORLD_GDF.empty:
        raise RuntimeError('Country boundaries must be loaded before rendering the map.')
    WORLD_GDF.boundary.plot(ax=ax,color=MAP_COUNTRY_COLOR,linewidth=.55,zorder=5)

    for collection in ax.collections:
        collection.set_clip_path(clip_circle)


    outer=Circle((0,0),radius,facecolor='none',edgecolor=MAP_BORDER_COLOR,linewidth=1.0,zorder=6)
    ax.add_patch(outer)
    ax.set_xlim(-radius,radius);ax.set_ylim(-radius,radius);ax.set_axis_off()
    save_figure(fig,base)

def style_axes(ax):
    ax.spines['top'].set_visible(False);ax.spines['right'].set_visible(False)
    ax.tick_params(direction='out')
    ax.set_xlabel('Snow-depth change magnitude (cm)')
    ax.set_xlim(0,62)


def draw_response_figures(bins,ratios,pred,out):
    for method in ['raw','adjusted']:
        fig,ax=plt.subplots(figsize=FIGSIZE)
        fig.subplots_adjust(left=.15,right=.95,bottom=.13,top=.91)
        for sign,color,label in [('increase',POSITIVE_COLOR,'Increased snowpack thickness'),
                                  ('decrease',NEGATIVE_COLOR,'Decreased snowpack thickness')]:
            sub=bins[(bins.method==method)&(bins.direction==sign)&bins.meets_min_n]
            ax.errorbar(sub.center_cm,sub.mean_response,yerr=1.96*sub.SE_descriptive,
                        fmt='o-',color=color,label=label,capsize=3)
        ax.axhline(0,color='#888888',lw=.8)
        ax.set_ylabel(RESPONSE_LABEL)
        ax.set_title('Observed response' if method=='raw' else 'Covariate- and FE-adjusted response')
        style_axes(ax);ax.legend(frameon=False,fontsize=10)
        save_figure(fig,out/f'{method}_binned_response')
    fig,ax=plt.subplots(figsize=FIGSIZE)
    fig.subplots_adjust(left=.15,right=.95,bottom=.13,top=.91)
    for method,color,label in [('raw','#777777','Observed'),('adjusted','#222222','Adjusted')]:
        sub=ratios[ratios.method==method]
        ax.plot(sub.center_cm,sub.ratio_abs_decrease_over_increase,'o-',color=color,label=label)
    ax.axhline(1,color='#AAAAAA',lw=.8,ls='--')
    ax.set_ylabel('|Decreased response| / |Increased response|')
    style_axes(ax);ax.legend(frameon=False)
    save_figure(fig,out/'raw_adjusted_binned_ratio')
    fig,ax=plt.subplots(figsize=FIGSIZE)
    fig.subplots_adjust(left=.15,right=.95,bottom=.13,top=.91)
    for label,color,key,prefix in [('Increased snowpack thickness',POSITIVE_COLOR,'predicted_increase_response','increase'),
                                    ('Decreased snowpack thickness',NEGATIVE_COLOR,'predicted_decrease_response','decrease')]:
        ax.plot(pred.amplitude_cm,pred[key],color=color,label=label)
        ax.fill_between(pred.amplitude_cm,pred[prefix+'_CI_low'],pred[prefix+'_CI_high'],color=color,alpha=.2)
    ax.axhline(0,color='#888888',lw=.8)
    ax.set_ylabel(RESPONSE_LABEL);ax.set_title('Model response at identical magnitudes')
    style_axes(ax);ax.legend(frameon=False,fontsize=10)
    save_figure(fig,out/'model_equal_amplitude_response')


def glacier_mask_path(requested):
    """Use a local glacier mask or download the ECMWF file."""
    path = Path(requested) if requested is not None else GLACIER_MASK_PATH
    if path.is_file():
        return path
    if requested is not None:
        raise FileNotFoundError(f'Specified glacier mask not found: {path}')
    from tempfile import NamedTemporaryFile
    from urllib.request import urlopen
    from shutil import copyfileobj

    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        log(f'Downloading the ERA5-Land glacier mask: {path}')
        with NamedTemporaryFile(dir=path.parent, suffix='.part', delete=False) as handle:
            temporary = Path(handle.name)
            with urlopen(GLACIER_MASK_URL, timeout=45) as source:
                copyfileobj(source, handle)
        with temporary.open('rb') as handle:
            magic = handle.read(8)
        if magic[:4] not in (b'CDF\x01', b'CDF\x02') and magic != b'\x89HDF\r\n\x1a\n':
            raise ValueError('The downloaded file is not NetCDF.')
        temporary.replace(path)
    except Exception as exc:
        raise RuntimeError(
            f'Could not download the glacier mask to {path}. '
            'Alternatively, specify a local file with --glacier-mask.'
        ) from exc
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
    return path


def glacier_fraction(path, template):
    """Align glacier fractions to the analysis grid by nearest neighbor."""
    if path.suffix.lower() in ('.tif', '.tiff'):
        values = read_aligned(path, template, categorical=True)
    else:
        record_source(path)
        with open_daily_netcdf(path) as ds:
            lat_name = next((key for key in ('latitude', 'lat') if key in ds.coords), None)
            lon_name = next((key for key in ('longitude', 'lon') if key in ds.coords), None)
            if lat_name is None or lon_name is None:
                raise ValueError(f'The glacier file has no latitude or longitude coordinates: {path}')
            candidates = [key for key, var in ds.data_vars.items()
                          if set(var.squeeze(drop=True).dims) == {lat_name, lon_name}]
            known = [key for key in candidates
                     if key.lower() in ('glm', 'glacier_mask', 'cicecap', 'var260294')
                     or str(ds[key].attrs.get('GRIB_shortName', '')).lower() == 'glm']
            fields = known or candidates
            if len(fields) != 1:
                raise ValueError(f'Could not identify a unique two-dimensional glacier fraction variable in {path}: {fields}')
            grid = ds[fields[0]].squeeze(drop=True).transpose(lat_name, lon_name)
            lat = np.asarray(grid[lat_name], dtype=float)
            lon = np.asarray(grid[lon_name], dtype=float)
            data = np.asarray(grid, dtype=np.float32)
        values = match_glacier_coordinates(data, lat, lon, template)
    finite = np.isfinite(values)
    if not finite.any() or np.any((values[finite] < -1e-6) | (values[finite] > 1 + 1e-6)):
        raise ValueError(f'Glacier fractions must be between 0 and 1: {path}')
    return np.clip(values, 0, 1)


def nearest_axis_cells(source, target):
    """Find nearest grid centers and flag targets outside the grid coverage."""
    right = np.clip(np.searchsorted(source, target), 0, len(source)-1)
    left = np.maximum(right-1, 0)
    indices = np.where(np.abs(target-source[left]) < np.abs(source[right]-target),
                       left, right)
    lower = source[0] - (source[1]-source[0])/2
    upper = source[-1] + (source[-1]-source[-2])/2
    return indices, (target >= lower) & (target <= upper)


def match_glacier_coordinates(data, lat, lon, template):
    """Match NetCDF grid centers using their actual coordinates and cyclic longitude."""
    if (lat.ndim != 1 or lon.ndim != 1 or min(len(lat), len(lon)) < 2
            or data.shape != (len(lat), len(lon))
            or not (np.isfinite(lat).all() and np.isfinite(lon).all())
            or np.any(np.abs(lat) > 90.000001)):
        raise ValueError('The glacier NetCDF requires valid one-dimensional latitude and longitude axes.')

    lat_order = np.argsort(lat, kind='stable')
    sorted_lat = lat[lat_order]
    if np.any(np.diff(sorted_lat) <= 0):
        raise ValueError('The glacier NetCDF contains duplicate latitude coordinates.')
    rows, valid_rows = nearest_axis_cells(sorted_lat, template['lat'])

    cyclic_lon = np.mod(lon, 360.0)
    lon_order = np.argsort(cyclic_lon, kind='stable')
    sorted_lon = cyclic_lon[lon_order]
    keep = np.r_[True, np.diff(sorted_lon) > 1e-6]
    lon_order, sorted_lon = lon_order[keep], sorted_lon[keep]
    if len(sorted_lon) > 1 and sorted_lon[0] + 360 - sorted_lon[-1] <= 1e-6:
        lon_order, sorted_lon = lon_order[:-1], sorted_lon[:-1]
    if len(sorted_lon) < 2:
        raise ValueError('The glacier NetCDF requires at least two distinct longitudes.')

    gaps = np.diff(np.r_[sorted_lon, sorted_lon[0] + 360])
    seam = int(np.argmax(gaps))
    if gaps[seam] <= 1.5 * np.median(gaps):
        search_lon = np.r_[sorted_lon[-1]-360, sorted_lon, sorted_lon[0]+360]
        search_order = np.r_[lon_order[-1], lon_order, lon_order[0]]
        cols, valid_cols = nearest_axis_cells(search_lon, np.mod(template['lon'], 360.0))
    else:
        # Unwrap regional coverage across the longitude seam.
        start = (seam + 1) % len(sorted_lon)
        search_order = np.roll(lon_order, -start)
        origin = sorted_lon[start]
        search_lon = origin + np.mod(cyclic_lon[search_order]-origin, 360)
        midpoint = (search_lon[0] + search_lon[-1])/2
        targets = midpoint + (template['lon']-midpoint+180) % 360 - 180
        cols, valid_cols = nearest_axis_cells(search_lon, targets)

    aligned = data[np.ix_(lat_order[rows], search_order[cols])].copy()
    aligned[~valid_rows, :] = np.nan
    aligned[:, ~valid_cols] = np.nan
    return aligned


def align_mask_grid(data, transform, template):
    """Align a static field while preserving the template longitude order."""
    crs = rasterio.crs.CRS.from_epsg(4326)
    data, transform = wrap_geographic(data, transform, crs)
    x = (template['lon'] + 180) % 360 - 180
    order = np.argsort(x)
    target_tr = Affine(template['transform'].a, 0, x[order[0]]-template['transform'].a/2,
                       0, template['transform'].e, template['transform'].f)
    dest = np.full((template['height'], template['width']), np.nan, np.float32)
    reproject(data, dest, src_transform=transform, src_crs=crs, src_nodata=np.nan,
              dst_transform=target_tr, dst_crs=template['crs'], dst_nodata=np.nan,
              resampling=Resampling.nearest)
    return dest[:, np.argsort(order)]


def analysis_domain(template, freq, ice):
    north = np.broadcast_to(template['lat'][:, None] >= ANALYSIS_MIN_LAT, freq.shape)
    snowy = north & np.isfinite(freq) & (freq >= SNOW_FREQUENCY_THRESHOLD)
    known = np.isfinite(ice)
    if not np.any(snowy & known):
        raise ValueError('The glacier mask does not overlap the Northern Hemisphere snow domain. Check its CRS and coverage.')
    excluded = snowy & known & (ice >= GLACIER_FRACTION_THRESHOLD)
    domain = snowy & known & ~excluded
    if not domain.any():
        raise ValueError('The common spatial mask retained no pixels.')
    counts = {
        'indicator': INDICATOR, 'n_northern_pixels': int(north.sum()),
        'n_snow_frequency_pass': int(snowy.sum()),
        'n_excluded_glacier': int(excluded.sum()),
        'n_missing_glacier_fraction': int((snowy & ~known).sum()),
        'n_domain_pixels': int(domain.sum()),
        'snow_frequency_threshold': SNOW_FREQUENCY_THRESHOLD,
        'glacier_fraction_exclusion_threshold': GLACIER_FRACTION_THRESHOLD,
    }
    return domain, counts


def partial_paths(run, period):
    return {key: run/period/'partial_rasters'/f'{INDICATOR}_{period}_{key}.tif'
            for key in PARTIAL_FIELDS}


def same_grid(path, template):
    with rasterio.open(path) as src:
        return (src.count == 1 and src.crs == template['crs']
                and (src.height, src.width) == (template['height'], template['width'])
                and src.transform.almost_equals(template['transform'], precision=1e-8))


def annual_cache(run):
    """Validate the previous annual cache against the current analysis settings."""
    config = json.loads((run/'run_config.json').read_text(encoding='utf-8-sig'))
    expected = {'indicator': INDICATOR, 'controls': CONTROLS, 'periods': PERIODS,
                'partial_years': PARTIAL_YEARS}
    if any(config.get(key) != value for key, value in expected.items()):
        raise ValueError('The cached indicator, covariates, periods, or partial correlation years differ.')
    folder = run/'00_annual_inputs'
    reference = next((partial_paths(run, period)['r'] for period in PERIODS
                      if partial_paths(run, period)['r'].is_file()), None)
    if reference is None:
        reference = next((path for path in (run/'analysis_domain.tif', folder/'snow_frequency.tif')
                          if path.is_file()), run/'analysis_domain.tif')
    template = template_from_path(reference)
    shape = template['height'], template['width']
    old_domain = np.load(folder/'domain.npy', mmap_mode='r', allow_pickle=False)
    if old_domain.shape != shape or old_domain.dtype != np.bool_:
        raise ValueError('The cached domain mask does not match the raster grid.')
    for period in PERIODS:
        for year in PARTIAL_YEARS:
            array = np.load(folder/period/f'{year}.npy', mmap_mode='r', allow_pickle=False)
            if array.shape != (2 + len(CONTROLS), *shape) or array.dtype != np.float32:
                raise ValueError(f'Incompatible cached annual stack: {period}/{year}')
            del array
    return folder, template, old_domain, config


def reusable_partials(run, template, config):
    expected = {
        'effective_n_formula': 'n*(1-phi_x*phi_y)/(1+phi_x*phi_y)',
        'cap_neff_at_n': CAP_NEFF_AT_N, 'min_adjacent_pairs': MIN_ADJACENT_PAIRS,
        'require_complete_partial': REQUIRE_COMPLETE_PARTIAL,
        'min_partial_years': MIN_PARTIAL_YEARS,
    }
    if any(config.get(key) != value for key, value in expected.items()):
        return {}
    reusable = {}
    completed = set()
    summary_path = run/'NH_all_periods_partial_summary.csv'
    if summary_path.is_file():
        summary = pd.read_csv(summary_path, encoding='utf-8-sig')
        if 'period' in summary:
            completed.update(summary['period'].dropna())
    for period in PERIODS:
        if period not in completed and not (run/period/'partial_summary_NH.csv').is_file():
            continue
        paths = partial_paths(run, period)
        try:
            if all(path.is_file() and same_grid(path, template) for path in paths.values()):
                reusable[period] = paths
        except (OSError, rasterio.errors.RasterioError) as exc:
            log(f'{period}: recomputing partial correlations from annual data because the cached rasters are unavailable: {exc}')
    return reusable


def resolve_analysis_inputs(ice_path, requested_run=None, rebuild=False):
    """Reuse compatible annual stacks that cover the current analysis domain."""
    runs = []
    if requested_run is not None:
        runs = [Path(requested_run)]
    elif not rebuild and LEGACY_OUT_ROOT.is_dir():
        runs = sorted((path for path in LEGACY_OUT_ROOT.iterdir() if path.is_dir()),
                      reverse=True)
    for run in runs:
        try:
            annual, template, old_domain, config = annual_cache(run)
        except (OSError, ValueError, rasterio.errors.RasterioError) as exc:
            log(f'Skipping cached run {run.name}：{exc}')
            continue
        freq = snow_frequency(template)
        ice = glacier_fraction(ice_path, template)
        domain, counts = analysis_domain(template, freq, ice)
        added = int((domain & ~old_domain).sum())
        del old_domain
        if added:
            log(f'Cached run {run.name} does not cover {added} pixels in the new domain.')
            continue
        partials = reusable_partials(run, template, config)
        log(f'Reusing annual inputs: {annual}')
        return template, {}, annual, domain, counts, partials
    template, inputs = preflight()
    freq = snow_frequency(template)
    ice = glacier_fraction(ice_path, template)
    domain, counts = analysis_domain(template, freq, ice)
    return template, inputs, None, domain, counts, {}


def export_vif(gram, scales, nobs, npix, period, out):
    """Calculate VIFs for the predictors after absorbing pixel and year fixed effects."""
    names = ['snow_increase_cm', 'snow_decrease_cm'] + CONTROLS
    active = np.flatnonzero(scales > 1e-14)
    corr = np.full(gram.shape, np.nan)
    if active.size:
        norm = np.sqrt(np.diag(gram)[active])
        idx = np.ix_(active, active)
        corr[idx] = gram[idx] / np.outer(norm, norm)
    rows = []
    for j, name in enumerate(names):
        tolerance = vif = r2 = np.nan
        flag = 'zero_or_negligible_within_variation'
        if j in active:
            others = active[active != j]
            r = corr[j, others]
            r2 = float(np.clip(r @ np.linalg.pinv(corr[np.ix_(others, others)], hermitian=True) @ r,
                               0, 1)) if others.size else 0.0
            tolerance = 1.0 - r2
            vif = 1.0 / tolerance if tolerance > 1e-12 else np.inf
            flag = ('linear_dependence' if not np.isfinite(vif) else
                    'VIF_ge_10' if vif >= 10 else 'VIF_ge_5' if vif >= 5 else 'VIF_lt_5')
        rows.append({'indicator': INDICATOR, 'Season': period, 'variable': name,
                     'variable_type': 'snow_direction' if j < 2 else 'covariate',
                     'stage': 'within_pixel_year_FE', 'n_obs': nobs, 'n_pixels': npix,
                     'auxiliary_R2': r2, 'tolerance': tolerance, 'VIF': vif, 'status': flag})
    pd.DataFrame(rows).to_csv(out/'FE_covariate_VIF.csv', index=False, encoding='utf-8-sig')


# 6. Period analysis and result tables
def run_period(period,annual_folder,template,domain,reuse_paths=None):
    out=RUN_DIR/period
    out.mkdir()
    panel_dir=out/'panel_blocks';panel_dir.mkdir()
    rasters=out/'partial_rasters';rasters.mkdir()
    years=PARTIAL_YEARS
    sources=[np.load(annual_folder/period/f'{year}.npy',mmap_mode='r') for year in years]
    panel_paths=[]
    # Track valid correlation estimates and testable pixels separately.
    summary={'indicator':INDICATOR,'period':period,'n_domain_pixels':0,'n_r_pixels':0,'n_testable_pixels':0,
             'n_positive':0,'n_negative':0,'n_zero':0,'n_significant_positive':0,'n_significant_negative':0}
    rsum=nesum=0.;nen=0
    with ExitStack() as stack:
        writers = {key: stack.enter_context(rasterio.open(
            rasters/f'{INDICATOR}_{period}_{key}.tif', 'w', **template['profile']))
            for key in PARTIAL_FIELDS}
        readers = {key: stack.enter_context(rasterio.open(path))
                   for key, path in (reuse_paths or {}).items()}
        for r0 in range(0,template['height'],BLOCK_ROWS):
            r1=min(r0+BLOCK_ROWS,template['height'])
            h=r1-r0;w=template['width']
            dom=np.asarray(domain[r0:r1]).reshape(-1)
            window=Window(0,r0,w,h)
            if dom.any():
                cube=np.stack([a[:,r0:r1,:] for a in sources]).reshape(len(years),2+len(CONTROLS),h*w)
                cube[:,:,~dom]=np.nan
                if readers:
                    result = {}
                    for key, src in readers.items():
                        values = src.read(1, window=window, masked=True).filled(np.nan).reshape(-1)
                        result[key] = np.where(dom, values, np.nan)
                else:
                    result=partial_block(cube,years,dom)
                panel=panel_chunk(cube,years,r0,w,template)
                if panel is not None:
                    path=panel_dir/f'rows_{r0:05d}_{r1:05d}.npz'
                    np.savez(path,**panel)
                    panel_paths.append(path)
                r,p=result['r'],result['p']
                rr=np.isfinite(r);pp=rr&np.isfinite(p)
                summary['n_domain_pixels']+=int(dom.sum())
                summary['n_r_pixels']+=int(rr.sum());summary['n_testable_pixels']+=int(pp.sum())
                for key,mask in [('n_positive',rr&(r>0)),('n_negative',rr&(r<0)),('n_zero',rr&(r==0)),
                                 ('n_significant_positive',pp&(p<ALPHA)&(r>0)),('n_significant_negative',pp&(p<ALPHA)&(r<0))]:
                    summary[key]+=int(mask.sum())
                rsum+=float(np.nansum(r))
                ng=rr&np.isfinite(result['n_eff'])
                nesum+=float(result['n_eff'][ng].sum());nen+=int(ng.sum())
            else:
                result={key:np.full(h*w,np.nan) for key in PARTIAL_FIELDS}
            for key,writer in writers.items():
                a=result[key].reshape(h,w)
                writer.write(np.where(np.isfinite(a),a,NODATA).astype('float32'),1,window=window)
            if r0 % (BLOCK_ROWS*20)==0:
                log(f'{period}: processed raster rows {r1}/{template["height"]}')
    sources.clear()
    n=summary['n_r_pixels'];nt=summary['n_testable_pixels']
    summary['partial_source']='reused_rasters' if reuse_paths else 'annual_inputs'
    summary['mean_r_unweighted']=rsum/n if n else np.nan
    summary['mean_n_eff']=nesum/nen if nen else np.nan
    for key in ['positive','negative','zero']:
        summary[f'percent_{key}_among_valid_r']=100*summary['n_'+key]/n if n else np.nan
    for key in ['significant_positive','significant_negative']:
        summary[f'percent_{key}_among_testable']=100*summary['n_'+key]/nt if nt else np.nan
    pd.DataFrame([summary]).to_csv(out/'partial_summary_NH.csv',index=False,encoding='utf-8-sig')
    if n:
        polar_map(rasters/f'{INDICATOR}_{period}_r.tif',template,out/'NH_partial_correlation_map')
    if not panel_paths:
        log(f'{period}: saved partial correlations. No observations met the panel criteria.')
        return summary,None
    model=analyze_panel(panel_paths,period,out)
    return summary,model


def self_test():
    """Check core calculations with synthetic data."""
    raw,used=effective_n_product(25,np.array([.5]),np.array([.4]))
    assert np.allclose(raw,25*(1-.2)/(1+.2))
    raw,_=effective_n_product(25,np.array([.5]),np.array([-.4]))
    assert np.allclose(raw,37.5)
    rng=np.random.default_rng(712)
    t=26;p=4;k=len(CONTROLS)
    z=rng.normal(size=(t,k,p));x=rng.normal(size=(t,p));y=.4*x+z[:,0]+rng.normal(size=(t,p))
    cube=np.concatenate([x[:,None,:],y[:,None,:],z],axis=1)
    result=partial_block(cube,list(range(2000,2000+t)),np.ones(p,bool))
    a=np.column_stack([np.ones(t),z[:,:,0]])
    rx=x[:,0]-a@np.linalg.lstsq(a,x[:,0],rcond=None)[0]
    ry=y[:,0]-a@np.linalg.lstsq(a,y[:,0],rcond=None)[0]
    assert np.allclose(result['r'][0],np.corrcoef(rx,ry)[0,1])
    assert np.isfinite(result['p']).all()
    assert np.isnan(lag1_columns(np.array([[1.,2.,3.,4.]]),[2000,2002,2004,2006])).all()
    print('Self-check passed.')


def aggregate_period_tables(run_dir, periods):
    """Combine result tables across periods; keep only figures/rasters in period folders."""
    names = [
        'partial_summary_NH.csv',
        'FE_coefficients.csv',
        'FE_model_summary.csv',
        'FE_model_diagnostics.csv',
        'FE_covariate_VIF.csv',
        'equal_amplitude_binned_response.csv',
        'equal_amplitude_binned_ratio.csv',
        'equal_amplitude_model_prediction_1_60cm.csv',
    ]
    gathered = {name: [] for name in names}
    covariance_rows = []
    for season in periods:
        folder = run_dir / season
        for name in names:
            path = folder / name
            if path.is_file():
                df = pd.read_csv(path, encoding='utf-8-sig')
                if 'Season' not in df.columns:
                    df.insert(0, 'Season', season)
                gathered[name].append(df)
                path.unlink()
        cov_path = folder / 'FE_covariance_cluster.csv'
        if cov_path.is_file():
            cov = pd.read_csv(cov_path, index_col=0, encoding='utf-8-sig')
            cov.index.name = 'row_variable'
            long_cov = cov.reset_index().melt(id_vars='row_variable',
                var_name='column_variable', value_name='covariance')
            long_cov.insert(0, 'Season', season)
            covariance_rows.append(long_cov)
            cov_path.unlink()
    output_names = {
        'partial_summary_NH.csv': 'NH_all_periods_partial_summary.csv',
        'FE_coefficients.csv': 'FE_coefficients.csv',
        'FE_model_summary.csv': 'NH_all_periods_FE_summary.csv',
        'FE_model_diagnostics.csv': 'FE_model_diagnostics.csv',
        'FE_covariate_VIF.csv': 'FE_covariate_VIF.csv',
        'equal_amplitude_binned_response.csv': 'equal_amplitude_binned_response.csv',
        'equal_amplitude_binned_ratio.csv': 'equal_amplitude_binned_ratio.csv',
        'equal_amplitude_model_prediction_1_60cm.csv': 'equal_amplitude_model_prediction_1_60cm.csv',
    }
    for source_name, frames in gathered.items():
        if frames:
            pd.concat(frames, ignore_index=True, sort=False).to_csv(
                run_dir / output_names[source_name], index=False, encoding='utf-8-sig')
    if covariance_rows:
        pd.concat(covariance_rows, ignore_index=True).to_csv(
            run_dir / 'FE_covariance_cluster.csv', index=False, encoding='utf-8-sig')
    # Retain the familiar summary filename as an alias of the pooled-period table.
    summary_path = run_dir / 'NH_all_periods_FE_summary.csv'
    if summary_path.is_file():
        pd.read_csv(summary_path, encoding='utf-8-sig').to_csv(
            run_dir / 'FE_model_summary.csv', index=False, encoding='utf-8-sig')


def main():
    global RUN_DIR
    parser=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--check-inputs',action='store_true',help='Check inputs and the common spatial mask.')
    parser.add_argument('--self-test',action='store_true',help='Run numerical checks with synthetic data.')
    parser.add_argument('--glacier-mask', type=Path, help='Use a glacier fraction GeoTIFF or NetCDF.')
    cache_group = parser.add_mutually_exclusive_group()
    cache_group.add_argument('--reuse-run', type=Path, help='Reuse a compatible previous run.')
    cache_group.add_argument('--rebuild-annual', action='store_true', help='Rebuild annual data from source rasters.')
    args=parser.parse_args()
    if args.self_test:
        self_test();return
    load_map_boundaries()
    ice_path = glacier_mask_path(args.glacier_mask)
    template, inputs, annual, domain, mask_counts, partials = resolve_analysis_inputs(
        ice_path, args.reuse_run, args.rebuild_annual)
    if args.check_inputs:
        log(f'Input checks passed. The spatial mask retained {mask_counts["n_domain_pixels"]} pixels.')
        return
    RUN_DIR=OUT_ROOT/datetime.now().strftime('%Y%m%d_%H%M%S_%f')
    RUN_DIR.mkdir(parents=True,exist_ok=False)
    log(f'Running {INDICATOR}. Output: {RUN_DIR}')
    config={'indicator':INDICATOR,'partial_years':PARTIAL_YEARS,'panel_years':PANEL_YEARS,
            'periods':PERIODS,'controls':CONTROLS,'effective_n_formula':'n*(1-phi_x*phi_y)/(1+phi_x*phi_y)',
            'cap_neff_at_n':CAP_NEFF_AT_N,'min_adjacent_pairs':MIN_ADJACENT_PAIRS,
            'snow_frequency_threshold':SNOW_FREQUENCY_THRESHOLD,
            'glacier_fraction_exclusion_threshold':GLACIER_FRACTION_THRESHOLD,
            'glacier_mask_path':str(ice_path.resolve()),
            'annual_cache_reused_from':str(annual.resolve()) if annual is not None else None,
            'partial_rasters_reused_from':{
                period:str(paths['r'].parent.resolve()) for period,paths in partials.items()},
            'vif_stage':'within_pixel_year_FE',
            'analysis_min_lat':ANALYSIS_MIN_LAT,'map_min_lat':MAP_MIN_LAT,
            'require_complete_partial':REQUIRE_COMPLETE_PARTIAL,'min_partial_years':MIN_PARTIAL_YEARS,
            'min_panel_observations':MIN_PANEL_OBSERVATIONS,'bin_edges_cm':AMP_EDGES_CM.tolist(),
            'bin_centers_cm':AMP_CENTERS_CM.tolist(),'min_bin_n':MIN_BIN_N,
            'adjustment_definition':ADJUSTED_RESPONSE_DEFINITION,'figsize_inches':FIGSIZE,
            'map_axes_rectangle':MAP_AX_RECT,'dpi':FIG_DPI,
            'map_background_color':MAP_BACKGROUND_COLOR,'map_grid_color':MAP_GRID_COLOR,
            'map_country_color':MAP_COUNTRY_COLOR,'map_border_color':MAP_BORDER_COLOR,
            'map_parallels':MAP_PARALLELS,'map_meridians':MAP_MERIDIANS,
            'map_colorbar':False,'map_external_coordinate_labels':False,
            'negative_color':NEGATIVE_COLOR,'positive_color':POSITIVE_COLOR,
            'cr1_count_absorbed_fe':CR1_COUNT_ABSORBED_FE}
    (RUN_DIR/'run_config.json').write_text(json.dumps(config,ensure_ascii=False,indent=2),encoding='utf-8')
    try:
        write_tif(RUN_DIR/'analysis_domain.tif', domain.astype(np.float32), template)
        pd.DataFrame([mask_counts]).to_csv(RUN_DIR/'analysis_domain_QA.csv',
                                          index=False, encoding='utf-8-sig')
        if annual is None:
            annual=prepare_annual(template,domain,inputs)
        summaries=[];models=[]
        for period in PERIODS:
            log(f'Processing {period}. Covariates: {CONTROLS}')
            s,m=run_period(period,annual,template,domain,partials.get(period))
            summaries.append(s)
            if m is not None:
                models.append(m)
        pd.DataFrame(summaries).to_csv(RUN_DIR/'NH_all_periods_partial_summary.csv',index=False,encoding='utf-8-sig')
        pd.DataFrame(models).to_csv(RUN_DIR/'NH_all_periods_FE_summary.csv',index=False,encoding='utf-8-sig')
        aggregate_period_tables(RUN_DIR, list(PERIODS))
        log('Analysis complete. See FE_covariate_VIF.csv for model predictor VIFs.')
    except Exception:
        aggregate_period_tables(RUN_DIR, list(PERIODS))
        raise
    finally:
        (RUN_DIR/'input_manifest.json').write_text(json.dumps(list(SOURCE_RECORDS.values()),ensure_ascii=False,indent=2),encoding='utf-8')


if __name__=='__main__':
    try:
        main()
    except Exception as exc:
        log(f'Run failed: {type(exc).__name__}: {exc}')
        raise
