# ============================================================
# Snowpack Meta-analysis: Funnel Plots and Subgroup Statistics
# ============================================================

# 1. Packages
packages <- c("readxl", "metafor", "ggplot2", "dplyr", "openxlsx")
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}
suppressPackageStartupMessages({
  library(readxl)
  library(metafor)
  library(ggplot2)
  library(dplyr)
  library(openxlsx)
  library(grid)
})

# 2. Paths and settings
input_file <- "E:/Snowpack_Meta_Analysis/Input/Meta_effect_sizes_2026-10-29_Nakagawa_Hamman_QC.xlsx"
output_dir <- "E:/Snowpack_Meta_Analysis/Output/Funnel_Plots"

single_dir <- file.path(output_dir, "01_Single_Funnel_Plots")
row_dir    <- file.path(output_dir, "02_Row_Combination_Plots")
full_dir   <- file.path(output_dir, "03_Full_Combination_Plot")
table_dir  <- file.path(output_dir, "04_Statistics_and_QC")
for (d in c(single_dir, row_dir, full_dir, table_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# Source column names
COL_TRAIT <- "Trait"
COL_TREATMENT <- "Treatment_1"
COL_YI <- "yi"
COL_VI <- "vi"
COL_STUDY <- "Study_id"
GROUP_COLS <- c("Eco_1", "Types_of_frozen_soil", "Sampling_season",
                "StandAge_year", "Soil_depth")

MIN_N_FOR_FUNNEL <- 10L
DPI <- 600
FONT <- "Arial"
BASE_SIZE <- 11
AXIS_TITLE_SIZE <- BASE_SIZE + 2

# Colors
COLOR_ADDITION <- "#DF7E21"
COLOR_REMOVAL  <- "#5A769E"
POINT_BORDER   <- "gray55"

# Single-panel dimensions
SINGLE_WIDTH_CM  <- 7.0
SINGLE_HEIGHT_CM <- 7.0

# Combination-panel dimensions
COMBO_PANEL_WIDTH_CM  <- 6.3
COMBO_PANEL_HEIGHT_CM <- 6.2
GROUP_GAP_WIDTH_RATIO <- 0.50 / 3

# 3. Helper functions
clean_filename <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("₂", "2", x, fixed = TRUE)
  x <- gsub("₄", "4", x, fixed = TRUE)
  x <- gsub("β", "beta", x, fixed = TRUE)
  x <- iconv(x, from = "", to = "ASCII//TRANSLIT", sub = "_")
  x <- gsub('[\\/:*?"<>|]', "_", x)
  x <- gsub("[[:cntrl:]]", "_", x)
  x
}

normalize_treatment <- function(x) {
  s <- tolower(trimws(as.character(x)))
  out <- rep(NA_character_, length(s))
  out[grepl("addition|increas|snow[._ -]?(add|plus)|thicken|deepen", s)] <- "Addition"
  out[grepl("remov|decreas|reduc|snow[._ -]?(minus|excl|remov)|thin", s)] <- "Removal"
  out
}

normalize_trait <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("[[:space:]_–—-]+", "", x)
  tolower(x)
}

match_trait <- function(all_traits, aliases) {
  hit <- all_traits[normalize_trait(all_traits) %in% normalize_trait(aliases)]
  unique(hit)
}

# Plot theme
apply_funnel_theme <- function(p, show_title = FALSE, show_subtitle = FALSE,
                               title_hjust = 0.5) {
  p +
    theme_classic(base_family = FONT, base_size = BASE_SIZE) +
    theme(
      text = element_text(family = FONT, colour = "black"),
      plot.title = if (show_title) {
        element_text(size = BASE_SIZE + 1, face = "bold", hjust = title_hjust)
      } else {
        element_blank()
      },
      plot.subtitle = if (show_subtitle) {
        element_text(size = BASE_SIZE, hjust = 0.5)
      } else {
        element_blank()
      },
      axis.title = element_text(size = AXIS_TITLE_SIZE),
      axis.text = element_text(size = BASE_SIZE - 1),
      axis.line = element_blank(),
      axis.ticks = element_line(linewidth = 0.4),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
      aspect.ratio = 1,
      plot.margin = margin(6, 6, 6, 6, "pt")
    )
}

make_placeholder <- function(msg, title = NULL, show_title = FALSE, title_hjust = 0.5) {
  p <- ggplot() +
    annotate("text", x = 0.5, y = 0.55, label = msg,
             size = 3.5, family = FONT, color = "gray35") +
    xlim(0, 1) + ylim(0, 1) +
    labs(title = if (show_title) title else NULL,
         x = "Effect size (Hedges’ g)", y = "Standard error (SE)") +
    theme_void(base_family = FONT, base_size = BASE_SIZE)

  p +
    theme(
      text = element_text(family = FONT, colour = "black"),
      plot.title = if (show_title) {
        element_text(size = BASE_SIZE + 1, face = "bold", hjust = title_hjust)
      } else {
        element_blank()
      },
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
      aspect.ratio = 1,
      plot.margin = margin(6, 6, 6, 6, "pt")
    )
}

make_funnel <- function(data, title = NULL, trt, show_title = FALSE,
                        show_subtitle = FALSE, title_hjust = 0.5) {
  n_effects <- nrow(data)
  n_studies <- if (COL_STUDY %in% names(data)) {
    ids <- trimws(as.character(data[[COL_STUDY]]))
    length(unique(ids[!is.na(ids) & nzchar(ids)]))
  } else {
    NA_integer_
  }

  result <- list(
    plot = NULL, n_effects = n_effects, n_studies = n_studies,
    pooled_effect = NA_real_, egger_p_exploratory = NA_real_, status = ""
  )

  if (n_effects < MIN_N_FOR_FUNNEL) {
    result$plot <- make_placeholder(
      msg = paste0("n = ", n_effects, " (< ", MIN_N_FOR_FUNNEL, ")"),
      title = title, show_title = show_title, title_hjust = title_hjust
    )
    result$status <- "Too few effects"
    return(result)
  }

  model <- tryCatch(
    metafor::rma(yi = yi, vi = vi, data = data, method = "REML"),
    error = function(e) {
      message("REML model failed: ", title %||% "", " / ", trt, " / ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(model)) {
    result$plot <- make_placeholder(
      msg = "REML fitting failed",
      title = title, show_title = show_title, title_hjust = title_hjust
    )
    result$status <- "Model failed"
    return(result)
  }

  pooled <- as.numeric(coef(model)[1])
  result$pooled_effect <- pooled
  result$status <- "OK"

  if (is.na(n_studies) || n_studies >= 10L) {
    egger <- tryCatch(
      metafor::regtest(model, model = "lm"),
      error = function(e) NULL
    )
    if (!is.null(egger)) result$egger_p_exploratory <- as.numeric(egger$pval)
  }

  df_plot <- data.frame(yi = data$yi, SE = sqrt(data$vi))
  se_max <- max(df_plot$SE, na.rm = TRUE) * 1.08
  se_seq <- seq(0, se_max, length.out = 200L)
  edges <- data.frame(
    SE = rep(se_seq, times = 2L),
    x = c(pooled + 1.96 * se_seq, pooled - 1.96 * se_seq),
    side = rep(c("Right", "Left"), each = length(se_seq))
  )

  xrange <- range(c(df_plot$yi, edges$x, pooled), finite = TRUE)
  padx <- diff(xrange) * 0.06
  if (!is.finite(padx) || padx <= 0) padx <- max(abs(pooled) * 0.10, 0.1)
  xrange <- xrange + c(-padx, padx)

  fill_col <- if (trt == "Addition") COLOR_ADDITION else COLOR_REMOVAL

  p <- ggplot(df_plot, aes(x = yi, y = SE)) +
    geom_vline(xintercept = pooled, color = "gray30",
               linetype = "dashed", linewidth = 0.50) +
    geom_line(data = edges, aes(x = x, y = SE, group = side),
              inherit.aes = FALSE, color = "gray30",
              linetype = "dashed", linewidth = 0.60) +
    geom_point(shape = 21, fill = fill_col, colour = POINT_BORDER,
               stroke = 0.45, size = 2.15, alpha = 0.85) +
    scale_y_reverse(limits = c(se_max, 0), expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = xrange, expand = FALSE) +
    labs(
      title = if (show_title) title else NULL,
      subtitle = NULL,
      x = "Effect size (Hedges’ g)",
      y = "Standard error (SE)"
    )

  result$plot <- apply_funnel_theme(
    p, show_title = show_title, show_subtitle = show_subtitle,
    title_hjust = title_hjust
  )
  result
}

save_tiff <- function(p, target, width_cm, height_cm) {
  ggplot2::ggsave(
    filename = target, plot = p,
    width = width_cm, height = height_cm, units = "cm",
    dpi = DPI, device = "tiff", compression = "lzw",
    limitsize = FALSE, bg = "white"
  )
}

`%||%` <- function(a, b) if (!is.null(a)) a else b


# 3.1 Subgroup statistics
# Map study duration to predefined intervals; unmatched values remain NA.
AGE_LEVELS <- c("≤1 year", "2-5 years", "6-10 years", "11-15 years", "16-25 years")
DEPTH_LEVELS <- c("0-10 cm", "11-30 cm", "31-60 cm", "61-100 cm")

extract_numbers <- function(x) {
  z <- as.character(x)
  m <- gregexpr("[0-9]+(?:\\.[0-9]+)?", z, perl = TRUE)
  suppressWarnings(as.numeric(unlist(regmatches(z, m))))
}

bin_stand_age <- function(x) {
  x_chr <- trimws(as.character(x))
  out <- rep(NA_character_, length(x_chr))
  for (i in seq_along(x_chr)) {
    if (is.na(x_chr[i]) || !nzchar(x_chr[i])) next
    nums <- extract_numbers(x_chr[i])
    nums <- nums[is.finite(nums)]
    if (length(nums) == 1L) {
      v <- nums[1]
      if (v <= 1) out[i] <- AGE_LEVELS[1]
      else if (v >= 2 && v <= 5) out[i] <- AGE_LEVELS[2]
      else if (v >= 6 && v <= 10) out[i] <- AGE_LEVELS[3]
      else if (v >= 11 && v <= 15) out[i] <- AGE_LEVELS[4]
      else if (v >= 16 && v <= 25) out[i] <- AGE_LEVELS[5]
    } else if (length(nums) >= 2L) {
      lo <- min(nums[1:2]); hi <- max(nums[1:2])
      if (lo >= 0 && hi <= 1) out[i] <- AGE_LEVELS[1]
      else if (lo >= 2 && hi <= 5) out[i] <- AGE_LEVELS[2]
      else if (lo >= 6 && hi <= 10) out[i] <- AGE_LEVELS[3]
      else if (lo >= 11 && hi <= 15) out[i] <- AGE_LEVELS[4]
      else if (lo >= 16 && hi <= 25) out[i] <- AGE_LEVELS[5]
    }
  }
  out
}

bin_soil_depth <- function(x) {
  x_chr <- trimws(as.character(x))
  out <- rep(NA_character_, length(x_chr))
  for (i in seq_along(x_chr)) {
    if (is.na(x_chr[i]) || !nzchar(x_chr[i])) next
    nums <- extract_numbers(x_chr[i])
    nums <- nums[is.finite(nums)]
    if (length(nums) == 1L) {
      v <- nums[1]
      if (v >= 0 && v <= 10) out[i] <- DEPTH_LEVELS[1]
      else if (v >= 11 && v <= 30) out[i] <- DEPTH_LEVELS[2]
      else if (v >= 31 && v <= 60) out[i] <- DEPTH_LEVELS[3]
      else if (v >= 61 && v <= 100) out[i] <- DEPTH_LEVELS[4]
    } else if (length(nums) >= 2L) {
      lo <- min(nums[1:2]); hi <- max(nums[1:2])
      if (lo >= 0 && hi <= 10) out[i] <- DEPTH_LEVELS[1]
      else if (lo >= 11 && hi <= 30) out[i] <- DEPTH_LEVELS[2]
      else if (lo >= 31 && hi <= 60) out[i] <- DEPTH_LEVELS[3]
      else if (lo >= 61 && hi <= 100) out[i] <- DEPTH_LEVELS[4]
    }
  }
  out
}

count_studies <- function(dat) {
  if (!COL_STUDY %in% names(dat)) return(NA_integer_)
  ids <- trimws(as.character(dat[[COL_STUDY]]))
  ids <- ids[!is.na(ids) & nzchar(ids)]
  length(unique(ids))
}

calculate_one_meta_row <- function(dat, group_type, group_value, trait, trt) {
  n_eff <- nrow(dat)
  study_n <- count_studies(dat)
  mu <- se <- ci_low <- ci_high <- tau2 <- i2 <- egger_stat <- egger_p <- NA_real_

  # Fit the REML model when at least two effect sizes are available.
  model <- NULL
  if (n_eff >= 2L) {
    model <- tryCatch(
      metafor::rma(yi = yi, vi = vi, data = dat, method = "REML"),
      error = function(e) {
        message("Subgroup model failed: ", group_type, " / ", group_value,
                " / ", trait, " / ", trt, "; ", conditionMessage(e))
        NULL
      }
    )
  }

  if (!is.null(model)) {
    mu <- as.numeric(coef(model)[1])
    se <- as.numeric(model$se[1])
    ci_low <- as.numeric(model$ci.lb[1])
    ci_high <- as.numeric(model$ci.ub[1])
    tau2 <- as.numeric(model$tau2[1])
    i2 <- as.numeric(model$I2[1])

    # Apply the same Egger regression used for funnel plots.
    # Require at least 10 effect sizes and 10 independent studies.
    enough_studies <- is.na(study_n) || study_n >= 10L
    if (n_eff >= 10L && enough_studies) {
      egger <- tryCatch(
        metafor::regtest(model, model = "lm"),
        error = function(e) NULL
      )
      if (!is.null(egger)) {
        if (!is.null(egger$zval)) egger_stat <- as.numeric(egger$zval[1])
        if (!is.null(egger$pval)) egger_p <- as.numeric(egger$pval[1])
      }
    }
  }

  data.frame(
    Group_type = group_type,
    Group_value = group_value,
    Trait = trait,
    Treatment_1 = if (trt == "Addition") "Increased snowpack thickness" else "Decreased snowpack thickness",
    n_effects = as.integer(n_eff),
    Study_n = as.integer(study_n),
    REML_mu = mu,
    REML_SE = se,
    REML_CI_low = ci_low,
    REML_CI_high = ci_high,
    tau2 = tau2,
    I2_percent = i2,
    Egger_statistic = egger_stat,
    Egger_p = egger_p,
    stringsAsFactors = FALSE
  )
}

build_all_group_statistics <- function(dat) {
  dat_stats <- dat
  dat_stats$StandAge_group <- bin_stand_age(dat_stats$StandAge_year)
  dat_stats$Soil_depth_group <- bin_soil_depth(dat_stats$Soil_depth)

  get_observed_values <- function(x) {
    vals <- trimws(as.character(x))
    unique(vals[!is.na(vals) & nzchar(vals)])
  }

  group_specs <- list(
    list(type = "Overall", column = NULL, values = "Overall"),
    list(type = "Eco_1", column = "Eco_1", values = get_observed_values(dat_stats$Eco_1)),
    list(type = "Types_of_frozen_soil", column = "Types_of_frozen_soil", values = get_observed_values(dat_stats$Types_of_frozen_soil)),
    list(type = "Sampling_season", column = "Sampling_season", values = get_observed_values(dat_stats$Sampling_season)),
    list(type = "StandAge_year", column = "StandAge_group", values = AGE_LEVELS),
    list(type = "Soil_depth", column = "Soil_depth_group", values = DEPTH_LEVELS)
  )

  traits <- sort(unique(dat_stats$Trait))
  treatments <- c("Addition", "Removal")
  out_rows <- list()
  row_id <- 0L

  for (spec in group_specs) {
    if (length(spec$values) == 0L) next
    for (group_value in spec$values) {
      if (is.null(spec$column)) {
        group_dat <- dat_stats
      } else {
        gv <- trimws(as.character(dat_stats[[spec$column]]))
        group_dat <- dat_stats[!is.na(gv) & gv == group_value, , drop = FALSE]
      }

      for (trait in traits) {
        for (trt in treatments) {
          one <- group_dat[group_dat$Trait == trait & group_dat$Treatment == trt, , drop = FALSE]
          row_id <- row_id + 1L
          out_rows[[row_id]] <- calculate_one_meta_row(
            dat = one,
            group_type = spec$type,
            group_value = group_value,
            trait = trait,
            trt = trt
          )
        }
      }
    }
  }

  if (length(out_rows) == 0L) return(data.frame())
  do.call(rbind, out_rows)
}

# 4. Read and validate data
if (!file.exists(input_file)) stop("Input file not found: ", input_file)
raw <- readxl::read_excel(input_file)
names(raw) <- trimws(names(raw))
required <- c(COL_TRAIT, COL_TREATMENT, COL_YI, COL_VI, GROUP_COLS)
missing_cols <- setdiff(required, names(raw))
if (length(missing_cols) > 0L) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "),
       "\nAvailable columns: ", paste(names(raw), collapse = " | "))
}

use_cols <- unique(c(required, intersect(COL_STUDY, names(raw))))
df <- as.data.frame(raw[, use_cols, drop = FALSE])
names(df)[names(df) == COL_TRAIT] <- "Trait"
names(df)[names(df) == COL_TREATMENT] <- "Treatment_raw"
names(df)[names(df) == COL_YI] <- "yi"
names(df)[names(df) == COL_VI] <- "vi"

df$Trait <- trimws(as.character(df$Trait))
df$Treatment_raw <- trimws(as.character(df$Treatment_raw))
df$Treatment <- normalize_treatment(df$Treatment_raw)
df$yi <- suppressWarnings(as.numeric(df$yi))
df$vi <- suppressWarnings(as.numeric(df$vi))

unknown <- sort(unique(df$Treatment_raw[!is.na(df$Treatment_raw) &
                                           nzchar(df$Treatment_raw) &
                                           is.na(df$Treatment)]))
if (length(unknown) > 0L) {
  stop("Unrecognized Treatment_1 values. Update normalize_treatment():\n",
       paste(unknown, collapse = " | "))
}

n_before <- nrow(df)
df <- df[!is.na(df$Trait) & nzchar(df$Trait) &
           !is.na(df$Treatment) & is.finite(df$yi) &
           is.finite(df$vi) & df$vi > 0, , drop = FALSE]
message("Input rows: ", n_before,
        "; valid rows: ", nrow(df),
        "; excluded rows: ", n_before - nrow(df))
if (nrow(df) == 0L) stop("No valid effect-size records are available.")

# 5. Save single funnel plots for all traits and treatments
all_traits <- sort(unique(df$Trait))
treatment_order <- c("Addition", "Removal")
plot_store <- list()
summary_rows <- list()
k <- 0L

for (trait in all_traits) {
  for (trt in treatment_order) {
    data_one <- df[df$Trait == trait & df$Treatment == trt, , drop = FALSE]
    res <- make_funnel(
      data = data_one, title = NULL, trt = trt,
      show_title = FALSE, show_subtitle = FALSE
    )
    key <- paste(trait, trt, sep = "|||")
    plot_store[[key]] <- res$plot

    out <- file.path(single_dir, paste0(clean_filename(trait), "_",
                                        if (trt == "Addition") "Increased_snowpack" else "Decreased_snowpack",
                                        "_funnel.tif"))
    if (nrow(data_one) > 0L) {
      save_tiff(res$plot, out, SINGLE_WIDTH_CM, SINGLE_HEIGHT_CM)
    } else {
      out <- NA_character_
    }

    k <- k + 1L
    summary_rows[[k]] <- data.frame(
      Trait = trait,
      Treatment = if (trt == "Addition") "Increased snowpack thickness" else "Decreased snowpack thickness",
      n_effects = res$n_effects,
      Study_n = res$n_studies,
      Pooled_effect_REML = res$pooled_effect,
      Egger_p_exploratory = res$egger_p_exploratory,
      Status = res$status,
      File = out,
      stringsAsFactors = FALSE
    )
  }
}
summary_df <- do.call(rbind, summary_rows)
write.csv(summary_df,
          file.path(table_dir, "All_traits_funnel_Egger_exploratory_results.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

# 6. Define the 15 traits used in combination plots
combo_traits <- list(
  "Soil temperature" = c("ST", "Soil temperature"),
  "Soil water content" = c("SWC", "Soil water content", "Soil moisture"),
  "Snowmelt date" = c("thaw-date", "Thaw date", "Thawing date", "Snow-off", "Snowmelt date"),
  "Number of freeze-thaw cycles" = c("FTC", "Freeze-thaw cycles", "Freeze-thaw cycle", "Number of freeze-thaw cycles"),
  "Active layer thickness" = c("Active layer thickness depth", "Active layer thickness", "ALT"),
  "CO2 flux" = c("CO2", "CO₂", "CO2 flux", "CO₂ flux"),
  "CH4 flux" = c("CH4", "CH₄", "CH4 flux", "CH₄ flux"),
  "Soil organic carbon" = c("SOC", "Soil organic carbon"),
  "β-glucosidase" = c("βG", "betaG", "BG", "β-glucosidase", "beta-glucosidase"),
  "β-n-acetylglucosaminidase" = c("NAG", "N-acetylglucosaminidase", "N-acetyl-β-glucosaminidase", "β-n-acetylglucosaminidase"),
  "Invertase" = c("Invertase"),
  "Aboveground biomass" = c("AGB", "Aboveground biomass"),
  "Belowground biomass" = c("BGB", "Belowground biomass"),
  "Plant cover (%)" = c("Plant cover", "Plant coverage", "Vegetation cover", "Vegetation coverage"),
  "Litter decomposition rate" = c("Litter decomposition rate", "Litter decomposition")
)
stopifnot(length(combo_traits) == 15L)

matching_rows <- list()
for (i in seq_along(combo_traits)) {
  hits <- match_trait(all_traits, combo_traits[[i]])
  matching_rows[[i]] <- data.frame(
    Order = i,
    Display = names(combo_traits)[i],
    Matched_Trait = if (length(hits) > 0L) hits[1] else NA_character_,
    Match_count = length(hits),
    All_matches = paste(hits, collapse = " | "),
    Status = if (length(hits) == 0L) "NOT FOUND" else
      if (length(hits) > 1L) "MULTIPLE: CHECK" else "OK",
    stringsAsFactors = FALSE
  )
}
matching_df <- do.call(rbind, matching_rows)
write.csv(matching_df,
          file.path(table_dir, "Combination_plot_trait_matching.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

if (any(matching_df$Status != "OK")) {
  message("Some combination-plot traits were not matched or matched multiple source traits. Check Combination_plot_trait_matching.csv.")
  print(matching_df[matching_df$Status != "OK", ], row.names = FALSE)
}


# 6.1 Export overall and subgroup statistics for all traits
# Include every trait in the source data, not only the 15 plotted traits.
# Group_value stores the observed category within each grouping variable.
all_stats_df <- build_all_group_statistics(df)

if (nrow(all_stats_df) == 0L) {
  warning("No summary statistics were generated. Check grouping columns and effect-size data.")
} else {
  stats_file <- file.path(table_dir, "All_traits_overall_and_subgroup_REML_Egger_summary.xlsx")

  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "All_statistics")
  openxlsx::writeDataTable(
    wb, sheet = "All_statistics", x = all_stats_df,
    tableStyle = "TableStyleMedium2", withFilter = TRUE
  )
  openxlsx::freezePane(wb, "All_statistics", firstRow = TRUE)
  openxlsx::setColWidths(wb, "All_statistics", cols = seq_len(ncol(all_stats_df)), widths = "auto")

  openxlsx::addWorksheet(wb, "Overall_only")
  overall_df <- all_stats_df[all_stats_df$Group_type == "Overall", , drop = FALSE]
  openxlsx::writeDataTable(
    wb, sheet = "Overall_only", x = overall_df,
    tableStyle = "TableStyleMedium2", withFilter = TRUE
  )
  openxlsx::freezePane(wb, "Overall_only", firstRow = TRUE)
  openxlsx::setColWidths(wb, "Overall_only", cols = seq_len(ncol(overall_df)), widths = "auto")

  openxlsx::addWorksheet(wb, "Grouping_check")
  check_df <- data.frame(
    Group_type = c("Eco_1", "Types_of_frozen_soil", "Sampling_season", "StandAge_year", "Soil_depth"),
    Source_column = c("Eco_1", "Types_of_frozen_soil", "Sampling_season", "StandAge_year", "Soil_depth"),
    Group_values_used = c(
      paste(unique(na.omit(trimws(as.character(df$Eco_1)))), collapse = " | "),
      paste(unique(na.omit(trimws(as.character(df$Types_of_frozen_soil)))), collapse = " | "),
      paste(unique(na.omit(trimws(as.character(df$Sampling_season)))), collapse = " | "),
      paste(unique(na.omit(bin_stand_age(df$StandAge_year))), collapse = " | "),
      paste(unique(na.omit(bin_soil_depth(df$Soil_depth))), collapse = " | ")
    ),
    stringsAsFactors = FALSE
  )
  openxlsx::writeDataTable(wb, "Grouping_check", check_df,
                           tableStyle = "TableStyleMedium2", withFilter = TRUE)
  openxlsx::setColWidths(wb, "Grouping_check", cols = 1:3, widths = "auto")

  openxlsx::addWorksheet(wb, "Plot_trait_matching")
  openxlsx::writeDataTable(wb, "Plot_trait_matching", matching_df,
                           tableStyle = "TableStyleMedium2", withFilter = TRUE)
  openxlsx::setColWidths(wb, "Plot_trait_matching", cols = 1:ncol(matching_df), widths = "auto")

  openxlsx::saveWorkbook(wb, stats_file, overwrite = TRUE)
  message("All-trait subgroup statistics saved: ", stats_file)
}

# 7. Combination-plot panels
# Draw trait titles separately above each addition/removal pair.
# This allows long titles to span the full pair width.
# Align the title with the left edge of the addition panel.
make_combo_single <- function(actual_trait, trt) {
  if (is.na(actual_trait)) {
    p <- make_placeholder("Trait not found")
  } else {
    key <- paste(actual_trait, trt, sep = "|||")
    p <- plot_store[[key]]
    if (is.null(p)) p <- make_placeholder("No data")
  }

  p +
    labs(title = NULL, subtitle = NULL) +
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.margin = margin(t = 2, r = 6, b = 6, l = 6, unit = "pt")
    )
}

# 8. Draw combination plots with grid
# Each trait has one title row spanning its addition and removal panels.
# The title is left-aligned with the addition panel and is not clipped.
# Each row contains three trait pairs separated by group spacers.
# Fifteen traits produce five rows and 30 funnel panels.
GROUPS_PER_ROW <- 3L
N_ROWS <- ceiling(nrow(matching_df) / GROUPS_PER_ROW)
TITLE_HEIGHT_CM <- 0.55
TITLE_X_IN_PAIR <- 0.105
ROW_WIDTHS <- c(1, 1, GROUP_GAP_WIDTH_RATIO,
                1, 1, GROUP_GAP_WIDTH_RATIO,
                1, 1)
GROUP_START_COL <- c(1L, 4L, 7L)

# Draw a centered legend with treatment colors.
# The legend is drawn directly in the grid viewport.
draw_grid_legend <- function() {
  grid::grid.circle(
    x = grid::unit(0.225, "npc"), y = grid::unit(0.50, "npc"),
    r = grid::unit(1.4, "mm"),
    gp = grid::gpar(fill = COLOR_ADDITION, col = POINT_BORDER, lwd = 0.7)
  )
  grid::grid.text(
    "Increased snowpack thickness",
    x = grid::unit(0.225, "npc") + grid::unit(3.5, "mm"),
    y = grid::unit(0.50, "npc"), just = "left",
    gp = grid::gpar(fontfamily = FONT, fontsize = AXIS_TITLE_SIZE, col = "black")
  )
  grid::grid.circle(
    x = grid::unit(0.555, "npc"), y = grid::unit(0.50, "npc"),
    r = grid::unit(1.4, "mm"),
    gp = grid::gpar(fill = COLOR_REMOVAL, col = POINT_BORDER, lwd = 0.7)
  )
  grid::grid.text(
    "Decreased snowpack thickness",
    x = grid::unit(0.555, "npc") + grid::unit(3.5, "mm"),
    y = grid::unit(0.50, "npc"), just = "left",
    gp = grid::gpar(fontfamily = FONT, fontsize = AXIS_TITLE_SIZE, col = "black")
  )
}

# Use plotmath for CO2 and CH4 subscripts.
get_display_title_label <- function(display_title) {
  if (identical(display_title, "CO2 flux")) {
    return(expression(bold(CO[2]~flux)))
  }
  if (identical(display_title, "CH4 flux")) {
    return(expression(bold(CH[4]~flux)))
  }
  return(display_title)
}

# Draw each panel directly into the TIFF grid layout.
# indices selects one row of traits or all traits.
save_combo_grid_tiff <- function(indices, outfile) {
  stopifnot(length(indices) > 0L)
  rows_here <- ceiling(length(indices) / GROUPS_PER_ROW)
  if (rows_here < 1L) stop("No traits were supplied for the combination plot.")

  legend_height_cm <- 1.15
  width_cm <- COMBO_PANEL_WIDTH_CM * sum(ROW_WIDTHS)
  # Allocate a title row across each addition/removal pair.
  height_cm <- (COMBO_PANEL_HEIGHT_CM + TITLE_HEIGHT_CM) * rows_here + legend_height_cm
  layout_heights <- c(rep(c(TITLE_HEIGHT_CM, COMBO_PANEL_HEIGHT_CM), rows_here),
                      legend_height_cm)

  grDevices::tiff(
    filename = outfile,
    width = width_cm, height = height_cm,
    units = "cm", res = DPI, compression = "lzw",
    bg = "white"
  )
  # Always close the graphics device, including after plotting errors.
  tryCatch({
    grid::grid.newpage()
    outer_vp <- grid::viewport(
      layout = grid::grid.layout(
        nrow = rows_here * 2L + 1L,
        ncol = length(ROW_WIDTHS),
        widths = grid::unit(ROW_WIDTHS, "null"),
        heights = grid::unit(layout_heights, "cm")
      )
    )
    grid::pushViewport(outer_vp)

    for (j in seq_along(indices)) {
      actual_index <- indices[j]
      actual_trait <- matching_df$Matched_Trait[actual_index]
      display_title <- matching_df$Display[actual_index]
      display_label <- get_display_title_label(display_title)
      row_pos <- (j - 1L) %/% GROUPS_PER_ROW + 1L
      group_pos <- (j - 1L) %% GROUPS_PER_ROW + 1L
      col_add <- GROUP_START_COL[group_pos]
      title_row <- (row_pos - 1L) * 2L + 1L
      plot_row <- title_row + 1L
      col_end <- if (group_pos < GROUPS_PER_ROW) col_add + 2L else col_add + 1L

      # Draw the title across the full treatment pair.
      # TITLE_X_IN_PAIR aligns the title with the addition-panel plotting region.
      grid::pushViewport(
        grid::viewport(layout.pos.row = title_row,
                       layout.pos.col = c(col_add, col_end))
      )
      grid::grid.text(
        label = display_label,
        x = grid::unit(TITLE_X_IN_PAIR, "npc"),
        y = grid::unit(0.50, "npc"),
        just = c("left", "centre"),
        gp = grid::gpar(
          fontfamily = FONT, fontsize = BASE_SIZE + 1,
          fontface = "bold", col = "black"
        )
      )
      grid::upViewport()

      # Keep panel titles empty; the pair title is drawn above both panels.
      for (trt in c("Addition", "Removal")) {
        col_pos <- col_add + if (trt == "Removal") 1L else 0L
        one_plot <- make_combo_single(actual_trait, trt)
        grid::pushViewport(
          grid::viewport(layout.pos.row = plot_row, layout.pos.col = col_pos)
        )
        grid::grid.draw(ggplot2::ggplotGrob(one_plot))
        grid::upViewport()
      }
    }

    # Draw one legend below each combination figure.
    grid::pushViewport(
      grid::viewport(layout.pos.row = rows_here * 2L + 1L,
                     layout.pos.col = c(1L, length(ROW_WIDTHS)))
    )
    draw_grid_legend()
    grid::upViewport(2L)
  }, finally = {
    grDevices::dev.off()
  })

  message("Combination plot saved: ", outfile)
  invisible(outfile)
}

# 9. Save one combination figure per row
for (row_i in seq_len(N_ROWS)) {
  idx <- seq.int(from = (row_i - 1L) * GROUPS_PER_ROW + 1L,
                 to = min(row_i * GROUPS_PER_ROW, nrow(matching_df)))
  row_file <- file.path(row_dir,
                        sprintf("Funnel_combination_row_%02d.tif", row_i))
  save_combo_grid_tiff(idx, row_file)
}

# 10. Save the complete 15-trait combination figure
full_file <- file.path(full_dir, "Funnel_combination_all_15_traits.tif")
save_combo_grid_tiff(seq_len(nrow(matching_df)), full_file)

cat("\nAnalysis completed.\n")
cat("Single funnel plots: ", single_dir, "\n", sep = "")
cat("Row combination plots: ", row_dir, "\n", sep = "")
cat("Full combination plot: ", full_file, "\n", sep = "")
cat("Statistics and QC: ", table_dir, "\n", sep = "")
cat("All-trait subgroup statistics: ",
    file.path(table_dir, "All_traits_overall_and_subgroup_REML_Egger_summary.xlsx"),
    "\n", sep = "")
cat("Egger results are exploratory and do not account for within-study effect-size dependence.\n")
