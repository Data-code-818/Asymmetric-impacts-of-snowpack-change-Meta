# =============================================================================
# Linear mixed-effects regression among paired ecological response indicators
# =============================================================================
# Effect sizes: Hedges' g (yi) from the processed meta-analysis dataset.
# Pairing identifier: Pair_ID.
# Model: Y ~ X + (1 | Study_id), fitted by restricted maximum likelihood.
# Inference: Satterthwaite t-tests and 95% confidence intervals for the fixed X slope.
# Sampling variances and sampling variance-covariance matrices are not used.
# Analyses: both X-Y directions, both snow treatments, Overall and three seasons.
# Figures: Arial; 2050 x 1890 px single plots; 4000 x 1840 px combined plots;
# orange/blue treatment colors; point size cex = 3; annotation cex = 1.00.
# =============================================================================

options(stringsAsFactors = FALSE, scipen = 100)

required_packages <- c("readxl", "dplyr", "tibble", "writexl", "lme4", "lmerTest")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Unable to load package: ", pkg)
  }
}

# Input and output paths. Rename the source folders/files to these English names,
# or edit the paths below to match the local English directory structure.
INPUT_DIR <- "E:/PhD_Research/Meta_Analysis/Final_Dataset_2026-10-29/Effect_Size_Results_2026-10-31"
FILE_QC <- file.path(INPUT_DIR, "Meta_Effect_Sizes_2026-10-29_Nakagawa_Hamman_QC.xlsx")
FILE_D2 <- file.path(INPUT_DIR, "Meta_Effect_Sizes_2026-10-29_d2_Nakagawa_main.rds")
OUTPUT_ROOT <- "E:/PhD_Research/Meta_Analysis/Final_Dataset_2026-10-29/Regression_Results_2026-11-04/04_LMM_Linear_Mixed_Effects_Model"

MIN_N <- 4L
MIN_STUDIES <- 2L
RUN_BOTH_DIRECTIONS <- TRUE
GENERATE_FIGURES <- TRUE
CHECKPOINT_EVERY <- 25L
PERIOD_LEVELS <- c("Overall", "Freezing", "Freeze_thaw", "Growing")
TREATMENTS <- c("Increased snowpack thickness", "Decreased snowpack thickness")
DPI <- 600L
SINGLE_WIDTH <- 2050L
SINGLE_HEIGHT <- 1890L
COMBINED_WIDTH <- 4000L
COMBINED_HEIGHT <- 1840L
ADDITION_POINT_COLOR <- "#DF7E21"
ADDITION_CI_COLOR <- "#FEF1E7"
REMOVAL_POINT_COLOR <- "#5A769E"
REMOVAL_CI_COLOR <- "#E9F3F9"
POINT_CEX <- 3.0
STATS_CEX <- 1.00

TABLE_ROOT <- file.path(OUTPUT_ROOT, "01_Results_Tables")
QC_ROOT <- file.path(OUTPUT_ROOT, "02_QC")
FIGURE_ROOT <- file.path(OUTPUT_ROOT, "03_Regression_Figures")
for (p in c(TABLE_ROOT, QC_ROOT, FIGURE_ROOT)) {
  dir.create(p, recursive = TRUE, showWarnings = FALSE)
}

period_dirs <- function(period) {
  stopifnot(period %in% PERIOD_LEVELS)
  paths <- list(
    TABLE = file.path(TABLE_ROOT, period),
    QC = file.path(QC_ROOT, period),
    SINGLE = file.path(FIGURE_ROOT, period, "Single_Plots"),
    COMBINED = file.path(FIGURE_ROOT, period, "Combined_Plots")
  )
  for (p in unlist(paths, use.names = FALSE)) {
    dir.create(p, recursive = TRUE, showWarnings = FALSE)
  }
  paths
}

LOG_FILE <- file.path(QC_ROOT, "LMM_Run_Log.txt")
writeLines(paste("Run started:", Sys.time()), LOG_FILE)
log_message <- function(...) {
  msg <- paste0(..., collapse = "")
  cat(msg, "\n")
  cat(msg, "\n", file = LOG_FILE, append = TRUE)
}

# -----------------------------------------------------------------------------
# Data import and preprocessing
# -----------------------------------------------------------------------------
if (!all(file.exists(c(FILE_QC, FILE_D2)))) {
  stop("One or more input files do not exist. Check INPUT_DIR, FILE_QC, and FILE_D2.")
}

qc <- as.data.frame(readxl::read_excel(FILE_QC, sheet = 1))
d2 <- as.data.frame(readRDS(FILE_D2))

if (!"Row_id" %in% names(qc) && "row_id" %in% names(qc)) {
  names(qc)[names(qc) == "row_id"] <- "Row_id"
}
if (!"Row_id" %in% names(d2) && "row_id" %in% names(d2)) {
  names(d2)[names(d2) == "row_id"] <- "Row_id"
}

required_d2 <- c("Row_id", "Trait", "Study_id", "Treatment_1", "yi")
missing_d2 <- setdiff(required_d2, names(d2))
if (length(missing_d2)) {
  stop("Required columns are missing from the effect-size data: ", paste(missing_d2, collapse = ", "))
}
if (!"Row_id" %in% names(qc)) stop("The QC workbook is missing Row_id.")
if (!"Pair_ID" %in% names(qc)) stop("The QC workbook is missing Pair_ID.")
if (anyDuplicated(qc$Row_id) || anyDuplicated(d2$Row_id)) {
  stop("Row_id values must be unique in both input tables.")
}

if ("Pair_ID" %in% names(d2)) d2$Pair_ID <- NULL
row_index <- match(d2$Row_id, qc$Row_id)
if (anyNA(row_index)) stop("Some effect-size Row_id values were not found in the QC workbook.")
d2$Pair_ID <- qc$Pair_ID[row_index]

if (!"Sampling_season" %in% names(d2)) {
  if (!"Sampling_season" %in% names(qc)) {
    stop("Sampling_season is missing from both the effect-size data and the QC workbook.")
  }
  d2$Sampling_season <- qc$Sampling_season[row_index]
}

d2$yi <- suppressWarnings(as.numeric(d2$yi))
if (any(!is.finite(d2$yi))) stop("The yi column contains missing or non-finite values.")

clean_id <- function(x) {
  z <- trimws(as.character(x))
  z[is.na(x) | !nzchar(z)] <- NA_character_
  z
}

d2$Row_id <- clean_id(d2$Row_id)
d2$Pair_ID <- clean_id(d2$Pair_ID)
d2$Study_id <- clean_id(d2$Study_id)
d2$Trait <- trimws(as.character(d2$Trait))
if (anyNA(d2$Study_id) || anyNA(d2$Trait) || any(d2$Trait == "")) {
  stop("Study_id or Trait contains missing values.")
}

normalize_treatment <- function(x) {
  z <- tolower(trimws(as.character(x)))
  out <- rep(NA_character_, length(z))
  out[z %in% c("snow addition", "increased snowpack thickness", "increased snow", "snow increase")] <-
    "Increased snowpack thickness"
  out[z %in% c("snow remove", "snow removal", "decreased snowpack thickness", "decreased snow", "snow decrease")] <-
    "Decreased snowpack thickness"
  out
}

normalize_season <- function(x) {
  z <- tolower(trimws(as.character(x)))
  z <- gsub("_", "-", z, fixed = TRUE)
  z <- gsub(" +", " ", z)
  out <- rep(NA_character_, length(z))
  out[z %in% c("annual", "annual average", "overall", "all season", "all seasons",
               "year round", "whole year", "yearly")] <- "Annual_average"
  out[z %in% c("freeze", "freezing", "freeze season", "freezing season", "frozen season")] <- "Freezing"
  out[z %in% c("freeze-thaw", "freeze-thaw season", "freeze thaw", "freeze thaw season",
               "freezethaw", "freezethaw season")] <- "Freeze_thaw"
  out[z %in% c("growing", "growing season")] <- "Growing"
  out
}

d2$Treatment_std <- normalize_treatment(d2$Treatment_1)
d2$Season_std <- normalize_season(d2$Sampling_season)
if (anyNA(d2$Treatment_std)) {
  bad <- unique(d2$Treatment_1[is.na(d2$Treatment_std)])
  stop("Unrecognized Treatment_1 values: ", paste(bad, collapse = " | "))
}
if (anyNA(d2$Season_std)) {
  log_message(
    "Unrecognized Sampling_season values will be included only in Overall: ",
    paste(unique(d2$Sampling_season[is.na(d2$Season_std)]), collapse = " | ")
  )
}

log_message(
  "Imported effect-size records: ", nrow(d2),
  "; traits: ", length(unique(d2$Trait)),
  "; studies: ", length(unique(d2$Study_id))
)

# -----------------------------------------------------------------------------
# Pairing by Pair_ID
# -----------------------------------------------------------------------------
make_pairs <- function(data_source, y_trait, x_trait) {
  y_data <- data_source[
    data_source$Trait == y_trait & !is.na(data_source$Pair_ID),
    c("Study_id", "Pair_ID", "Treatment_std", "yi", "Row_id"),
    drop = FALSE
  ]
  x_data <- data_source[
    data_source$Trait == x_trait & !is.na(data_source$Pair_ID),
    c("Study_id", "Pair_ID", "Treatment_std", "yi", "Row_id"),
    drop = FALSE
  ]
  if (!nrow(y_data) || !nrow(x_data)) return(NULL)

  y_data <- dplyr::add_count(y_data, Study_id, Pair_ID, Treatment_std, name = "nY")
  x_data <- dplyr::add_count(x_data, Study_id, Pair_ID, Treatment_std, name = "nX")

  paired <- dplyr::inner_join(
    y_data, x_data,
    by = c("Study_id", "Pair_ID", "Treatment_std"),
    suffix = c("_Y", "_X")
  )
  paired <- paired[paired$nY == 1L & paired$nX == 1L, , drop = FALSE]
  if (!nrow(paired)) return(NULL)

  tibble::tibble(
    Study_id = as.character(paired$Study_id),
    Pair_ID = as.character(paired$Pair_ID),
    Treatment = as.character(paired$Treatment_std),
    Y_trait = y_trait,
    X_trait = x_trait,
    Y = as.numeric(paired$yi_Y),
    X = as.numeric(paired$yi_X),
    Row_id_Y = as.character(paired$Row_id_Y),
    Row_id_X = as.character(paired$Row_id_X)
  )
}

pairing_diagnostic <- function(data_source, y_trait, x_trait) {
  y_data <- data_source[
    data_source$Trait == y_trait & !is.na(data_source$Pair_ID),
    c("Study_id", "Pair_ID", "Treatment_std"),
    drop = FALSE
  ]
  x_data <- data_source[
    data_source$Trait == x_trait & !is.na(data_source$Pair_ID),
    c("Study_id", "Pair_ID", "Treatment_std"),
    drop = FALSE
  ]

  y_keys <- if (nrow(y_data)) {
    dplyr::count(y_data, Study_id, Pair_ID, Treatment_std, name = "nY")
  } else {
    data.frame(Study_id = character(), Pair_ID = character(), Treatment_std = character(), nY = integer())
  }
  x_keys <- if (nrow(x_data)) {
    dplyr::count(x_data, Study_id, Pair_ID, Treatment_std, name = "nX")
  } else {
    data.frame(Study_id = character(), Pair_ID = character(), Treatment_std = character(), nX = integer())
  }

  common_keys <- dplyr::inner_join(
    y_keys, x_keys,
    by = c("Study_id", "Pair_ID", "Treatment_std")
  )
  valid_keys <- if (nrow(common_keys)) {
    sum(common_keys$nY == 1L & common_keys$nX == 1L)
  } else 0L
  ambiguous_keys <- if (nrow(common_keys)) {
    sum(!(common_keys$nY == 1L & common_keys$nX == 1L))
  } else 0L

  paired <- tryCatch(make_pairs(data_source, y_trait, x_trait), error = function(e) NULL)
  if (is.null(paired)) paired <- data.frame()

  data.frame(
    Y_trait = y_trait,
    X_trait = x_trait,
    Y_records_with_Pair_ID = nrow(y_data),
    X_records_with_Pair_ID = nrow(x_data),
    Common_Pair_ID_keys = nrow(common_keys),
    Valid_unique_Pair_ID_keys = valid_keys,
    Ambiguous_Pair_ID_keys = ambiguous_keys,
    Final_paired_observations = nrow(paired),
    Increased_snow_n = if (nrow(paired)) sum(paired$Treatment == "Increased snowpack thickness") else 0L,
    Decreased_snow_n = if (nrow(paired)) sum(paired$Treatment == "Decreased snowpack thickness") else 0L,
    stringsAsFactors = FALSE
  )
}

# -----------------------------------------------------------------------------
# Linear mixed-effects model
# -----------------------------------------------------------------------------
fit_lmm <- function(data) {
  if (nrow(data) < MIN_N) stop("Insufficient paired observations for LMM.")
  if (any(!is.finite(data$X)) || any(!is.finite(data$Y))) {
    stop("Non-finite X or Y values were found in the model data.")
  }
  if (length(unique(data$X)) < 2L || !is.finite(stats::var(data$X)) || stats::var(data$X) < 1e-12) {
    stop("The predictor has insufficient variation.")
  }

  model_data <- as.data.frame(data)
  model_data$Study_id <- factor(model_data$Study_id)
  n_studies <- nlevels(model_data$Study_id)
  if (n_studies < MIN_STUDIES) {
    stop("At least two independent Study_id levels are required; found ", n_studies, ".")
  }
  if (n_studies >= nrow(model_data)) {
    stop("The number of Study_id levels must be smaller than the number of paired observations.")
  }

  control <- lme4::lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 100000)
  )
  model <- lmerTest::lmer(
    Y ~ X + (1 | Study_id),
    data = model_data,
    REML = TRUE,
    control = control
  )

  model_summary <- summary(model)
  coefficient_table <- as.matrix(model_summary$coefficients)
  if (!"X" %in% rownames(coefficient_table)) {
    stop("The fitted LMM did not return the fixed-effect slope for X.")
  }
  required_statistics <- c("Estimate", "Std. Error", "df", "Pr(>|t|)")
  if (!all(required_statistics %in% colnames(coefficient_table))) {
    stop("The expected Satterthwaite coefficient statistics are unavailable.")
  }

  beta <- as.numeric(coefficient_table["X", "Estimate"])
  standard_error <- as.numeric(coefficient_table["X", "Std. Error"])
  degrees_freedom <- as.numeric(coefficient_table["X", "df"])
  p_value <- as.numeric(coefficient_table["X", "Pr(>|t|)"])
  if (!all(is.finite(c(beta, standard_error, degrees_freedom, p_value))) ||
      standard_error <= 0 || degrees_freedom <= 0) {
    stop("The LMM returned an invalid slope estimate or inferential statistic.")
  }

  critical_value <- stats::qt(0.975, df = degrees_freedom)
  slope_ci <- beta + c(-1, 1) * critical_value * standard_error

  variance_table <- as.data.frame(lme4::VarCorr(model))
  study_variance <- variance_table$vcov[
    variance_table$grp == "Study_id" &
      variance_table$var1 == "(Intercept)" &
      is.na(variance_table$var2)
  ]
  study_variance <- if (length(study_variance)) as.numeric(study_variance[1]) else NA_real_
  residual_variance <- stats::sigma(model)^2

  fixed_predictions <- as.numeric(stats::predict(model, re.form = NA))
  fixed_effect_variance <- if (length(fixed_predictions) > 1L) stats::var(fixed_predictions) else 0
  total_variance <- fixed_effect_variance + study_variance + residual_variance
  marginal_r2_approx <- if (is.finite(total_variance) && total_variance > 0) {
    fixed_effect_variance / total_variance
  } else {
    NA_real_
  }
  conditional_r2_approx <- if (is.finite(total_variance) && total_variance > 0) {
    (fixed_effect_variance + study_variance) / total_variance
  } else {
    NA_real_
  }

  convergence_messages <- model@optinfo$conv$lme4$messages
  if (is.null(convergence_messages)) convergence_messages <- character(0)

  list(
    Model = "LMM",
    Beta = beta,
    SE = standard_error,
    CI_low = slope_ci[1],
    CI_high = slope_ci[2],
    p = p_value,
    DF_Satterthwaite = degrees_freedom,
    Intercept = as.numeric(lme4::fixef(model)["(Intercept)"]),
    Study_variance = study_variance,
    Residual_variance = residual_variance,
    R2_marginal_approx = marginal_r2_approx,
    R2_conditional_approx = conditional_r2_approx,
    Is_singular = lme4::isSingular(model, tol = 1e-4),
    Convergence_message = paste(convergence_messages, collapse = " | "),
    n = nrow(model_data),
    Study_n = n_studies,
    Fit = model,
    CI_method = "Satterthwaite t-based 95% confidence interval",
    P_method = "Satterthwaite two-sided t-test"
  )
}

# -----------------------------------------------------------------------------
# Figure utilities
# -----------------------------------------------------------------------------
if (.Platform$OS.type == "windows") {
  try(grDevices::windowsFonts(Arial = grDevices::windowsFont("Arial")), silent = TRUE)
}

safe_name <- function(x) {
  z <- as.character(x)
  z <- gsub("βG", "BetaG", z, fixed = TRUE)
  z <- gsub("β", "Beta", z, fixed = TRUE)
  z <- gsub("[\\\\/:*?\"<>|]", "_", z)
  z <- gsub("[[:space:]]+", "_", z)
  z <- gsub("[. ]+$", "", z)
  substr(ifelse(nzchar(z), z, "unknown"), 1, 85)
}

format_p <- function(p) {
  if (!is.finite(p)) return("NA")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

compact_axis <- function(values, pad = c(0.08, 0.08)) {
  values <- values[is.finite(values)]
  value_range <- if (length(values)) range(values) else c(-1, 1)
  if (diff(value_range) <= 0) {
    span <- max(abs(value_range[1]) * 0.2, 0.1)
    value_range <- value_range + c(-0.5, 0.5) * span
  }
  target <- value_range + c(-pad[1], pad[2]) * diff(value_range)
  raw_step <- diff(target) / 3
  quantum <- 10^(floor(log10(raw_step)) - 1)
  if (raw_step >= 0.01) quantum <- max(0.01, quantum)
  step <- ceiling(raw_step / quantum) * quantum
  lower <- floor(target[1] / quantum) * quantum
  if (lower + 3 * step < target[2]) step <- step + quantum
  ticks <- lower + (0:3) * step
  list(ticks = ticks, range = ticks[c(1, 4)])
}

format_axis_labels <- function(x) {
  labels <- format(round(x, 2), trim = TRUE, scientific = FALSE)
  if (anyDuplicated(labels)) labels <- formatC(x, format = "e", digits = 1)
  labels
}

x_limits_with_label_space <- function(ticks, cex = 1) {
  label_widths <- graphics::strwidth(
    format_axis_labels(ticks[c(1, 4)]), units = "inches", cex = cex, family = "Arial"
  )
  fraction <- max(0.065, (max(label_widths) / 2 + 0.035) / graphics::par("pin")[1])
  fraction <- min(fraction, 0.22)
  extra <- diff(range(ticks)) * fraction / (1 - 2 * fraction)
  range(ticks) + c(-extra, extra)
}

draw_x_axis <- function(ticks, cex = 1) {
  labels <- format_axis_labels(ticks)
  widths <- graphics::strwidth(labels, cex = cex, family = "Arial")
  required_spacing <- (widths[-4] + widths[-1]) / 2
  available_edges <- c(ticks[1] - graphics::par("usr")[1], graphics::par("usr")[2] - ticks[4])
  edge_fit <- min(available_edges / pmax(widths[c(1, 4)] * 0.55, .Machine$double.eps))
  label_cex <- cex * min(
    1,
    edge_fit,
    min(diff(ticks) / pmax(required_spacing * 1.15, .Machine$double.eps))
  )
  graphics::axis(1, at = ticks, labels = FALSE, tck = -0.02, lwd = 1.2)
  for (i in seq_along(ticks)) {
    graphics::mtext(
      labels[i], side = 1, at = ticks[i], line = 0.65, adj = 0.5,
      cex = label_cex, family = "Arial", font = 1, xpd = NA
    )
  }
}

draw_y_axis <- function(ticks, cex = 1, side = 2, labels = format_axis_labels(ticks)) {
  graphics::axis(side, at = ticks, labels = FALSE, tck = -0.02, lwd = 1.2)
  for (i in seq_along(ticks)) {
    graphics::mtext(
      labels[i], side = side, at = ticks[i], line = 0.6,
      las = 1, cex = cex, font = 1, family = "Arial", xpd = NA
    )
  }
}

point_size <- function(n) rep(POINT_CEX, n)

draw_statistics <- function(x, y, fit, n, cex = STATS_CEX) {
  if (is.null(fit)) return(invisible(NULL))
  p_label <- format_p(fit$p)
  items <- c(
    "β", " = ", sprintf("%.3f", fit$Beta), ", ", "p",
    if (p_label == "<0.001") " " else " = ", p_label,
    ", ", "n", " = ", as.character(n)
  )
  fonts <- c(3, 1, 1, 1, 3, 1, 1, 1, 3, 1, 1)
  widths <- mapply(
    function(label, font) graphics::strwidth(label, units = "user", cex = cex,
                                               family = "Arial", font = font),
    items, fonts
  )
  x_start <- x - sum(widths) / 2
  for (i in seq_along(items)) {
    graphics::text(
      x_start, y, items[i], adj = c(0, 1), family = "Arial",
      font = fonts[i], cex = cex, xpd = FALSE
    )
    x_start <- x_start + widths[i]
  }
  invisible(NULL)
}

make_predictions <- function(data, fit) {
  x_grid <- seq(min(data$X), max(data$X), length.out = 150)
  predicted <- fit$Intercept + fit$Beta * x_grid
  fixed_covariance <- as.matrix(stats::vcov(fit$Fit))[
    c("(Intercept)", "X"), c("(Intercept)", "X"), drop = FALSE
  ]
  design_matrix <- cbind("(Intercept)" = rep(1, length(x_grid)), X = x_grid)
  standard_error_mean <- sqrt(pmax(0, rowSums((design_matrix %*% fixed_covariance) * design_matrix)))
  critical_value <- stats::qt(0.975, df = fit$DF_Satterthwaite)
  list(
    x = x_grid,
    y = as.numeric(predicted),
    lower = as.numeric(predicted - critical_value * standard_error_mean),
    upper = as.numeric(predicted + critical_value * standard_error_mean)
  )
}

plot_panel <- function(data, fit, treatment, shared_y_axis = NULL, show_left_axis = TRUE) {
  if (is.null(fit) || is.null(data) || !nrow(data)) {
    graphics::plot.new()
    graphics::text(0.5, 0.5, "Insufficient data", family = "Arial")
    return(invisible(NULL))
  }

  predictions <- make_predictions(data, fit)
  x_axis <- compact_axis(data$X, pad = c(0.12, 0.12))
  x_limits <- x_limits_with_label_space(x_axis$ticks, cex = 1)
  y_values <- c(data$Y, predictions$y, predictions$lower, predictions$upper)
  y_values <- y_values[is.finite(y_values)]
  y_axis <- if (is.null(shared_y_axis)) compact_axis(y_values, pad = c(0.08, 0.20)) else shared_y_axis

  is_addition <- treatment == "Increased snowpack thickness"
  point_color <- if (is_addition) ADDITION_POINT_COLOR else REMOVAL_POINT_COLOR
  confidence_color <- if (is_addition) ADDITION_CI_COLOR else REMOVAL_CI_COLOR

  graphics::plot(
    NA, xlim = x_limits, ylim = y_axis$range,
    xlab = "", ylab = "", axes = FALSE, xaxs = "i", yaxs = "i"
  )
  if (all(is.finite(predictions$lower)) && all(is.finite(predictions$upper))) {
    graphics::polygon(
      c(predictions$x, rev(predictions$x)),
      c(predictions$lower, rev(predictions$upper)),
      col = confidence_color, border = NA
    )
    graphics::lines(predictions$x, predictions$lower, col = "gray60", lty = 2, lwd = 0.8)
    graphics::lines(predictions$x, predictions$upper, col = "gray60", lty = 2, lwd = 0.8)
  }
  graphics::lines(predictions$x, predictions$y, col = "gray10", lwd = 1.4)
  graphics::abline(h = 0, col = "gray40", lty = 3, lwd = 0.8)
  graphics::points(
    data$X, data$Y, pch = 21, bg = point_color, col = "gray20", lwd = 0.5,
    cex = point_size(nrow(data)), xpd = FALSE
  )
  draw_x_axis(x_axis$ticks, cex = 1)
  if (show_left_axis) draw_y_axis(y_axis$ticks, cex = 1, side = 2)
  graphics::box(lwd = 1.2)
  draw_statistics(
    mean(x_axis$range), y_axis$range[2] - 0.05 * diff(y_axis$range),
    fit, fit$n, cex = STATS_CEX
  )
  invisible(NULL)
}

open_png <- function(file, width, height) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  png_type <- if (isTRUE(capabilities("cairo"))) "cairo" else if (.Platform$OS.type == "windows") "windows" else "Xlib"
  grDevices::png(
    file, width = width, height = height, units = "px",
    res = DPI, bg = "white", type = png_type
  )
}

plot_single <- function(data, fit, treatment, path) {
  if (is.null(fit) || is.null(data) || !nrow(data)) return(invisible(FALSE))
  open_png(path, SINGLE_WIDTH, SINGLE_HEIGHT)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(
    family = "Arial", oma = c(1.5, 1.5, 0, 0), mar = c(1.5, 2.5, 1, 0.5),
    mgp = c(1.5, 0.2, 0), tcl = -0.02, cex.axis = 1,
    las = 1, lwd = 1.2, font = 1, pty = "s"
  )
  plot_panel(data, fit, treatment)
  invisible(TRUE)
}

plot_beta_panel <- function(addition_fit, removal_fit) {
  fits <- list(addition_fit, removal_fit)
  beta_values <- vapply(fits, function(x) if (is.null(x)) NA_real_ else x$Beta, numeric(1))
  valid <- which(is.finite(beta_values))
  if (!length(valid)) {
    graphics::plot.new()
    graphics::text(0.5, 0.5, "No slope estimates", family = "Arial")
    return(invisible(NULL))
  }

  lower <- min(c(0, beta_values[valid]))
  upper <- max(c(0, beta_values[valid]))
  span <- upper - lower
  if (!is.finite(span) || span < 1e-8) {
    lower <- -0.1
    upper <- 0.1
    span <- 0.2
  }
  padding <- 0.08 * span
  if (lower == 0) {
    y_range <- c(0, upper + padding)
  } else if (upper == 0) {
    y_range <- c(lower - padding, 0)
  } else if (lower < 0 && upper > 0) {
    limit <- max(abs(c(lower, upper))) * 1.08
    y_range <- c(-limit, limit)
  } else {
    y_range <- c(lower - padding, upper + padding)
  }

  ticks <- seq(y_range[1], y_range[2], length.out = 3)
  graphics::plot(
    NA, xlim = c(0.5, 2.5), ylim = y_range,
    axes = FALSE, xlab = "", ylab = "", xaxs = "i", yaxs = "i"
  )
  graphics::abline(h = 0, lty = 3, col = "gray40", lwd = 0.8)
  bar_colors <- c(ADDITION_POINT_COLOR, REMOVAL_POINT_COLOR)
  for (i in valid) {
    graphics::rect(
      i - 0.3, min(0, beta_values[i]), i + 0.3, max(0, beta_values[i]),
      col = bar_colors[i], border = gray(0.2), lwd = 0.5
    )
  }
  draw_y_axis(ticks, side = 4, labels = formatC(ticks, format = "f", digits = 1), cex = 1)
  graphics::box(lwd = 1.2)
  invisible(NULL)
}

plot_combined <- function(addition_data, addition_fit, removal_data, removal_fit, path) {
  if (is.null(addition_fit) && is.null(removal_fit)) return(invisible(FALSE))

  open_png(path, COMBINED_WIDTH, COMBINED_HEIGHT)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(family = "Arial", oma = c(0, 0, 0, 0), font = 1,
                font.axis = 1, font.lab = 1, pty = "m")

  y_values <- numeric(0)
  for (item in list(
    list(data = addition_data, fit = addition_fit),
    list(data = removal_data, fit = removal_fit)
  )) {
    if (!is.null(item$fit) && !is.null(item$data) && nrow(item$data)) {
      prediction <- make_predictions(item$data, item$fit)
      y_values <- c(y_values, item$data$Y, prediction$y, prediction$lower, prediction$upper)
    }
  }
  y_values <- y_values[is.finite(y_values)]
  if (!length(y_values)) stop("No valid response values are available for the combined plot.")
  shared_y_axis <- compact_axis(y_values, pad = c(0.10, 0.22))

  panel_edges <- c(0.09, 0.43, 0.77, 0.90)
  for (panel_index in 1:2) {
    graphics::par(
      fig = c(panel_edges[panel_index], panel_edges[panel_index + 1], 0, 1),
      new = panel_index > 1,
      mai = c(0.62, 0, 0.18, 0), mgp = c(1.5, 0.6, 0),
      tcl = -0.02, las = 1, lwd = 1.2
    )
    if (panel_index == 1) {
      plot_panel(
        addition_data, addition_fit, TREATMENTS[1],
        shared_y_axis = shared_y_axis, show_left_axis = TRUE
      )
    } else {
      plot_panel(
        removal_data, removal_fit, TREATMENTS[2],
        shared_y_axis = shared_y_axis, show_left_axis = FALSE
      )
    }
    if (abs(diff(graphics::par("pin"))) > 1 / 600) {
      stop("A regression panel does not match the specified square-panel geometry.")
    }
  }

  graphics::par(
    fig = c(panel_edges[3], panel_edges[4], 0, 1), new = TRUE,
    mai = c(0.62, 0, 0.18, 0), mgp = c(1.5, 0.6, 0),
    tcl = -0.02, las = 1, lwd = 1.2
  )
  plot_beta_panel(addition_fit, removal_fit)
  invisible(TRUE)
}

# -----------------------------------------------------------------------------
# Analysis workflow
# -----------------------------------------------------------------------------
get_period_data <- function(period) {
  if (period == "Overall") return(d2)
  d2[!is.na(d2$Season_std) & d2$Season_std == period, , drop = FALSE]
}

results <- list()
failures <- list()
pairing_audits <- list()
pairing_diagnostics <- list()
result_index <- 0L
failure_index <- 0L
audit_index <- 0L
diagnostic_index <- 0L
n_single_figures <- 0L
n_combined_figures <- 0L

append_failure <- function(period, y_trait, x_trait, treatment, stage, reason, n = NA_integer_) {
  failure_index <<- failure_index + 1L
  failures[[failure_index]] <<- data.frame(
    Period = period, Y_trait = y_trait, X_trait = x_trait,
    Treatment = treatment, Stage = stage, n = n,
    Reason = as.character(reason), stringsAsFactors = FALSE
  )
}

append_result <- function(period, y_trait, x_trait, treatment, data, fit) {
  result_index <<- result_index + 1L
  results[[result_index]] <<- data.frame(
    Period = period,
    Y_trait = y_trait,
    X_trait = x_trait,
    Treatment = treatment,
    Model = fit$Model,
    Formula = "Y ~ X + (1 | Study_id)",
    Beta = fit$Beta,
    SE = fit$SE,
    CI_low = fit$CI_low,
    CI_high = fit$CI_high,
    p = fit$p,
    DF_Satterthwaite = fit$DF_Satterthwaite,
    n = nrow(data),
    Study_n = fit$Study_n,
    Pair_ID_n = nrow(data),
    Intercept = fit$Intercept,
    Study_variance = fit$Study_variance,
    Residual_variance = fit$Residual_variance,
    R2_marginal_approx = fit$R2_marginal_approx,
    R2_conditional_approx = fit$R2_conditional_approx,
    Is_singular = fit$Is_singular,
    Convergence_message = fit$Convergence_message,
    CI_method = fit$CI_method,
    P_method = fit$P_method,
    stringsAsFactors = FALSE
  )
}

write_csv <- function(x, path) {
  utils::write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
}

write_checkpoint <- function() {
  if (length(results)) {
    write_csv(dplyr::bind_rows(results), file.path(TABLE_ROOT, "LMM_Results_Checkpoint.csv"))
  }
  if (length(failures)) {
    write_csv(dplyr::bind_rows(failures), file.path(QC_ROOT, "LMM_Failures_Checkpoint.csv"))
  }
}

for (period in PERIOD_LEVELS) {
  output_dirs <- period_dirs(period)
  period_data <- get_period_data(period)
  log_message("Processing period: ", period, "; source effect-size records: ", nrow(period_data))

  if (!nrow(period_data)) {
    append_failure(period, "ALL", "ALL", "All", "Period", "No effect-size records in this period", 0L)
    next
  }

  traits <- sort(unique(period_data$Trait))
  if (length(traits) < 2L) {
    append_failure(period, "ALL", "ALL", "All", "Period", "Fewer than two traits are available", length(traits))
    next
  }

  trait_pairs <- utils::combn(traits, 2, simplify = FALSE)
  for (pair_index in seq_along(trait_pairs)) {
    trait_pair <- trait_pairs[[pair_index]]
    directions <- if (RUN_BOTH_DIRECTIONS) {
      list(c(trait_pair[1], trait_pair[2]), c(trait_pair[2], trait_pair[1]))
    } else {
      list(trait_pair)
    }

    for (direction in directions) {
      y_trait <- direction[1]
      x_trait <- direction[2]

      diagnostic <- tryCatch(
        pairing_diagnostic(period_data, y_trait, x_trait),
        error = function(e) data.frame(
          Y_trait = y_trait, X_trait = x_trait,
          Diagnostic_error = conditionMessage(e), stringsAsFactors = FALSE
        )
      )
      diagnostic$Period <- period
      diagnostic_index <- diagnostic_index + 1L
      pairing_diagnostics[[diagnostic_index]] <- diagnostic

      paired_data <- tryCatch(
        make_pairs(period_data, y_trait, x_trait),
        error = function(e) {
          append_failure(period, y_trait, x_trait, "All", "Pairing", conditionMessage(e))
          NULL
        }
      )

      if (is.null(paired_data) || !nrow(paired_data)) {
        append_failure(
          period, y_trait, x_trait, "All", "Pairing",
          "No unique matched records were found for the selected Pair_ID values.", 0L
        )
        next
      }

      treatment_data <- list()
      treatment_fits <- list()
      for (treatment in TREATMENTS) {
        model_data <- paired_data[paired_data$Treatment == treatment, , drop = FALSE]
        model_data$Period <- period
        treatment_data[[treatment]] <- model_data
        n_observations <- nrow(model_data)
        n_studies <- length(unique(model_data$Study_id))

        audit_index <- audit_index + 1L
        pairing_audits[[audit_index]] <- data.frame(
          Period = period,
          Y_trait = y_trait,
          X_trait = x_trait,
          Treatment = treatment,
          n = n_observations,
          Study_n = n_studies,
          Pair_ID_n = n_observations,
          stringsAsFactors = FALSE
        )

        if (n_observations < MIN_N) {
          append_failure(
            period, y_trait, x_trait, treatment, "LMM_fit",
            paste0("Fewer than ", MIN_N, " paired observations."), n_observations
          )
          next
        }
        if (!is.finite(stats::var(model_data$X)) || stats::var(model_data$X) < 1e-12) {
          append_failure(
            period, y_trait, x_trait, treatment, "LMM_fit",
            "Insufficient variation in the predictor X.", n_observations
          )
          next
        }

        fit <- tryCatch(
          fit_lmm(model_data),
          error = function(e) {
            append_failure(period, y_trait, x_trait, treatment, "LMM_fit",
                           conditionMessage(e), n_observations)
            NULL
          }
        )
        if (is.null(fit)) next

        append_result(period, y_trait, x_trait, treatment, model_data, fit)
        treatment_fits[[treatment]] <- fit

        if (fit$Is_singular || nzchar(fit$Convergence_message)) {
          append_failure(
            period, y_trait, x_trait, treatment, "LMM_diagnostic",
            paste0("Model diagnostic: singular = ", fit$Is_singular,
                   "; convergence message = ", fit$Convergence_message),
            n_observations
          )
        }

        if (GENERATE_FIGURES) {
          treatment_suffix <- if (treatment == TREATMENTS[1]) "Increased" else "Decreased"
          figure_name <- paste0(safe_name(y_trait), "__", safe_name(x_trait),
                                "__", treatment_suffix, ".png")
          tryCatch({
            created <- plot_single(
              model_data, fit, treatment,
              file.path(output_dirs$SINGLE, figure_name)
            )
            if (isTRUE(created)) n_single_figures <- n_single_figures + 1L
          }, error = function(e) {
            append_failure(period, y_trait, x_trait, treatment, "Single_plot",
                           conditionMessage(e), n_observations)
          })
        }
      }

      if (GENERATE_FIGURES && any(!vapply(treatment_fits, is.null, logical(1)))) {
        combined_name <- paste0(safe_name(y_trait), "__", safe_name(x_trait), "__Combined.png")
        tryCatch({
          created <- plot_combined(
            treatment_data[[TREATMENTS[1]]], treatment_fits[[TREATMENTS[1]]],
            treatment_data[[TREATMENTS[2]]], treatment_fits[[TREATMENTS[2]]],
            file.path(output_dirs$COMBINED, combined_name)
          )
          if (isTRUE(created)) n_combined_figures <- n_combined_figures + 1L
        }, error = function(e) {
          append_failure(period, y_trait, x_trait, "Both", "Combined_plot", conditionMessage(e))
        })
      }
    }

    if (pair_index %% CHECKPOINT_EVERY == 0L || pair_index == length(trait_pairs)) {
      log_message(
        "Period ", period, ": completed trait pairs ", pair_index, "/", length(trait_pairs),
        "; successful LMMs = ", length(results),
        "; single figures = ", n_single_figures,
        "; combined figures = ", n_combined_figures
      )
      write_checkpoint()
      gc(verbose = FALSE)
    }
  }
}

# -----------------------------------------------------------------------------
# Result processing and export
# -----------------------------------------------------------------------------
all_results <- if (length(results)) dplyr::bind_rows(results) else data.frame()
all_failures <- if (length(failures)) dplyr::bind_rows(failures) else data.frame()
all_pairing_audits <- if (length(pairing_audits)) dplyr::bind_rows(pairing_audits) else data.frame()
all_pairing_diagnostics <- if (length(pairing_diagnostics)) dplyr::bind_rows(pairing_diagnostics) else data.frame()

if (nrow(all_results)) {
  all_results <- all_results |>
    dplyr::group_by(Period) |>
    dplyr::mutate(q_FDR = stats::p.adjust(p, method = "BH")) |>
    dplyr::ungroup() |>
    dplyr::arrange(Period, Y_trait, X_trait, Treatment)
} else {
  all_results <- data.frame(Note = "No LMM models were successfully fitted.")
}
if (!nrow(all_failures)) all_failures <- data.frame(Note = "No failed or flagged models.")
if (!nrow(all_pairing_audits)) all_pairing_audits <- data.frame(Note = "No eligible paired datasets.")
if (!nrow(all_pairing_diagnostics)) all_pairing_diagnostics <- data.frame(Note = "No pairing diagnostics available.")

safe_write_xlsx <- function(sheets, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary_file <- tempfile(pattern = ".lmm_tmp_", tmpdir = dirname(path), fileext = ".xlsx")
  on.exit(if (file.exists(temporary_file)) unlink(temporary_file), add = TRUE)
  writexl::write_xlsx(sheets, path = temporary_file)

  if (!file.exists(path) && isTRUE(file.rename(temporary_file, path))) {
    return(invisible(path))
  }
  if (file.exists(path) && isTRUE(suppressWarnings(file.remove(path))) &&
      isTRUE(file.rename(temporary_file, path))) {
    return(invisible(path))
  }

  file_stem <- sub("\\.xlsx$", "", basename(path), ignore.case = TRUE)
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  backup <- file.path(dirname(path), paste0(file_stem, "_", timestamp, ".xlsx"))
  suffix <- 1L
  while (file.exists(backup)) {
    backup <- file.path(dirname(path), paste0(file_stem, "_", timestamp, "_", suffix, ".xlsx"))
    suffix <- suffix + 1L
  }
  if (!isTRUE(file.rename(temporary_file, backup)) &&
      !isTRUE(file.copy(temporary_file, backup, overwrite = FALSE))) {
    stop("Unable to save the Excel output: ", backup)
  }
  log_message("The target workbook could not be replaced; a timestamped copy was saved: ", backup)
  invisible(backup)
}

period_table <- function(data, period, empty_note) {
  out <- if ("Period" %in% names(data)) data[data$Period == period, , drop = FALSE] else data[FALSE, , drop = FALSE]
  if (!nrow(out)) out <- data.frame(Note = empty_note)
  as.data.frame(out)
}

write_csv(all_results, file.path(TABLE_ROOT, "All_LMM_Results.csv"))
write_csv(all_failures, file.path(QC_ROOT, "LMM_Failed_or_Flagged.csv"))
write_csv(all_pairing_audits, file.path(QC_ROOT, "Pair_ID_Pairing_QC.csv"))
write_csv(all_pairing_diagnostics, file.path(QC_ROOT, "Pair_ID_Pairing_Diagnostics.csv"))

input_summary <- data.frame(
  Item = c("Effect_size_records", "Number_of_traits", "Number_of_studies",
           "Pairing_identifier", "Model", "Estimation", "Fixed_effect_inference",
           "Sampling_variance_weighting"),
  Value = c(
    as.character(nrow(d2)),
    as.character(length(unique(d2$Trait))),
    as.character(length(unique(d2$Study_id))),
    "Pair_ID",
    "Y ~ X + (1 | Study_id)",
    "Restricted maximum likelihood (REML)",
    "Satterthwaite approximation",
    "Not used"
  ),
  stringsAsFactors = FALSE
)

safe_write_xlsx(
  list(
    LMM_Regression = as.data.frame(all_results),
    Failed_or_Flagged = as.data.frame(all_failures),
    Pair_ID_Pairing_QC = as.data.frame(all_pairing_audits),
    Pair_ID_Pairing_Diagnostics = as.data.frame(all_pairing_diagnostics),
    Input_and_Methods = input_summary
  ),
  file.path(TABLE_ROOT, "All_Traits_LMM_Results_and_QC.xlsx")
)

for (period in PERIOD_LEVELS) {
  output_dirs <- period_dirs(period)
  period_results <- period_table(all_results, period, "No successful LMM models in this period.")
  period_failures <- period_table(all_failures, period, "No failed or flagged models in this period.")
  period_pairing <- period_table(all_pairing_audits, period, "No eligible paired datasets in this period.")
  period_diagnostics <- period_table(all_pairing_diagnostics, period, "No Pair_ID diagnostics in this period.")

  write_csv(period_results, file.path(output_dirs$TABLE, "LMM_Regression_Results.csv"))
  write_csv(period_failures, file.path(output_dirs$QC, "LMM_Failed_or_Flagged.csv"))
  write_csv(period_pairing, file.path(output_dirs$QC, "Pair_ID_Pairing_QC.csv"))
  write_csv(period_diagnostics, file.path(output_dirs$QC, "Pair_ID_Pairing_Diagnostics.csv"))

  safe_write_xlsx(
    list(
      LMM_Regression = period_results,
      Failed_or_Flagged = period_failures,
      Pair_ID_Pairing_QC = period_pairing,
      Pair_ID_Pairing_Diagnostics = period_diagnostics
    ),
    file.path(output_dirs$TABLE, paste0("LMM_", period, "_Results_and_QC.xlsx"))
  )
}

saveRDS(
  list(
    regression = all_results,
    failures = all_failures,
    pairing_qc = all_pairing_audits,
    pairing_diagnostics = all_pairing_diagnostics,
    formula = "Y ~ X + (1 | Study_id)",
    estimation = "lmerTest::lmer, REML = TRUE",
    inference = "Satterthwaite t-test and 95% slope confidence interval",
    pairing_identifier = "Pair_ID",
    sampling_variances_used = FALSE,
    input_files = c(FILE_QC, FILE_D2)
  ),
  file.path(TABLE_ROOT, "All_Traits_LMM_Results.rds")
)

writeLines(
  c(
    "Model: Y ~ X + (1 | Study_id).",
    "Pairing was performed using Pair_ID.",
    "A paired record was retained only when the corresponding Study_id, Pair_ID, and treatment combination contained one record for each indicator.",
    "The model includes a study-level random intercept.",
    "Sampling variances and sampling variance-covariance matrices are not used.",
    "The predictor effect size is treated as an observed covariate without an explicit measurement-error model.",
    "The slope p-value and 95% confidence interval use the Satterthwaite degrees-of-freedom approximation.",
    "Marginal and conditional R-squared values are approximate descriptive quantities.",
    "Review singular-fit and convergence diagnostics before interpreting individual models."
  ),
  file.path(OUTPUT_ROOT, "Methods_and_Data_Notes.txt")
)

capture.output(utils::sessionInfo(), file = file.path(OUTPUT_ROOT, "Session_Info.txt"))

log_message("Completed LMM models: ", if ("Model" %in% names(all_results)) sum(all_results$Model == "LMM") else 0L)
log_message("Single figures generated: ", n_single_figures)
log_message("Combined figures generated: ", n_combined_figures)
log_message("Output directory: ", OUTPUT_ROOT)
log_message("Model formula: Y ~ X + (1 | Study_id); REML = TRUE; Satterthwaite inference.")
log_message("Pairing identifier: Pair_ID; sampling variances and covariance matrices were not used.")
