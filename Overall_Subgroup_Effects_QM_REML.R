# =============================================================================
# Overall and subgroup meta-analysis: treatment effects, heterogeneity, and QM
#
# English companion of the final Chinese v4 analysis, revised to include QE.
# Inputs: Hedges' g (yi), sampling variance (vi), and the complete Nakagawa
# shared-control sampling variance-covariance matrix V from Part 1.
# Effect models: separate REML intercept-only rma.mv() fits for snow addition
# and removal, with random = ~ 1 | Study_id / Row_id.
# For each trait and subgroup, the treatment moderator model fits
# mods = ~ Treatment_group and tests removal minus addition with QM (df = 1).
# The residual heterogeneity test QE, its df, and p-value are reported from
# the treatment-specific and joint treatment models.
# Included: overall, sampling season, frozen soil, duration, and depth.
# Ecosystem classification is excluded. Original StandAge_year duration
# cutoffs from the Chinese v4 analysis are preserved unchanged.
# At least three effect sizes per treatment are required for model fitting.
# QM in this script tests the treatment contrast within each subgroup, not
# omnibus moderation by season, soil type, experimental duration or depth.
# =============================================================================

rm(list = ls())
gc()


# =============================================================================
# 1. Dependencies
# =============================================================================

required_packages <- c(
  "metafor",
  "Matrix",
  "dplyr",
  "writexl"
)

new_packages <- required_packages[
  !(required_packages %in% rownames(installed.packages()))
]

if (length(new_packages) > 0) {
  install.packages(new_packages)
}

library(metafor)
library(Matrix)
library(dplyr)
library(writexl)

MIN_EFFECTS_PER_TREATMENT <- 3L


# =============================================================================
# 2. Input and output paths
# =============================================================================

input_dir <- "results"

data_file <- file.path(
  input_dir,
  "Meta_effect_sizes_2026-10-29_d2_Nakagawa_main.rds"
)

V_file <- file.path(
  input_dir,
  "Meta_effect_sizes_2026-10-29_V_Nakagawa_main.rds"
)

output_dir <- file.path(
  input_dir,
  "overall_subgroup_QM_QE_no_ecosystem"
)

if (!dir.exists(output_dir)) {
  dir.create(
    output_dir,
    recursive = TRUE
  )
}

output_excel <- file.path(
  output_dir,
  "meta_overall_subgroup_QM_QE.xlsx"
)

output_rds <- file.path(
  output_dir,
  "meta_overall_subgroup_QM_QE.rds"
)


cat("\n")
cat(strrep("=", 90), "\n")
cat("Overall and subgroup treatment effects, QM, and QE (original StandAge_year bins)\n")
cat(strrep("=", 90), "\n\n")

cat("Input data file:\n", data_file, "\n\n")
cat("Input covariance matrix:\n", V_file, "\n\n")
cat("Output directory:\n", output_dir, "\n\n")


# =============================================================================
# 3. Check input files
# =============================================================================

if (!file.exists(data_file)) {
  stop(
    paste0(
      "Input data file not found:\n",
      data_file
    )
  )
}

if (!file.exists(V_file)) {
  stop(
    paste0(
      "Covariance matrix file not found:\n",
      V_file
    )
  )
}


# =============================================================================
# 4. Read effect-size data and covariance matrix
# =============================================================================

cat("Step 1/8: Loading Nakagawa primary results...\n")

dat <- as.data.frame(
  readRDS(data_file)
)

V_full <- readRDS(
  V_file
)

N_MAIN <- NROW(dat)

cat(
  "Input rows:",
  N_MAIN,
  " rows\n"
)

cat(
  "V：",
  NROW(V_full),
  " × ",
  NCOL(V_full),
  "\n\n"
)


# =============================================================================
# 5. Validate data and covariance matrix
# =============================================================================

if (
  NROW(dat) != NROW(V_full) ||
  NROW(dat) != NCOL(V_full)
) {
  stop(
    "Data dimensions do not match the covariance matrix."
  )
}


# -----------------------------------------------------------------------------
# V_index
# -----------------------------------------------------------------------------

if (!("V_index" %in% names(dat))) {
  dat$V_index <- seq_len(
    NROW(dat)
  )
}


if (
  anyNA(dat$V_index) ||
  anyDuplicated(dat$V_index) > 0
) {
  stop(
    "V_index contains missing or duplicate entries."
  )
}


if (!setequal(
  dat$V_index,
  seq_len(NROW(dat))
)) {
  stop(
    "V_index is not a complete sequence from 1 to n."
  )
}


# =============================================================================
# 6. Required variables
# =============================================================================

required_cols <- c(
  "Trait",
  "Treatment_1",
  "Study_id",
  "Row_id",
  "yi",
  "vi",
  "V_index",
  "Sampling_season",
  "Types_of_frozen_soil",
  "StandAge_year",
  "Soil_depth"
)

missing_cols <- setdiff(
  required_cols,
  names(dat)
)

if (length(missing_cols) > 0) {
  stop(
    paste0(
      "Required columns are missing:",
      paste(
        missing_cols,
        collapse = ", "
      ),
      "\nAvailable columns:",
      paste(names(dat), collapse = " | ")
    )
  )
}


# =============================================================================
# 7. Normalize treatment groups
# =============================================================================

cat("Step 2/8: Standardizing treatment labels...\n")


dat$Treatment_raw <- trimws(
  as.character(
    dat$Treatment_1
  )
)


normalize_treatment <- function(x) {

  z <- trimws(
    tolower(
      as.character(x)
    )
  )

  out <- rep(
    NA_character_,
    length(z)
  )


  # ---------------------------------------------------------------------------
  # Snow addition
  # ---------------------------------------------------------------------------

  add_flag <- (
    grepl("increase", z) |
    grepl("addition", z) |
    grepl("added", z) |
    grepl("add snow", z) |
    grepl("\u589e\u96ea", z) |
    grepl("\u52a0\u96ea", z)
  )


  # ---------------------------------------------------------------------------
  # Snow removal
  # ---------------------------------------------------------------------------

  remove_flag <- (
    grepl("decrease", z) |
    grepl("reduced", z) |
    grepl("reduction", z) |
    grepl("remove", z) |
    grepl("removal", z) |
    grepl("\u51cf\u96ea", z) |
    grepl("\u9664\u96ea", z)
  )

  add_flag[is.na(add_flag)] <- FALSE
  remove_flag[is.na(remove_flag)] <- FALSE


  out[
    add_flag &
    !remove_flag
  ] <- "Snow addition"


  out[
    remove_flag &
    !add_flag
  ] <- "Snow removal"


  return(out)
}


dat$Treatment_group <- normalize_treatment(
  dat$Treatment_raw
)


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

unmatched_treatment <- sort(
  unique(
    dat$Treatment_raw[
      is.na(
        dat$Treatment_group
      )
    ]
  )
)


if (length(unmatched_treatment) > 0) {

  cat("\nUnrecognized Treatment_1 labels:\n")

  print(
    unmatched_treatment
  )

  stop(
    "Unrecognized values in Treatment_1."
  )
}


dat$Treatment_group <- factor(
  dat$Treatment_group,
  levels = c(
    "Snow addition",
    "Snow removal"
  )
)


cat("\nTreatment counts:\n")

print(
  table(
    dat$Treatment_group
  )
)

cat("\n")


# =============================================================================
# 8. Normalize sampling season and frozen-soil groups
# =============================================================================

normalize_simple_key <- function(x) {
  z <- tolower(trimws(as.character(x)))
  z <- gsub("[[:space:]]+", "", z)
  z <- gsub("[‐‑‒–—−]", "-", z)
  z
}



# -----------------------------------------------------------------------------
# Sample season
# -----------------------------------------------------------------------------

normalize_season <- function(x) {

  z <- normalize_simple_key(x)
  z_no_dash <- gsub("-", "", z)

  out <- rep(NA_character_, length(z))

  out[!is.na(z_no_dash) & z_no_dash == "annualaverage"] <- "Annual average"

  out[!is.na(z_no_dash) & z_no_dash %in% c(
    "freezeseason",
    "freeze"
  )] <- "Freeze season"

  out[!is.na(z_no_dash) & z_no_dash %in% c(
    "freezethawseason",
    "freezethawing",
    "freezethaw",
    "freezethawingseason"
  )] <- "Freeze-thaw season"

  out[!is.na(z_no_dash) & z_no_dash == "growingseason"] <- "Growing season"

  out
}


dat$Sampling_season_group <- normalize_season(
  dat$Sampling_season
)


# -----------------------------------------------------------------------------
# Frozen soil types
# -----------------------------------------------------------------------------

normalize_frozen_soil <- function(x) {

  z <- normalize_simple_key(x)
  z_no_space <- gsub("[-_[:space:]]+", "", z)

  out <- rep(NA_character_, length(z))

  out[!is.na(z_no_space) & z_no_space %in% c(
    "permafrost",
    "p"
  )] <- "Permafrost"

  out[!is.na(z_no_space) & z_no_space %in% c(
    "seasonalfrostsoil",
    "seasonalfrost",
    "nonpermafrost",
    "np"
  )] <- "Seasonal frost soil"

  out
}


dat$Frozen_soil_group <- normalize_frozen_soil(
  dat$Types_of_frozen_soil
)


# =============================================================================
# 9. Classify experimental duration using StandAge_year
# =============================================================================
#   == 1       -> 1
#   > 1 & <= 5 -> 2-5
#   > 5 & <=10 -> 6-10
#   >10 & <=15 -> 11-15
#   >15 & <=25 -> 16-25

stand_age_num <- suppressWarnings(
  as.numeric(as.character(dat$StandAge_year))
)

dat$StandAge_group <- NA_character_

dat$StandAge_group[
  is.finite(stand_age_num) & stand_age_num == 1
] <- "1"

dat$StandAge_group[
  is.finite(stand_age_num) & stand_age_num > 1 & stand_age_num <= 5
] <- "2-5"

dat$StandAge_group[
  is.finite(stand_age_num) & stand_age_num > 5 & stand_age_num <= 10
] <- "6-10"

dat$StandAge_group[
  is.finite(stand_age_num) & stand_age_num > 10 & stand_age_num <= 15
] <- "11-15"

dat$StandAge_group[
  is.finite(stand_age_num) & stand_age_num > 15 & stand_age_num <= 25
] <- "16-25"

dat$StandAge_group <- factor(
  dat$StandAge_group,
  levels = c("1", "2-5", "6-10", "11-15", "16-25")
)

duration_label_map <- c(
  "1" = "≤1 year",
  "2-5" = "2-5 years",
  "6-10" = "6-10 years",
  "11-15" = "11-15 years",
  "16-25" = "16-25 years"
)

dat$Experimental_duration_group <- unname(
  duration_label_map[as.character(dat$StandAge_group)]
)

dat$Experimental_duration_group <- factor(
  dat$Experimental_duration_group,
  levels = c(
    "≤1 year",
    "2-5 years",
    "6-10 years",
    "11-15 years",
    "16-25 years"
  )
)

cat("\nExperimental-duration counts (original StandAge_year bins):\n")
print(table(dat$Experimental_duration_group, useNA = "ifany"))
cat("\n")


# =============================================================================
# 10. Classify soil depth
# =============================================================================
normalize_band_string <- function(x) {
  raw <- as.character(x)
  s <- tolower(trimws(raw))
  missing_flag <- is.na(raw) | is.na(s) | s == "" |
    s %in% c("na", "n/a", "null", "none", "missing", "nan")
  s <- gsub("[‐‑‒–—−]", "-", s)
  s <- gsub("[～〜~]", "-", s)
  s <- gsub("[\u81f3\u5230]", "-", s)
  s <- gsub("[[:space:]]+", "", s)
  s <- gsub("centimeters|centimeter|centimetres|centimetre|cm|\u5398\u7c73", "", s)
  s[missing_flag] <- NA_character_
  s
}



DEPTH_LEVELS <- c(
  "0-10 cm",
  "11-30 cm",
  "31-60 cm",
  "61-100 cm"
)


classify_depth_one <- function(x) {

  s <- normalize_band_string(x)

  if (is.na(s)) {
    return(NA_character_)
  }

  matched <- regmatches(
    s,
    regexec(
      "^([0-9]+(?:\\.[0-9]+)?)-([0-9]+(?:\\.[0-9]+)?)$",
      s,
      perl = TRUE
    )
  )[[1]]

  if (length(matched) == 3L) {

    low <- as.numeric(matched[2])
    high <- as.numeric(matched[3])

    if (
      !is.finite(low) ||
      !is.finite(high) ||
      low > high
    ) {
      return(NA_character_)
    }

    if (low >= 0 && high <= 10) {
      return(DEPTH_LEVELS[1])
    }

    if (low >= 10 && high <= 30) {
      return(DEPTH_LEVELS[2])
    }

    if (low >= 30 && high <= 60) {
      return(DEPTH_LEVELS[3])
    }

    if (low >= 60 && high <= 100) {
      return(DEPTH_LEVELS[4])
    }

    return(NA_character_)
  }

  value <- suppressWarnings(as.numeric(s))

  if (!is.finite(value) || value < 0) {
    return(NA_character_)
  }

  if (value <= 10) {
    return(DEPTH_LEVELS[1])
  }

  if (value <= 30) {
    return(DEPTH_LEVELS[2])
  }

  if (value <= 60) {
    return(DEPTH_LEVELS[3])
  }

  if (value <= 100) {
    return(DEPTH_LEVELS[4])
  }

  NA_character_
}


dat$Soil_depth_group <- vapply(
  dat$Soil_depth,
  classify_depth_one,
  character(1)
)

dat$Soil_depth_group <- factor(
  dat$Soil_depth_group,
  levels = DEPTH_LEVELS
)


# =============================================================================
# 11. Group-mapping quality control
# =============================================================================

make_group_mapping_qc <- function(
    raw_values,
    mapped_values,
    group_name,
    source_column
) {

  raw_label <- trimws(
    as.character(raw_values)
  )

  raw_label[
    is.na(raw_label) |
    raw_label == ""
  ] <- "<NA>"

  mapped_label <- as.character(
    mapped_values
  )

  raw_is_missing <- is.na(raw_values) |
    is.na(raw_label) |
    raw_label == "<NA>"

  mapped_is_missing <- is.na(mapped_label) |
    mapped_label == ""

  mapped_label[
    mapped_is_missing & raw_is_missing
  ] <- "<MISSING>"

  mapped_label[
    mapped_is_missing & !raw_is_missing
  ] <- "<UNMATCHED>"

  tmp <- data.frame(
    Group_type = group_name,
    Source_column = source_column,
    Original_value = raw_label,
    Mapped_value = mapped_label,
    stringsAsFactors = FALSE
  )

  tmp %>%
    dplyr::group_by(
      Group_type,
      Source_column,
      Original_value,
      Mapped_value
    ) %>%
    dplyr::summarise(
      n_records = dplyr::n(),
      .groups = "drop"
    )
}


group_mapping_qc <- dplyr::bind_rows(

  make_group_mapping_qc(
    dat$Sampling_season,
    dat$Sampling_season_group,
    "Sample season",
    "Sampling_season"
  ),

  make_group_mapping_qc(
    dat$Types_of_frozen_soil,
    dat$Frozen_soil_group,
    "Frozen soil types",
    "Types_of_frozen_soil"
  ),

  make_group_mapping_qc(
    dat$StandAge_year,
    dat$Experimental_duration_group,
    "Experimental duration",
    "StandAge_year -> original StandAge_group bins"
  ),

  make_group_mapping_qc(
    dat$Soil_depth,
    dat$Soil_depth_group,
    "Soil depth",
    "Soil_depth"
  )
)


# =============================================================================
# 12. Define prespecified groups
# =============================================================================

group_definitions <- list()


# -----------------------------------------------------------------------------
# Overall
# -----------------------------------------------------------------------------

group_definitions[["Overall"]] <- list(
  variable = NULL,
  levels = "Overall"
)



# -----------------------------------------------------------------------------
# Sample season
# -----------------------------------------------------------------------------

group_definitions[["Sample season"]] <- list(
  variable = "Sampling_season_group",
  levels = c(
    "Annual average",
    "Freeze season",
    "Freeze-thaw season",
    "Growing season"
  )
)


# -----------------------------------------------------------------------------
# Frozen soil types
# -----------------------------------------------------------------------------

group_definitions[["Frozen soil types"]] <- list(
  variable = "Frozen_soil_group",
  levels = c(
    "Permafrost",
    "Seasonal frost soil"
  )
)


# -----------------------------------------------------------------------------
# Experimental duration
# -----------------------------------------------------------------------------

group_definitions[["Experimental duration"]] <- list(
  variable = "Experimental_duration_group",
  levels = c(
    "≤1 year",
    "2-5 years",
    "6-10 years",
    "11-15 years",
    "16-25 years"
  )
)


# -----------------------------------------------------------------------------
# Soil depth
# -----------------------------------------------------------------------------

group_definitions[["Soil depth"]] <- list(
  variable = "Soil_depth_group",
  levels = c(
    "0-10 cm",
    "11-30 cm",
    "31-60 cm",
    "61-100 cm"
  )
)


group_level_order <- lapply(
  group_definitions,
  function(x) x$levels
)


cat(
  "Predefined subgroup levels:\n"
)

for (nm in names(group_definitions)) {
  cat(
    "  ",
    nm,
    ": ",
    paste(group_definitions[[nm]]$levels, collapse = " | "),
    "\n",
    sep = ""
  )
}

cat("\n")


# =============================================================================
# 13. Subset rows and V consistently
# =============================================================================

subset_data_V <- function(
    dat_input,
    V_input,
    row_positions
) {

  if (
    length(row_positions) == 0
  ) {

    d_sub <- dat_input[
      integer(0),
      ,
      drop = FALSE
    ]

    V_sub <- V_input[
      integer(0),
      integer(0),
      drop = FALSE
    ]

    return(
      list(
        data = d_sub,
        V = V_sub
      )
    )
  }


  d_sub <- dat_input[
    row_positions,
    ,
    drop = FALSE
  ]


  original_V_index <- as.integer(
    d_sub$V_index
  )


  V_sub <- V_input[
    original_V_index,
    original_V_index,
    drop = FALSE
  ]


  if (
    NROW(d_sub) != NROW(V_sub) ||
    NROW(d_sub) != NCOL(V_sub)
  ) {

    stop(
      "Subset data and covariance matrix dimensions do not match."
    )
  }


  return(
    list(
      data = d_sub,
      V = V_sub
    )
  )
}


# =============================================================================
# 14. Estimate treatment-specific effects and QE
#
#
# =============================================================================

extract_QE <- function(fit) {
  qe <- if (length(fit$QE) == 1L) as.numeric(fit$QE) else NA_real_
  qe_df <- if (length(fit$QEdf) == 1L) {
    as.numeric(fit$QEdf)
  } else if (length(fit$k) == 1L && length(fit$p) == 1L) {
    as.numeric(fit$k - fit$p)
  } else NA_real_
  qe_p <- if (length(fit$QEp) == 1L) as.numeric(fit$QEp) else NA_real_
  list(QE = qe, QE_df = qe_df, QE_p = qe_p)
}

fit_one_treatment <- function(
    dat_sub,
    V_sub,
    treatment_name,
    trait_name,
    group_type,
    group_level
) {

  local_idx <- which(
    as.character(
      dat_sub$Treatment_group
    ) ==
    treatment_name
  )


  if (length(local_idx) == 0) {

    return(
      data.frame(
        Trait = trait_name,
        Group_type = group_type,
        Group_level = group_level,
        Treatment = treatment_name,
        Estimate = NA_real_,
        SE = NA_real_,
        CI_lower = NA_real_,
        CI_upper = NA_real_,
        p_value = NA_real_,
        QE = NA_real_, QE_df = NA_real_, QE_p = NA_real_,
        n = 0L,
        Study_n = 0L,
        Model_status = "No data",
        stringsAsFactors = FALSE
      )
    )
  }


  d_trt <- dat_sub[
    local_idx,
    ,
    drop = FALSE
  ]


  V_trt <- V_sub[
    local_idx,
    local_idx,
    drop = FALSE
  ]


  n_effect <- NROW(
    d_trt
  )


  n_study <- length(
    unique(
      d_trt$Study_id
    )
  )


  if (n_effect < MIN_EFFECTS_PER_TREATMENT) {
    return(data.frame(
      Trait = trait_name, Group_type = group_type, Group_level = group_level,
      Treatment = treatment_name, Estimate = NA_real_, SE = NA_real_,
      CI_lower = NA_real_, CI_upper = NA_real_, p_value = NA_real_,
      QE = NA_real_, QE_df = NA_real_, QE_p = NA_real_,
      n = n_effect, Study_n = n_study,
      Model_status = paste0("Insufficient effects (n < ", MIN_EFFECTS_PER_TREATMENT, ")"),
      stringsAsFactors = FALSE
    ))
  }

  fit <- tryCatch(

    metafor::rma.mv(

      yi = yi,

      V = V_trt,

      random =
        ~ 1 | Study_id / Row_id,

      data = d_trt,

      method = "REML",

      sparse = TRUE
    ),

    error = function(e) {

      message(
        paste0(
          "Treatment-specific model failed:",
          trait_name,
          " | ",
          group_type,
          " | ",
          group_level,
          " | ",
          treatment_name,
          " | ",
          conditionMessage(e)
        )
      )

      NULL
    }
  )


  if (is.null(fit)) {

    return(
      data.frame(
        Trait = trait_name,
        Group_type = group_type,
        Group_level = group_level,
        Treatment = treatment_name,
        Estimate = NA_real_,
        SE = NA_real_,
        CI_lower = NA_real_,
        CI_upper = NA_real_,
        p_value = NA_real_,
        QE = NA_real_, QE_df = NA_real_, QE_p = NA_real_,
        n = n_effect,
        Study_n = n_study,
        Model_status = "Model failed",
        stringsAsFactors = FALSE
      )
    )
  }


  qe_one <- extract_QE(fit)

  return(
    data.frame(
      Trait = trait_name,
      Group_type = group_type,
      Group_level = group_level,
      Treatment = treatment_name,

      Estimate =
        as.numeric(
          fit$b[1]
        ),

      SE =
        as.numeric(
          fit$se[1]
        ),

      CI_lower =
        as.numeric(
          fit$ci.lb[1]
        ),

      CI_upper =
        as.numeric(
          fit$ci.ub[1]
        ),

      p_value =
        as.numeric(
          fit$pval[1]
        ),

      QE = qe_one$QE,
      QE_df = qe_one$QE_df,
      QE_p = qe_one$QE_p,

      n =
        n_effect,

      Study_n =
        n_study,

      Model_status =
        "rma.mv",

      stringsAsFactors =
        FALSE
    )
  )
}


# =============================================================================
# 15. Test treatment difference (QM) and residual QE
#
# Reference：
#
# Snow addition
#
# coefficient：
#
# Snow removal - Snow addition
# =============================================================================

fit_treatment_QM <- function(
    dat_sub,
    V_sub,
    trait_name,
    group_type,
    group_level
) {


  trt_char <- as.character(
    dat_sub$Treatment_group
  )


  add_idx <- which(
    trt_char ==
    "Snow addition"
  )


  remove_idx <- which(
    trt_char ==
    "Snow removal"
  )


  n_add <- length(
    add_idx
  )


  n_remove <- length(
    remove_idx
  )


  study_add <- length(
    unique(
      dat_sub$Study_id[
        add_idx
      ]
    )
  )


  study_remove <- length(
    unique(
      dat_sub$Study_id[
        remove_idx
      ]
    )
  )


  n_total <- NROW(
    dat_sub
  )


  study_total <- length(
    unique(
      dat_sub$Study_id
    )
  )


  # ===========================================================================
  # ===========================================================================

  if (
    n_total == 0 ||
    n_add < MIN_EFFECTS_PER_TREATMENT ||
    n_remove < MIN_EFFECTS_PER_TREATMENT
  ) {

    status_missing <- if (n_total == 0) {
      "No data"
    } else if (n_add == 0 || n_remove == 0) {
      "Only one treatment"
    } else {
      "Insufficient effects per treatment (n < 3)"
    }

    return(
      data.frame(

        Trait =
          trait_name,

        Group_type =
          group_type,

        Group_level =
          group_level,

        n_total =
          n_total,

        Study_n_total =
          study_total,

        n_add =
          n_add,

        Study_n_add =
          study_add,

        n_remove =
          n_remove,

        Study_n_remove =
          study_remove,

        Beta_remove_minus_add =
          NA_real_,

        SE_difference =
          NA_real_,

        CI_lower_difference =
          NA_real_,

        CI_upper_difference =
          NA_real_,

        p_difference =
          NA_real_,

        QE = NA_real_, QE_df = NA_real_, QE_p = NA_real_,

        QM =
          NA_real_,

        QM_df =
          NA_real_,

        QM_p =
          NA_real_,

        QM_status =
          status_missing,

        stringsAsFactors =
          FALSE
      )
    )
  }


  # ===========================================================================
  # ===========================================================================

  dat_qm <- dat_sub


  dat_qm$Treatment_group <- factor(

    as.character(
      dat_qm$Treatment_group
    ),

    levels = c(
      "Snow addition",
      "Snow removal"
    )
  )


  fit_qm <- tryCatch(

    metafor::rma.mv(

      yi = yi,

      V = V_sub,

      mods =
        ~ Treatment_group,

      random =
        ~ 1 | Study_id / Row_id,

      data = dat_qm,

      method = "REML",

      sparse = TRUE
    ),

    error = function(e) {

      message(
        paste0(
          "Treatment QM model failed:",
          trait_name,
          " | ",
          group_type,
          " | ",
          group_level,
          " | ",
          conditionMessage(e)
        )
      )

      NULL
    }
  )


  if (is.null(fit_qm)) {

    return(
      data.frame(

        Trait = trait_name,
        Group_type = group_type,
        Group_level = group_level,

        n_total = n_total,
        Study_n_total = study_total,

        n_add = n_add,
        Study_n_add = study_add,

        n_remove = n_remove,
        Study_n_remove = study_remove,

        Beta_remove_minus_add = NA_real_,
        SE_difference = NA_real_,
        CI_lower_difference = NA_real_,
        CI_upper_difference = NA_real_,
        p_difference = NA_real_,

        QE = NA_real_, QE_df = NA_real_, QE_p = NA_real_,
        QM = NA_real_,
        QM_df = NA_real_,
        QM_p = NA_real_,

        QM_status = "Model failed",

        stringsAsFactors = FALSE
      )
    )
  }


  # ===========================================================================
  # ===========================================================================

  coef_names <- rownames(
    fit_qm$b
  )


  coef_idx <- grep(
    "Treatment_group",
    coef_names,
    fixed = TRUE
  )


  if (length(coef_idx) != 1) {

    if (NROW(fit_qm$b) >= 2) {

      coef_idx <- 2

    } else {

      coef_idx <- NA_integer_
    }
  }


  if (is.na(coef_idx)) {

    beta_diff <- NA_real_
    se_diff <- NA_real_
    ci_lb_diff <- NA_real_
    ci_ub_diff <- NA_real_
    p_diff <- NA_real_

  } else {

    beta_diff <- as.numeric(
      fit_qm$b[
        coef_idx
      ]
    )

    se_diff <- as.numeric(
      fit_qm$se[
        coef_idx
      ]
    )

    ci_lb_diff <- as.numeric(
      fit_qm$ci.lb[
        coef_idx
      ]
    )

    ci_ub_diff <- as.numeric(
      fit_qm$ci.ub[
        coef_idx
      ]
    )

    p_diff <- as.numeric(
      fit_qm$pval[
        coef_idx
      ]
    )
  }


  QM_value <- if (
    length(
      fit_qm$QM
    ) == 1
  ) {

    as.numeric(
      fit_qm$QM
    )

  } else {

    NA_real_
  }


  QM_p_value <- if (
    length(
      fit_qm$QMp
    ) == 1
  ) {

    as.numeric(
      fit_qm$QMp
    )

  } else {

    p_diff
  }


  qe_joint <- extract_QE(fit_qm)

  return(
    data.frame(

      Trait =
        trait_name,

      Group_type =
        group_type,

      Group_level =
        group_level,

      n_total =
        n_total,

      Study_n_total =
        study_total,

      n_add =
        n_add,

      Study_n_add =
        study_add,

      n_remove =
        n_remove,

      Study_n_remove =
        study_remove,

      Beta_remove_minus_add =
        beta_diff,

      SE_difference =
        se_diff,

      CI_lower_difference =
        ci_lb_diff,

      CI_upper_difference =
        ci_ub_diff,

      p_difference =
        p_diff,

      QE = qe_joint$QE,
      QE_df = qe_joint$QE_df,
      QE_p = qe_joint$QE_p,

      QM =
        QM_value,

      QM_df =
        1,

      QM_p =
        QM_p_value,

      QM_status =
        "Calculated",

      stringsAsFactors =
        FALSE
    )
  )
}


# =============================================================================
# 16. Analyze one trait-subgroup combination
# =============================================================================

analyse_one_subgroup <- function(
    dat_sub,
    V_sub,
    trait_name,
    group_type,
    group_level
) {


  add_result <- fit_one_treatment(

    dat_sub = dat_sub,

    V_sub = V_sub,

    treatment_name =
      "Snow addition",

    trait_name =
      trait_name,

    group_type =
      group_type,

    group_level =
      group_level
  )


  remove_result <- fit_one_treatment(

    dat_sub = dat_sub,

    V_sub = V_sub,

    treatment_name =
      "Snow removal",

    trait_name =
      trait_name,

    group_type =
      group_type,

    group_level =
      group_level
  )


  effect_result <- rbind(

    add_result,

    remove_result
  )


  qm_result <- fit_treatment_QM(

    dat_sub = dat_sub,

    V_sub = V_sub,

    trait_name =
      trait_name,

    group_type =
      group_type,

    group_level =
      group_level
  )


  return(
    list(
      effect = effect_result,
      qm = qm_result
    )
  )
}


# =============================================================================
# 17. Loop over traits and groups
# =============================================================================

cat(
  "Step 3/8: Analyzing traits and subgroups...\n"
)


traits <- sort(
  unique(
    as.character(
      dat$Trait
    )
  )
)


effect_result_list <- list()

qm_result_list <- list()


effect_counter <- 0

qm_counter <- 0


for (trait_name in traits) {


  cat("\n")
  cat(strrep("-", 80), "\n")

  cat(
    "Trait：",
    trait_name,
    "\n"
  )

  cat(strrep("-", 80), "\n")


  # ===========================================================================
  # Trait rows
  # ===========================================================================

  trait_positions <- which(
    as.character(
      dat$Trait
    ) ==
    trait_name
  )


  if (length(trait_positions) == 0) {
    next
  }


  # ===========================================================================
  # ===========================================================================

  for (group_type in names(group_definitions)) {


    group_info <- group_definitions[[group_type]]


    group_var <- group_info$variable

    group_levels <- group_info$levels


    # =========================================================================
    # Overall
    # =========================================================================

    if (is.null(group_var)) {


      sub_obj <- subset_data_V(

        dat_input =
          dat,

        V_input =
          V_full,

        row_positions =
          trait_positions
      )


      if (is.null(sub_obj)) {
        next
      }


      result_one <- analyse_one_subgroup(

        dat_sub =
          sub_obj$data,

        V_sub =
          sub_obj$V,

        trait_name =
          trait_name,

        group_type =
          "Overall",

        group_level =
          "Overall"
      )


      effect_counter <- effect_counter + 1

      effect_result_list[[effect_counter]] <-
        result_one$effect


      qm_counter <- qm_counter + 1

      qm_result_list[[qm_counter]] <-
        result_one$qm


      next
    }


    # =========================================================================
    # =========================================================================

    if (!(group_var %in% names(dat))) {

      message(
        paste0(
          "Skipping missing subgroup column:",
          group_var
        )
      )

      next
    }


    group_values <- as.character(
      dat[[group_var]]
    )


    for (group_level in group_levels) {


      row_positions <- which(

        as.character(
          dat$Trait
        ) ==
          trait_name &

        !is.na(
          dat[[group_var]]
        ) &

        group_values ==
          as.character(
            group_level
          )
      )



      cat(
        "  ",
        group_type,
        " | ",
        group_level,
        " | n=",
        length(row_positions),
        "\n",
        sep = ""
      )


      sub_obj <- subset_data_V(

        dat_input =
          dat,

        V_input =
          V_full,

        row_positions =
          row_positions
      )


      if (is.null(sub_obj)) {
        next
      }


      result_one <- analyse_one_subgroup(

        dat_sub =
          sub_obj$data,

        V_sub =
          sub_obj$V,

        trait_name =
          trait_name,

        group_type =
          group_type,

        group_level =
          as.character(
            group_level
          )
      )


      effect_counter <- effect_counter + 1

      effect_result_list[[effect_counter]] <-
        result_one$effect


      qm_counter <- qm_counter + 1

      qm_result_list[[qm_counter]] <-
        result_one$qm
    }
  }
}


# =============================================================================
# 18. Assemble results
# =============================================================================

cat(
  "\nStep 4/8: Combining results...\n"
)


if (length(effect_result_list) > 0) {

  effect_results <- do.call(
    rbind,
    effect_result_list
  )

} else {

  effect_results <- data.frame()
}


if (length(qm_result_list) > 0) {

  qm_results <- do.call(
    rbind,
    qm_result_list
  )

} else {

  qm_results <- data.frame()
}


rownames(
  effect_results
) <- NULL


rownames(
  qm_results
) <- NULL


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

if (NROW(effect_results) > 0) {

  effect_results$Treatment <- dplyr::recode(
    as.character(effect_results$Treatment),
    "Snow addition" = "Increased snowpack thickness",
    "Snow removal" = "Decreased snowpack thickness"
  )
}


# =============================================================================
# 19. Significance flags
# =============================================================================

cat(
  "Step 5/8: Adding significance indicators...\n"
)


if (NROW(effect_results) > 0) {

  effect_results$Significant <- ifelse(
    !is.finite(effect_results$p_value),
    NA_character_,
    ifelse(
      effect_results$p_value < 0.05,
      "Yes",
      "No"
    )
  )
}


if (NROW(effect_results) > 0) {
  effect_results$QE_significant <- ifelse(
    is.finite(effect_results$QE_p),
    ifelse(effect_results$QE_p < 0.05, "Yes", "No"), NA_character_
  )
}

if (NROW(qm_results) > 0) {
  qm_results$QE_significant <- ifelse(
    is.finite(qm_results$QE_p),
    ifelse(qm_results$QE_p < 0.05, "Yes", "No"), NA_character_
  )
}

if (NROW(qm_results) > 0) {

  qm_results$QM_significant <- ifelse(
    !is.finite(qm_results$QM_p),
    NA_character_,
    ifelse(
      qm_results$QM_p < 0.05,
      "Yes",
      "No"
    )
  )
}


# =============================================================================
# 20. Sort results
# =============================================================================

group_order <- c(
  "Overall",
  "Sample season",
  "Frozen soil types",
  "Experimental duration",
  "Soil depth"
)


treatment_order <- c(
  "Increased snowpack thickness",
  "Decreased snowpack thickness"
)


if (NROW(effect_results) > 0) {

  effect_results$Group_type_sort <- match(
    effect_results$Group_type,
    group_order
  )


  effect_results$Treatment_sort <- match(
    effect_results$Treatment,
    treatment_order
  )

  effect_results$Group_level_sort <- mapply(
    function(gt, gl) {
      match(
        gl,
        group_level_order[[gt]]
      )
    },
    as.character(effect_results$Group_type),
    as.character(effect_results$Group_level)
  )


  effect_results <- effect_results %>%

    arrange(
      Trait,
      Group_type_sort,
      Group_level_sort,
      Treatment_sort
    )


  effect_results$Group_type_sort <- NULL
  effect_results$Group_level_sort <- NULL
  effect_results$Treatment_sort <- NULL
}


if (NROW(qm_results) > 0) {

  qm_results$Group_type_sort <- match(
    qm_results$Group_type,
    group_order
  )


  qm_results$Group_level_sort <- mapply(
    function(gt, gl) {
      match(
        gl,
        group_level_order[[gt]]
      )
    },
    as.character(qm_results$Group_type),
    as.character(qm_results$Group_level)
  )


  qm_results <- qm_results %>%

    arrange(
      Trait,
      Group_type_sort,
      Group_level_sort
    )


  qm_results$Group_type_sort <- NULL
  qm_results$Group_level_sort <- NULL
}


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

if (NROW(qm_results) > 0) {

  qm_results <- qm_results %>%
    dplyr::rename(
      n_increased = n_add,
      Study_n_increased = Study_n_add,
      n_decreased = n_remove,
      Study_n_decreased = Study_n_remove,
      Beta_decreased_minus_increased = Beta_remove_minus_add
    )
}


# =============================================================================
# 21. Combine treatment effects and tests
# =============================================================================

cat(
  "Step 6/8: Creating the combined table...\n"
)


if (
  NROW(effect_results) > 0 &&
  NROW(qm_results) > 0
) {

  combined_results <- effect_results %>%

    left_join(

      qm_results %>%

        select(

          Trait,
          Group_type,
          Group_level,

          n_total,
          Study_n_total,

          n_increased,
          Study_n_increased,

          n_decreased,
          Study_n_decreased,

          Beta_decreased_minus_increased,
          SE_difference,
          CI_lower_difference,
          CI_upper_difference,
          p_difference,

          QE_joint = QE,
          QE_df_joint = QE_df,
          QE_p_joint = QE_p,
          QE_significant_joint = QE_significant,
          QM,
          QM_df,
          QM_p,
          QM_significant,
          QM_status
        ),

      by = c(
        "Trait",
        "Group_type",
        "Group_level"
      )
    )

} else {

  combined_results <- effect_results
}


# =============================================================================
# 22. Significant QM results
# =============================================================================

if (NROW(qm_results) > 0) {

  significant_qm <- qm_results %>%

    filter(
      is.finite(QM_p),
      QM_p < 0.05
    )

} else {

  significant_qm <- data.frame()
}


# =============================================================================
# 23. Significant treatment effect results
# =============================================================================

if (NROW(effect_results) > 0) {

  significant_effects <- effect_results %>%

    filter(
      is.finite(p_value),
      p_value < 0.05
    )

} else {

  significant_effects <- data.frame()
}


# =============================================================================
# =============================================================================
significant_QE_effects <- effect_results %>%
  dplyr::filter(is.finite(QE_p), QE_p < 0.05)
significant_QE_joint <- qm_results %>%
  dplyr::filter(is.finite(QE_p), QE_p < 0.05)

# =============================================================================
# =============================================================================

RUN_GROUP_MODERATOR_QM <- TRUE

fit_group_moderator_QM <- function(trait_name, treatment_name, group_type) {
  info <- group_definitions[[group_type]]
  group_var <- info$variable
  predefined <- info$levels
  positions <- which(
    as.character(dat$Trait) == trait_name &
    as.character(dat$Treatment_group) == treatment_name &
    !is.na(dat[[group_var]]) &
    as.character(dat[[group_var]]) %in% predefined
  )
  group_values <- as.character(dat[[group_var]][positions])
  counts <- table(factor(group_values, levels = predefined))
  retained <- names(counts)[counts >= MIN_EFFECTS_PER_TREATMENT]
  kept_positions <- positions[group_values %in% retained]
  status <- NULL
  if (length(retained) < 2L) {
    status <- "Fewer than two eligible subgroup levels"
  }
  result <- data.frame(
    Trait = trait_name, Treatment = treatment_name,
    Group_type = group_type,
    Levels_used = paste(retained, collapse = " | "),
    n_all_mapped = length(positions), n_included = length(kept_positions),
    Study_n_included = length(unique(dat$Study_id[kept_positions])),
    QM_group = NA_real_, QM_group_df = NA_real_, QM_group_p = NA_real_,
    QE_group_model = NA_real_, QE_group_df = NA_real_, QE_group_p = NA_real_,
    Status = if (is.null(status)) "Pending" else status,
    stringsAsFactors = FALSE
  )
  if (!is.null(status)) return(result)
  obj <- subset_data_V(dat, V_full, kept_positions)
  d_mod <- obj$data
  d_mod$Group_mod <- factor(as.character(d_mod[[group_var]]), levels = retained)
  fit_mod <- tryCatch(
    metafor::rma.mv(
      yi = yi, V = obj$V, mods = ~ Group_mod,
      random = ~ 1 | Study_id / Row_id,
      data = d_mod, method = "REML", sparse = TRUE
    ),
    error = function(e) {
      message(paste0("Group moderator QM model failed:", trait_name, " | ",
                     treatment_name, " | ", group_type, " | ", conditionMessage(e)))
      NULL
    }
  )
  if (is.null(fit_mod)) {
    result$Status <- "Model fit failed"
    return(result)
  }
  qe_mod <- extract_QE(fit_mod)
  result$QM_group <- if (length(fit_mod$QM) == 1L) as.numeric(fit_mod$QM) else NA_real_
  result$QM_group_df <- length(retained) - 1L
  result$QM_group_p <- if (length(fit_mod$QMp) == 1L) as.numeric(fit_mod$QMp) else NA_real_
  result$QE_group_model <- qe_mod$QE
  result$QE_group_df <- qe_mod$QE_df
  result$QE_group_p <- qe_mod$QE_p
  result$Status <- "Calculated"
  result
}

if (RUN_GROUP_MODERATOR_QM) {
  group_mod_result_list <- list()
  counter_group_mod <- 0L
  for (trait_name in traits) {
    for (treatment_name in c("Snow addition", "Snow removal")) {
      for (gt in setdiff(names(group_definitions), "Overall")) {
        counter_group_mod <- counter_group_mod + 1L
        group_mod_result_list[[counter_group_mod]] <- fit_group_moderator_QM(
          trait_name, treatment_name, gt
        )
      }
    }
  }
  group_moderator_results <- dplyr::bind_rows(group_mod_result_list)
  group_moderator_results$Treatment <- dplyr::recode(
    group_moderator_results$Treatment,
    "Snow addition" = "Increased snowpack thickness",
    "Snow removal" = "Decreased snowpack thickness"
  )
  group_moderator_results$QM_group_significant <- ifelse(
    is.finite(group_moderator_results$QM_group_p),
    ifelse(group_moderator_results$QM_group_p < 0.05, "Yes", "No"), NA_character_
  )
} else {
  group_moderator_results <- data.frame(
    Note = "Group moderator QM analysis was disabled", stringsAsFactors = FALSE
  )
}

# =============================================================================
# 24. Coverage and QC
# =============================================================================

cat(
  "Step 7/8: Creating quality-control tables...\n"
)


group_coverage <- qm_results %>%
  dplyr::select(
    Trait,
    Group_type,
    Group_level,
    n_total,
    Study_n_total,
    n_increased,
    Study_n_increased,
    n_decreased,
    Study_n_decreased,
    QM_status
  )


unmapped_group_rows <- sum(
  group_mapping_qc$n_records[
    group_mapping_qc$Mapped_value == "<UNMATCHED>"
  ],
  na.rm = TRUE
)

missing_group_rows <- sum(
  group_mapping_qc$n_records[
    group_mapping_qc$Mapped_value == "<MISSING>"
  ],
  na.rm = TRUE
)


no_data_combinations <- sum(
  qm_results$n_total == 0,
  na.rm = TRUE
)


qc_summary <- data.frame(

  Item = c(
    "Main_data_rows",
    "V_rows",
    "V_columns",
    "Trait_number",
    "Minimum_effects_per_treatment",
    "Effect_result_rows",
    "QM_result_rows",
    "Significant_effect_rows",
    "Significant_QM_rows",
    "Significant_QE_effect_rows",
    "Significant_QE_joint_rows",
    "Group_moderator_QM_rows",
    "Expected_Trait_Group_combinations",
    "No_data_Trait_Group_combinations",
    "Experimental_duration_grouping_source",
    "Unmapped_group_records_across_variables",
    "Missing_group_records_across_variables"
  ),

  Value = c(
    NROW(dat),
    NROW(V_full),
    NCOL(V_full),
    length(traits),
    MIN_EFFECTS_PER_TREATMENT,
    NROW(effect_results),
    NROW(qm_results),
    NROW(significant_effects),
    NROW(significant_qm),
    NROW(significant_QE_effects),
    NROW(significant_QE_joint),
    NROW(group_moderator_results),
    NROW(qm_results),
    no_data_combinations,
    "Yes (StandAge_year; exact first-version StandAge_group bins)",
    unmapped_group_rows,
    missing_group_rows
  ),

  stringsAsFactors = FALSE
)


# =============================================================================
# 25. Export workbook
# =============================================================================

cat(
  "Step 8/8: Exporting results...\n"
)


writexl::write_xlsx(

  list(

    Combined =
      combined_results,

    Effect_sizes =
      effect_results,

    QM_treatment_difference =
      qm_results,

    Significant_effects =
      significant_effects,

    Significant_QM =
      significant_qm,

    Significant_QE_effects = significant_QE_effects,
    Significant_QE_joint = significant_QE_joint,
    QM_group_moderator = group_moderator_results,

    Group_coverage =
      group_coverage,

    Group_mapping_QC =
      group_mapping_qc,

    QC_summary =
      qc_summary

  ),

  output_excel
)


# =============================================================================
# 26. Save RDS
# =============================================================================

saveRDS(

  list(

    Effect_sizes =
      effect_results,

    QM =
      qm_results,

    Combined =
      combined_results,

    Significant_effects =
      significant_effects,

    Significant_QM =
      significant_qm,

    Significant_QE_effects = significant_QE_effects,
    Significant_QE_joint = significant_QE_joint,
    QM_group_moderator = group_moderator_results,

    Group_coverage =
      group_coverage,

    Group_mapping_QC =
      group_mapping_qc,

    QC =
      qc_summary

  ),

  output_rds
)


# =============================================================================
# 27. Summary
# =============================================================================

cat("\n")
cat(strrep("=", 90), "\n")
cat("Analysis completed.\n")
cat(strrep("=", 90), "\n\n")


cat(
  "Input data rows:",
  NROW(dat),
  "\n"
)


cat(
  "Number of traits:",
  length(traits),
  "\n"
)


cat(
  "Treatment-effect rows:",
  NROW(effect_results),
  " rows\n"
)


cat(
  "QM rows:",
  NROW(qm_results),
  " rows\n"
)


cat(
  "Significant treatment-effect rows:",
  NROW(significant_effects),
  " rows\n"
)


cat(
  "Significant treatment QM rows:",
  NROW(significant_qm),
  " rows\n"
)


cat(
  "Trait-subgroup combinations without data:",
  no_data_combinations,
  "\n"
)


cat(
  "Unmatched subgroup records (across variables):",
  unmapped_group_rows,
  "\n"
)


cat(
  "Missing subgroup records (across variables):",
  missing_group_rows,
  "\n\n"
)


cat(
  "Excel：\n",
  output_excel,
  "\n\n"
)


cat(
  "RDS：\n",
  output_rds,
  "\n\n"
)


cat(strrep("=", 90), "\n")
