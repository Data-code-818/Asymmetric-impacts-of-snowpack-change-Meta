# Nakagawa shared-control SMD variance-covariance matrix
#
# Shared-control blocks are defined by Study_id + Trait + Common_id.
# The first record in each block supplies the shared control mean, SD, and sample size.
# Control statistics are not used to define or split shared-control blocks.
#
# For each block:
#   k    = number of treatment groups
#   N    = n_control + sum(n_treatment)
#   df   = N - (k + 1)
#   sp^2 = [(nc - 1) * sdC^2 + sum((nt_i - 1) * sdT_i^2)] / df
#   d_i  = (meanT_i - meanC) / sqrt(sp^2)
#   J    = Gamma(df/2) / [sqrt(df/2) * Gamma((df - 1)/2)]
#   g_i  = J * d_i
#   V_ij = J^2 * [1/nc + d_i * d_j/(2*N)]
#   V_ii = V_ii + J^2/nt_i
#
# Place the input workbook in ./data or update input_file below.

# 0. Environment

rm(list = ls())

gc()

# 1. Packages

required_packages <- c(
  "readxl",
  "metafor",
  "Matrix",
  "writexl"
)

new_packages <- required_packages[
  !(required_packages %in% rownames(installed.packages()))
]

if (length(new_packages) > 0) {

  install.packages(
    new_packages
  )
}

library(readxl)

library(metafor)

library(Matrix)

library(writexl)

# 2. File paths

input_file <- file.path("data", "Meta_data_2026-10-30_renamed.xlsx")

output_dir <- "results"

if (!dir.exists(output_dir)) {

  dir.create(
    output_dir,
    recursive = TRUE
  )
}

prefix <- file.path(
  output_dir,
  "Meta_effect_sizes_2026-10-29"
)

cat("\n")

cat(
  strrep(
    "=",
    90
  ),
  "\n"
)

cat(
  "Nakagawa shared-control SMD VCV\n"
)

cat(
  strrep(
    "=",
    90
  ),
  "\n\n"
)

cat(
  "Input file:\n",
  input_file,
  "\n\n"
)

cat(
  "Output directory:\n",
  output_dir,
  "\n\n"
)

# 3. Exact Hedges' J correction
# Use the exact gamma-function correction.

hedges_J_exact <- function(df) {

  if (
    length(df) != 1 ||
    !is.finite(df) ||
    df <= 1
  ) {

    return(
      NA_real_
    )
  }

  J <- exp(

    lgamma(
      df / 2
    ) -

    0.5 *
    log(
      df / 2
    ) -

    lgamma(
      (df - 1) / 2
    )
  )

  return(
    J
  )
}

# 4. Input-file validation

if (!file.exists(input_file)) {

  stop(
    paste0(
      "Input file not found:\n",
      input_file
    )
  )
}

# 5. Data import

cat(
  "Step 1/12: Importing the data...\n"
)

raw_dat <- as.data.frame(

  readxl::read_excel(
    input_file
  )
)

RAW_ROWS <- NROW(
  raw_dat
)

if (
  length(RAW_ROWS) != 1 ||
  is.na(RAW_ROWS) ||
  RAW_ROWS <= 0
) {

  stop(
    "Failed to read the Excel workbook."
  )
}

# Preserve original row indices

raw_dat$original_row <- seq_len(
  RAW_ROWS
)

cat(
  "Number of input rows: ",
  RAW_ROWS,
  "\n"
)

if (RAW_ROWS == 7183) {

  cat(
    "Input row count matches the expected 7,183 records.\n\n"
  )

} else {

  warning(

    paste0(

      "Read ",
      RAW_ROWS,

      " rows; ",

      "the expected count is 7,183 after excluding three records."
    )
  )
}

# 6. Required-column validation

required_cols <- c(

  "Study_id",

  "Common_id",

  "Row_id",

  "Trait",

  "mean_treatment",

  "sd_treatment",

  "n_treatment",

  "mean_control",

  "sd_control",

  "n_control"
)

missing_cols <- setdiff(

  required_cols,

  names(
    raw_dat
  )
)

if (length(missing_cols) > 0) {

  stop(

    paste0(

      "Missing required columns: ",

      paste(
        missing_cols,
        collapse = ", "
      )
    )
  )
}

# 7. Numeric conversion

numeric_cols <- c(

  "mean_treatment",

  "sd_treatment",

  "n_treatment",

  "mean_control",

  "sd_control",

  "n_control"
)

for (v in numeric_cols) {

  raw_dat[[v]] <- suppressWarnings(

    as.numeric(

      as.character(
        raw_dat[[v]]
      )
    )
  )
}

# 8. Data validation

cat(
  "Step 2/12: Validating the input data...\n"
)

required_complete <- c(

  "Study_id",

  "Common_id",

  "Row_id",

  "Trait",

  "mean_treatment",

  "sd_treatment",

  "n_treatment",

  "mean_control",

  "sd_control",

  "n_control"
)

missing_rows <- which(

  !complete.cases(

    raw_dat[
      ,
      required_complete,
      drop = FALSE
    ]
  )
)

if (length(missing_rows) > 0) {

  stop(

    paste0(

      "Found ",

      length(
        missing_rows
      ),

      " row(s) with missing values in required fields."
    )
  )
}

invalid_rows <- which(

  raw_dat$n_treatment <= 0 |

  raw_dat$n_control <= 0 |

  raw_dat$sd_treatment <= 0 |

  raw_dat$sd_control <= 0
)

if (length(invalid_rows) > 0) {

  stop(

    paste0(

      "Found ",

      length(
        invalid_rows
      ),

      " row(s) with non-positive sample sizes or standard deviations."
    )
  )
}

# 9. Define shared-control blocks

cat(
  "Step 3/12: Defining shared-control blocks...\n"
)

raw_dat$block <- interaction(

  raw_dat$Study_id,

  raw_dat$Trait,

  raw_dat$Common_id,

  drop = TRUE,

  lex.order = TRUE
)

groups <- split(

  seq_len(
    RAW_ROWS
  ),

  raw_dat$block
)

N_BLOCKS <- length(
  groups
)

block_sizes <- lengths(
  groups
)

N_SHARED_BLOCKS <- sum(
  block_sizes >= 2
)

N_EFFECTS_SHARED <- sum(

  block_sizes[
    block_sizes >= 2
  ]
)

cat(
  "Total blocks: ",
  N_BLOCKS,
  "\n"
)

cat(
  "Shared-control blocks (k >= 2): ",
  N_SHARED_BLOCKS,
  "\n"
)

cat(
  "Effect sizes in shared-control blocks: ",
  N_EFFECTS_SHARED,
  "\n\n"
)

# 10. Initialize the primary analysis

cat(
  "Step 4/12: Initializing the Nakagawa analysis...\n"
)

dat_main <- raw_dat

dat_main$yi <- NA_real_

dat_main$vi <- NA_real_

dat_main$Cohen_d <- NA_real_

dat_main$Hedges_J <- NA_real_

dat_main$pooled_SD <- NA_real_

dat_main$block_N <- NA_real_

dat_main$block_df <- NA_real_

dat_main$block_k <- NA_integer_

# Record the control values used for each block

dat_main$shared_control_n <- NA_real_

dat_main$shared_control_mean <- NA_real_

dat_main$shared_control_sd <- NA_real_

# Preserve indices for covariance-matrix construction

dat_main$V_index <- seq_len(
  RAW_ROWS
)

# 11. Initialize sparse covariance-matrix storage

V_i_list <- vector(
  "list",
  N_BLOCKS
)

V_j_list <- vector(
  "list",
  N_BLOCKS
)

V_x_list <- vector(
  "list",
  N_BLOCKS
)

block_qc_list <- vector(
  "list",
  N_BLOCKS
)

covariance_qc_list <- list()

group_counter <- 0

cov_counter <- 0

# 12. Calculate effect sizes and within-block covariance

cat(
  "Step 5/12: Calculating Hedges' g and the covariance matrix by block...\n"
)

for (block_name in names(groups)) {

  group_counter <- group_counter + 1

  # Current block indices

  idx <- groups[[block_name]]

  x <- dat_main[
    idx,
    ,
    drop = FALSE
  ]

  k <- NROW(
    x
  )

  # 12.1 Verify block identifiers

  if (
    length(
      unique(
        x$Study_id
      )
    ) != 1
  ) {

    stop(
      paste(
        "Multiple Study_id values found within block: ",
        block_name
      )
    )
  }

  if (
    length(
      unique(
        x$Trait
      )
    ) != 1
  ) {

    stop(
      paste(
        "Multiple Trait values found within block: ",
        block_name
      )
    )
  }

  if (
    length(
      unique(
        x$Common_id
      )
    ) != 1
  ) {

    stop(
      paste(
        "Multiple Common_id values found within block: ",
        block_name
      )
    )
  }

  # 12.2 Select the shared control

  nc <- as.numeric(
    x$n_control[1]
  )

  mean_C <- as.numeric(
    x$mean_control[1]
  )

  sd_C <- as.numeric(
    x$sd_control[1]
  )

  # 12.3 Treatment sample sizes

  nt <- as.numeric(
    x$n_treatment
  )

  # 12.4 Total sample size

  N_block <- (

    nc +

    sum(
      nt
    )
  )

  # 12.5 Degrees of freedom

  df_block <- (

    N_block -

    (
      k +
      1
    )
  )

  if (
    !is.finite(
      df_block
    ) ||

    df_block <= 1
  ) {

    stop(
      paste(
        "Insufficient degrees of freedom in block: ",
        block_name
      )
    )
  }

  # 12.6 Pooled variance

  sp2 <- (

    (

      (nc - 1) *

      sd_C^2
    ) +

    sum(

      (nt - 1) *

      x$sd_treatment^2
    )

  ) /

  df_block

  if (
    !is.finite(
      sp2
    ) ||

    sp2 <= 0
  ) {

    stop(
      paste(
        "Invalid pooled variance in block: ",
        block_name
      )
    )
  }

  sp <- sqrt(
    sp2
  )

  # 12.7 Cohen's d

  cohen_d <- (

    x$mean_treatment -

    mean_C

  ) /

  sp

  # 12.8 Exact Hedges' J correction

  J_block <- hedges_J_exact(
    df_block
  )

  if (
    length(
      J_block
    ) != 1 ||

    !is.finite(
      J_block
    )
  ) {

    stop(
      paste(
        "Failed to calculate Hedges' J in block: ",
        block_name
      )
    )
  }

  # 12.9 Hedges' g

  hedges_g <- (

    J_block *

    cohen_d
  )

  # 12.10 Nakagawa shared-control variance-covariance matrix

  B <- J_block^2 * (

    matrix(

      1 / nc,

      nrow = k,

      ncol = k

    ) +

    outer(

      cohen_d,

      cohen_d

    ) /

    (
      2 *
      N_block
    )
  )

  # 12.11 Add treatment sampling variance to the diagonal

  diag(
    B
  ) <- (

    diag(
      B
    ) +

    J_block^2 /
    nt
  )

  # 12.12 Validate the block covariance matrix

  if (
    any(
      !is.finite(
        B
      )
    )
  ) {

    stop(
      paste(
        "Non-finite values in the block covariance matrix: ",
        block_name
      )
    )
  }

  if (
    any(
      diag(
        B
      ) <= 0
    )
  ) {

    stop(
      paste(
        "Non-positive sampling variances in the block covariance matrix: ",
        block_name
      )
    )
  }

  # 12.13 Store effect-size results

  dat_main$yi[
    idx
  ] <- hedges_g

  dat_main$vi[
    idx
  ] <- diag(
    B
  )

  dat_main$Cohen_d[
    idx
  ] <- cohen_d

  dat_main$Hedges_J[
    idx
  ] <- J_block

  dat_main$pooled_SD[
    idx
  ] <- sp

  dat_main$block_N[
    idx
  ] <- N_block

  dat_main$block_df[
    idx
  ] <- df_block

  dat_main$block_k[
    idx
  ] <- k

  dat_main$shared_control_n[
    idx
  ] <- nc

  dat_main$shared_control_mean[
    idx
  ] <- mean_C

  dat_main$shared_control_sd[
    idx
  ] <- sd_C

  # 12.14 Store sparse covariance entries

  B_row <- row(
    B
  )

  B_col <- col(
    B
  )

  V_i_list[[group_counter]] <- idx[
    as.vector(
      B_row
    )
  ]

  V_j_list[[group_counter]] <- idx[
    as.vector(
      B_col
    )
  ]

  V_x_list[[group_counter]] <- as.vector(
    B
  )

  # 12.15 Block-level quality control

  block_qc_list[[group_counter]] <- data.frame(

    block =
      as.character(
        block_name
      ),

    Study_id =
      as.character(
        x$Study_id[1]
      ),

    Trait =
      as.character(
        x$Trait[1]
      ),

    Common_id =
      as.character(
        x$Common_id[1]
      ),

    k =
      k,

    Shared_control_original_row =
      x$original_row[1],

    Shared_control_n =
      nc,

    Shared_control_mean =
      mean_C,

    Shared_control_sd =
      sd_C,

    Sum_n_treatment =
      sum(
        nt
      ),

    N_total =
      N_block,

    df =
      df_block,

    pooled_SD =
      sp,

    Hedges_J =
      J_block,

    stringsAsFactors =
      FALSE
  )

  # 12.16 Pairwise covariance diagnostics

  if (k >= 2) {

    pair_index <- combn(

      seq_len(
        k
      ),

      2
    )

    for (
      p in seq_len(
        ncol(
          pair_index
        )
      )
    ) {

      local_i <- pair_index[
        1,
        p
      ]

      local_j <- pair_index[
        2,
        p
      ]

      global_i <- idx[
        local_i
      ]

      global_j <- idx[
        local_j
      ]

      covariance_value <- B[
        local_i,
        local_j
      ]

      rho_value <- (

        covariance_value /

        sqrt(

          B[
            local_i,
            local_i
          ] *

          B[
            local_j,
            local_j
          ]
        )
      )

      cov_counter <- cov_counter + 1

      covariance_qc_list[[cov_counter]] <- data.frame(

        Index_1 =
          global_i,

        Index_2 =
          global_j,

        Original_row_1 =
          x$original_row[
            local_i
          ],

        Original_row_2 =
          x$original_row[
            local_j
          ],

        Study_id =
          as.character(
            x$Study_id[1]
          ),

        Trait =
          as.character(
            x$Trait[1]
          ),

        Common_id =
          as.character(
            x$Common_id[1]
          ),

        k =
          k,

        N_total =
          N_block,

        df =
          df_block,

        n_control =
          nc,

        n_treatment_1 =
          nt[
            local_i
          ],

        n_treatment_2 =
          nt[
            local_j
          ],

        Cohen_d_1 =
          cohen_d[
            local_i
          ],

        Cohen_d_2 =
          cohen_d[
            local_j
          ],

        Hedges_g_1 =
          hedges_g[
            local_i
          ],

        Hedges_g_2 =
          hedges_g[
            local_j
          ],

        Covariance =
          covariance_value,

        Sampling_correlation =
          rho_value,

        stringsAsFactors =
          FALSE
      )
    }
  }
}

# 13. Assemble the full sparse covariance matrix

cat(
  "Step 6/12: Assembling the full covariance matrix...\n"
)

V_i <- unlist(

  V_i_list,

  use.names = FALSE
)

V_j <- unlist(

  V_j_list,

  use.names = FALSE
)

V_x <- unlist(

  V_x_list,

  use.names = FALSE
)

V_main <- Matrix::sparseMatrix(

  i =
    V_i,

  j =
    V_j,

  x =
    V_x,

  dims = c(

    RAW_ROWS,

    RAW_ROWS
  ),

  giveCsparse =
    TRUE
)

# 14. Validate effect-size estimates

if (
  any(
    !is.finite(
      dat_main$yi
    )
  )
) {

  stop(
    "Some effect sizes (yi) were not calculated."
  )
}

if (
  any(
    !is.finite(
      dat_main$vi
    )
  )
) {

  stop(
    "Some sampling variances (vi) were not calculated."
  )
}

if (
  any(
    dat_main$vi <= 0
  )
) {

  stop(
    "Non-positive sampling variances (vi) detected."
  )
}

# 15. Assemble block-level diagnostics

block_qc_df <- do.call(

  rbind,

  block_qc_list
)

rownames(
  block_qc_df
) <- NULL

# 16. Assemble covariance diagnostics

if (
  length(
    covariance_qc_list
  ) > 0
) {

  covariance_qc_df <- do.call(

    rbind,

    covariance_qc_list
  )

  rownames(
    covariance_qc_df
  ) <- NULL

} else {

  covariance_qc_df <- data.frame(

    Index_1 =
      integer(0),

    Index_2 =
      integer(0),

    Original_row_1 =
      integer(0),

    Original_row_2 =
      integer(0),

    Study_id =
      character(0),

    Trait =
      character(0),

    Common_id =
      character(0),

    k =
      integer(0),

    N_total =
      numeric(0),

    df =
      numeric(0),

    n_control =
      numeric(0),

    n_treatment_1 =
      numeric(0),

    n_treatment_2 =
      numeric(0),

    Cohen_d_1 =
      numeric(0),

    Cohen_d_2 =
      numeric(0),

    Hedges_g_1 =
      numeric(0),

    Hedges_g_2 =
      numeric(0),

    Covariance =
      numeric(0),

    Sampling_correlation =
      numeric(0),

    stringsAsFactors =
      FALSE
  )
}

# 17. Validate the Nakagawa covariance matrix

cat("\n")

cat(
  strrep(
    "=",
    90
  ),
  "\n"
)

cat(
  "NAKAGAWA V QC\n"
)

cat(
  strrep(
    "=",
    90
  ),
  "\n"
)

MAX_DIAG_DIFF <- max(

  abs(

    as.numeric(
      diag(
        V_main
      )
    ) -

    dat_main$vi
  ),

  na.rm = TRUE
)

V_MAIN_SYMMETRIC <- isSymmetric(
  V_main
)

N_NONZERO <- Matrix::nnzero(
  V_main
)

N_OFFDIAGONAL <- (

  N_NONZERO -

  RAW_ROWS
)

MAX_ABS_RHO <- if (
  NROW(
    covariance_qc_df
  ) > 0
) {

  max(

    abs(
      covariance_qc_df$Sampling_correlation
    ),

    na.rm = TRUE
  )

} else {

  NA_real_
}

N_INVALID_RHO <- if (
  NROW(
    covariance_qc_df
  ) > 0
) {

  sum(

    abs(
      covariance_qc_df$Sampling_correlation
    ) >= 1,

    na.rm = TRUE
  )

} else {

  0
}

cat(
  "Covariance matrix dimensions: ",
  NROW(
    V_main
  ),
  " x ",
  NCOL(
    V_main
  ),
  "\n"
)

cat(
  "Maximum absolute difference between diag(V) and vi: ",
  format(
    MAX_DIAG_DIFF,
    scientific = TRUE
  ),
  "\n"
)

cat(
  "Is V symmetric: ",
  V_MAIN_SYMMETRIC,
  "\n"
)

cat(
  "Non-zero off-diagonal elements: ",
  N_OFFDIAGONAL,
  "\n"
)

cat(
  "Covariance pairs:",
  NROW(
    covariance_qc_df
  ),
  "\n"
)

cat(
  "Maximum |rho|: ",
  MAX_ABS_RHO,
  "\n"
)

cat(
  "Number of pairs with |rho| >= 1: ",
  N_INVALID_RHO,
  "\n\n"
)

if (
  MAX_DIAG_DIFF > 1e-10
) {

  stop(
    "The diagonal of V does not match vi."
  )
}

if (
  !V_MAIN_SYMMETRIC
) {

  stop(
    "The covariance matrix V is not symmetric."
  )
}

# 18. Check block-level positive definiteness

cat(
  "Step 7/12: Checking block-level positive definiteness...\n"
)

pd_list <- vector(

  "list",

  N_BLOCKS
)

pd_counter <- 0

for (block_name in names(groups)) {

  pd_counter <- pd_counter + 1

  idx <- groups[[block_name]]

  V_block <- as.matrix(

    V_main[
      idx,
      idx,
      drop = FALSE
    ]
  )

  eig_values <- eigen(

    V_block,

    symmetric = TRUE,

    only.values = TRUE

  )$values

  pd_list[[pd_counter]] <- data.frame(

    block =
      as.character(
        block_name
      ),

    Study_id =
      as.character(
        dat_main$Study_id[
          idx[1]
        ]
      ),

    Trait =
      as.character(
        dat_main$Trait[
          idx[1]
        ]
      ),

    Common_id =
      as.character(
        dat_main$Common_id[
          idx[1]
        ]
      ),

    n_effects =
      length(
        idx
      ),

    Min_eigenvalue =
      min(
        eig_values
      ),

    Positive_definite =
      all(
        eig_values > 1e-10
      ),

    stringsAsFactors =
      FALSE
  )
}

block_pd_df <- do.call(

  rbind,

  pd_list
)

rownames(
  block_pd_df
) <- NULL

N_NON_PD_BLOCKS <- sum(

  !block_pd_df$Positive_definite,

  na.rm = TRUE
)

cat(
  "Non-positive-definite blocks: ",
  N_NON_PD_BLOCKS,
  "\n\n"
)

# 19. Stop if covariance-matrix QC fails

precheck_file <- paste0(

  prefix,

  "_Nakagawa_V_precheck.xlsx"
)

writexl::write_xlsx(

  list(

    Block_summary =
      block_qc_df,

    Covariance_pairs =
      covariance_qc_df,

    Block_PD =
      block_pd_df

  ),

  precheck_file
)

if (
  N_INVALID_RHO > 0 ||
  N_NON_PD_BLOCKS > 0
) {

  stop(

    paste0(

      "The Nakagawa covariance matrix failed quality control.\n",

      "|rho|>=1 = ",
      N_INVALID_RHO,
      "\n",

      "Non-positive-definite blocks = ",
      N_NON_PD_BLOCKS,
      "\n\n",

      "See the diagnostic file:\n",

      precheck_file
    )
  )
}

# 20. Fit the primary Nakagawa model

cat(
  "Step 8/12: Fitting the primary Nakagawa model...\n"
)

model_main <- metafor::rma.mv(

  yi =
    yi,

  V =
    V_main,

  random =
    ~ 1 | Study_id / Row_id,

  data =
    dat_main,

  method =
    "REML",

  sparse =
    TRUE
)

print(
  summary(
    model_main
  )
)

# 21. Extract model estimates

cat(
  "\nStep 9/12: Extracting model results...\n"
)

extract_model_result <- function(
    model,
    analysis_name
) {

  data.frame(

    Analysis =
      analysis_name,

    k =
      as.numeric(
        model$k
      ),

    Estimate =
      as.numeric(
        model$b[1]
      ),

    SE =
      as.numeric(
        model$se[1]
      ),

    CI_lower =
      as.numeric(
        model$ci.lb[1]
      ),

    CI_upper =
      as.numeric(
        model$ci.ub[1]
      ),

    p_value =
      as.numeric(
        model$pval[1]
      ),

    QE =
      if (
        length(
          model$QE
        ) == 1
      ) {

        as.numeric(
          model$QE
        )

      } else {

        NA_real_
      },

    QE_p =
      if (
        length(
          model$QEp
        ) == 1
      ) {

        as.numeric(
          model$QEp
        )

      } else {

        NA_real_
      },

    stringsAsFactors =
      FALSE
  )
}

main_results <- extract_model_result(
  model_main,
  "Main_Nakagawa_shared_control_VCV"
)

# 22. Compile the quality-control summary

cat(
  "Step 10/12: Compiling quality-control summaries...\n"
)

make_qc_row <- function(
    item,
    value
) {

  if (
    length(
      value
    ) == 0
  ) {

    value <- NA_character_
  }

  if (
    length(
      value
    ) > 1
  ) {

    value <- paste(

      value,

      collapse = ";"
    )
  }

  data.frame(

    Item =
      as.character(
        item
      ),

    Value =
      as.character(
        value
      ),

    stringsAsFactors =
      FALSE
  )
}

QC_summary <- do.call(

  rbind,

  list(

    make_qc_row(
      "Raw_rows",
      RAW_ROWS
    ),

    make_qc_row(
      "Study_number",
      length(
        unique(
          dat_main$Study_id
        )
      )
    ),

    make_qc_row(
      "Trait_number",
      length(
        unique(
          dat_main$Trait
        )
      )
    ),

    make_qc_row(
      "Total_blocks",
      N_BLOCKS
    ),

    make_qc_row(
      "Shared_control_blocks_k_ge_2",
      N_SHARED_BLOCKS
    ),

    make_qc_row(
      "Effects_in_shared_blocks",
      N_EFFECTS_SHARED
    ),

    make_qc_row(
      "Covariance_pairs",
      NROW(
        covariance_qc_df
      )
    ),

    make_qc_row(
      "Nonzero_offdiagonal_elements",
      N_OFFDIAGONAL
    ),

    make_qc_row(
      "Max_abs_sampling_correlation",
      MAX_ABS_RHO
    ),

    make_qc_row(
      "Invalid_rho_pairs_abs_ge_1",
      N_INVALID_RHO
    ),

    make_qc_row(
      "Non_PD_blocks",
      N_NON_PD_BLOCKS
    ),

    make_qc_row(
      "Max_diag_V_minus_vi",
      MAX_DIAG_DIFF
    ),

    make_qc_row(
      "V_main_symmetric",
      V_MAIN_SYMMETRIC
    )
  )
)

print(
  QC_summary
)

# 23. Define output files

qc_excel <- paste0(prefix, "_Nakagawa_main_QC.xlsx")
main_data_rds <- paste0(prefix, "_d2_Nakagawa_main.rds")
main_V_rds <- paste0(prefix, "_V_Nakagawa_main.rds")
main_model_rds <- paste0(prefix, "_model_Nakagawa_main.rds")
summary_txt <- paste0(prefix, "_Nakagawa_main_summary.txt")

# 24. Export primary analysis results

cat("Step 11/12: Exporting Excel and RDS files...\n")

export_dat <- as.data.frame(dat_main)
export_dat$V_diag <- as.numeric(diag(V_main))

writexl::write_xlsx(
  list(
    Effect_sizes = export_dat,
    Main_model = main_results,
    QC_summary = QC_summary,
    Block_summary = block_qc_df,
    Covariance_pairs = covariance_qc_df,
    Block_PD = block_pd_df
  ),
  qc_excel
)

saveRDS(dat_main, main_data_rds)
saveRDS(V_main, main_V_rds)
saveRDS(model_main, main_model_rds)

# 25. Export the primary model summary

cat("Step 12/12: Exporting the primary model summary...\n")

sink(summary_txt)
cat("Nakagawa shared-control SMD VCV: primary analysis\n\n")
cat("Run time: ", as.character(Sys.time()), "\n\n", sep = "")
cat("Shared-control definition: Study_id + Trait + Common_id\n\n")
print(QC_summary)
cat("\nPRIMARY NAKAGAWA MODEL\n\n")
print(summary(model_main))
cat("\nPRIMARY MODEL ESTIMATES\n\n")
print(main_results)
sink()

# 26. Completion

cat("\nPrimary analysis completed.\n\n")
print(QC_summary)
print(main_results)
cat("\nPrimary analysis files:\n")
cat(main_data_rds, "\n", main_V_rds, "\n", main_model_rds, "\n", sep = "")
cat("\nQuality-control workbook:\n", qc_excel, "\n", sep = "")
cat("\nModel summary:\n", summary_txt, "\n", sep = "")
