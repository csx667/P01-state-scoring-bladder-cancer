options(
  stringsAsFactors = FALSE,
  warn = 1
)

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
})

set.seed(20260801)

BASE <- "/root/autodl-tmp/scrna_project"

EXPR_FILE <- file.path(
  BASE,
  "data",
  "BLCA_expression_matrix.csv"
)

CLINICAL_FILE <- file.path(
  BASE,
  "data",
  "TCGA_BLCA_survival_01_for_analysis.tsv"
)

REFERENCE_FILE <- file.path(
  BASE,
  "result",
  "37_BayesPrism_rebuild_v5",
  "37B_reference_3state",
  "37B_reference_counts_main_LAMC2_excluded.rds"
)

LABEL_FILE <- file.path(
  BASE,
  "result",
  "37_BayesPrism_rebuild_v5",
  "37B_reference_3state",
  "37B_reference_cell_labels.rds"
)

OUT <- file.path(
  BASE,
  "result",
  "37_BayesPrism_rebuild_v5",
  "37K1_P01_LAMC2_excluded_ITGB4_eligible"
)

dir.create(
  OUT,
  recursive = TRUE,
  showWarnings = FALSE
)

STATUS_FILE <- file.path(
  OUT,
  "37I_RUN_STATUS.txt"
)

REVIEW_FILE <- file.path(
  OUT,
  "37I_REVIEW_CN.txt"
)

LOG_FILE <- file.path(
  OUT,
  "37I_run_log.txt"
)

BULK_QC_FILE <- file.path(
  OUT,
  "37I_bulk_matrix_QC.tsv"
)

REFERENCE_QC_FILE <- file.path(
  OUT,
  "37I_reference_structure_QC.tsv"
)

STATE_MARKER_FILE <- file.path(
  OUT,
  "37I_LAMC2_high_vs_non_high_marker_table.tsv"
)

EPI_MARKER_FILE <- file.path(
  OUT,
  "37I_epithelial_vs_non_epithelial_marker_table.tsv"
)

MARKER_QC_FILE <- file.path(
  OUT,
  "37I_marker_selection_QC.tsv"
)

SCORE_FILE <- file.path(
  OUT,
  "37I_TCGA_processed_bulk_rank_scores_408.tsv"
)

CLINICAL_SCORE_FILE <- file.path(
  OUT,
  "37I_TCGA_rank_scores_clinical_matched.tsv"
)

CORRELATION_FILE <- file.path(
  OUT,
  "37I_score_feature_correlations.tsv"
)

SCORE_RDS <- file.path(
  OUT,
  "37I_TCGA_processed_bulk_rank_scores_408.rds"
)

EXPECTED_TUMORS <- 408L
MAX_MARKERS_PER_DIRECTION <- 120L
MIN_MARKERS_PER_DIRECTION <- 25L
MAX_BULK_MISSING_FRACTION <- 0.20

log_con <- file(
  LOG_FILE,
  open = "wt"
)

on.exit(
  try(
    close(log_con),
    silent = TRUE
  ),
  add = TRUE
)

log_msg <- function(...) {
  text <- paste0(
    ...,
    collapse = ""
  )

  cat(
    text,
    "\n"
  )

  cat(
    text,
    "\n",
    file = log_con
  )

  flush(
    log_con
  )
}

fail <- function(...) {
  error_text <- paste0(
    ...,
    collapse = ""
  )

  log_msg(
    "[FATAL] ",
    error_text
  )

  stop(
    error_text,
    call. = FALSE
  )
}

write_tsv <- function(
    x,
    path
) {
  write.table(
    x,
    path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = "NA"
  )
}

normalize_tcga_id <- function(x) {
  x <- toupper(
    trimws(
      as.character(x)
    )
  )

  gsub(
    "\\.",
    "-",
    x
  )
}

patient_id_from_barcode <- function(x) {
  x <- normalize_tcga_id(x)

  result <- rep(
    NA_character_,
    length(x)
  )

  valid <- !is.na(x) &
    nchar(x) >= 12L

  result[valid] <- substr(
    x[valid],
    1L,
    12L
  )

  result
}

sample_type_from_barcode <- function(x) {
  x <- normalize_tcga_id(x)

  result <- rep(
    NA_character_,
    length(x)
  )

  valid <- !is.na(x) &
    nchar(x) >= 15L

  result[valid] <- substr(
    x[valid],
    14L,
    15L
  )

  result
}

safe_z <- function(x) {
  result <- rep(
    NA_real_,
    length(x)
  )

  ok <- is.finite(x)

  if (sum(ok) < 3L) {
    return(result)
  }

  standard_deviation <- sd(
    x[ok]
  )

  if (
    !is.finite(standard_deviation) ||
    standard_deviation <= 0
  ) {
    return(result)
  }

  result[ok] <- (
    x[ok] -
      mean(
        x[ok]
      )
  ) / standard_deviation

  result
}

collapse_duplicate_genes <- function(
    matrix_value,
    gene_names
) {
  gene_names <- toupper(
    trimws(
      as.character(gene_names)
    )
  )

  valid <- !is.na(gene_names) &
    nzchar(gene_names) &
    gene_names != "NA"

  matrix_value <- matrix_value[
    valid,
    ,
    drop = FALSE
  ]

  gene_names <- gene_names[
    valid
  ]

  if (!anyDuplicated(gene_names)) {
    rownames(matrix_value) <- gene_names

    return(matrix_value)
  }

  zero_filled <- matrix_value

  zero_filled[
    !is.finite(
      zero_filled
    )
  ] <- 0

  observed <- matrix(
    as.numeric(
      is.finite(
        matrix_value
      )
    ),
    nrow = nrow(matrix_value),
    ncol = ncol(matrix_value)
  )

  sum_matrix <- rowsum(
    zero_filled,
    group = gene_names,
    reorder = FALSE,
    na.rm = TRUE
  )

  count_matrix <- rowsum(
    observed,
    group = gene_names,
    reorder = FALSE,
    na.rm = TRUE
  )

  result <- sum_matrix /
    count_matrix

  result[
    count_matrix == 0
  ] <- NA_real_

  result
}

read_processed_bulk <- function() {
  log_msg(
    "Reading the original processed expression matrix: ",
    EXPR_FILE
  )

  bulk_table <- data.table::fread(
    EXPR_FILE,
    data.table = FALSE,
    check.names = FALSE,
    na.strings = c(
      "NA",
      "NaN",
      ""
    ),
    showProgress = TRUE
  )

  if (
    nrow(bulk_table) < 10000L ||
    ncol(bulk_table) < 400L
  ) {
    fail(
      "Unexpected BLCA expression matrix dimensions: ",
      nrow(bulk_table),
      " × ",
      ncol(bulk_table)
    )
  }

  gene_column <- colnames(
    bulk_table
  )[1L]

  gene_names <- as.character(
    bulk_table[[1L]]
  )

  value_table <- bulk_table[
    ,
    -1L,
    drop = FALSE
  ]

  non_numeric <- colnames(value_table)[
    !vapply(
      value_table,
      is.numeric,
      logical(1)
    )
  ]

  if (length(non_numeric) > 0L) {
    log_msg(
      "Converting non-numeric expression columns to numeric: ",
      paste(
        non_numeric,
        collapse = " | "
      )
    )

    value_table[] <- lapply(
      value_table,
      function(x) {
        suppressWarnings(
          as.numeric(
            as.character(x)
          )
        )
      }
    )
  }

  expression <- as.matrix(
    value_table
  )

  storage.mode(expression) <- "double"

  sample_ids_original <- colnames(
    expression
  )

  sample_ids_normalized <- normalize_tcga_id(
    sample_ids_original
  )

  sample_type <- sample_type_from_barcode(
    sample_ids_normalized
  )

  tumor_columns <- which(
    sample_type == "01"
  )

  normal_columns <- which(
    sample_type == "11"
  )

  sample_type_qc <- data.frame(
    sample_type_code = sort(
      unique(
        sample_type
      )
    ),
    n_samples = as.integer(
      table(
        sample_type
      )[
        sort(
          unique(
            sample_type
          )
        )
      ]
    ),
    stringsAsFactors = FALSE
  )

  write_tsv(
    sample_type_qc,
    file.path(
      OUT,
      "37I_bulk_sample_type_QC.tsv"
    )
  )

  if (length(tumor_columns) != EXPECTED_TUMORS) {
    fail(
      "Expected the number of primary tumor (01) samples in the original matrix to be ",
      EXPECTED_TUMORS,
      ", but found ",
      length(tumor_columns),
      "."
    )
  }

  expression <- expression[
    ,
    tumor_columns,
    drop = FALSE
  ]

  colnames(expression) <-
    sample_ids_normalized[
      tumor_columns
    ]

  expression <- collapse_duplicate_genes(
    matrix_value = expression,
    gene_names = gene_names
  )

  if (anyDuplicated(colnames(expression))) {
    fail(
      "Duplicate sample names in the primary tumor matrix."
    )
  }

  patient_ids <- patient_id_from_barcode(
    colnames(expression)
  )

  if (anyDuplicated(patient_ids)) {
    fail(
      "Duplicate patients in the primary tumor matrix."
    )
  }

  finite_values <- expression[
    is.finite(
      expression
    )
  ]

  if (length(finite_values) == 0L) {
    fail(
      "The expression matrix contains no finite values."
    )
  }

  bulk_qc <- data.frame(
    source_file =
      EXPR_FILE,
    first_column_name =
      gene_column,
    original_feature_rows =
      nrow(bulk_table),
    original_sample_columns =
      ncol(bulk_table) - 1L,
    primary_tumor_samples =
      length(tumor_columns),
    normal_samples =
      length(normal_columns),
    final_unique_gene_symbols =
      nrow(expression),
    finite_values =
      length(finite_values),
    missing_values =
      sum(
        !is.finite(
          expression
        )
      ),
    missing_fraction =
      mean(
        !is.finite(
          expression
        )
      ),
    minimum_value =
      min(finite_values),
    first_quartile =
      as.numeric(
        quantile(
          finite_values,
          0.25,
          names = FALSE
        )
      ),
    median_value =
      median(finite_values),
    third_quartile =
      as.numeric(
        quantile(
          finite_values,
          0.75,
          names = FALSE
        )
      ),
    maximum_value =
      max(finite_values),
    negative_fraction =
      mean(
        finite_values < 0
      ),
    exact_zero_fraction =
      mean(
        finite_values == 0
      ),
    integer_fraction =
      mean(
        abs(
          finite_values -
            round(
              finite_values
            )
        ) <
          1e-8
      ),
    stringsAsFactors = FALSE
  )

  write_tsv(
    bulk_qc,
    BULK_QC_FILE
  )

  list(
    expression = expression,
    qc = bulk_qc
  )
}

orient_reference <- function(
    reference,
    labels
) {
  required_label_columns <- c(
    "cell_id",
    "Patient",
    "cell_type",
    "cell_state"
  )

  missing_label_columns <- setdiff(
    required_label_columns,
    colnames(labels)
  )

  if (length(missing_label_columns) > 0L) {
    fail(
      "Missing columns in 37B labels: ",
      paste(
        missing_label_columns,
        collapse = ", "
      )
    )
  }

  labels$cell_id <- as.character(
    labels$cell_id
  )

  column_overlap <- if (
    is.null(
      colnames(reference)
    )
  ) {
    0L
  } else {
    sum(
      labels$cell_id %in%
        colnames(reference)
    )
  }

  row_overlap <- if (
    is.null(
      rownames(reference)
    )
  ) {
    0L
  } else {
    sum(
      labels$cell_id %in%
        rownames(reference)
    )
  }

  if (
    column_overlap == 0L &&
    row_overlap == 0L
  ) {
    fail(
      "Neither row names nor column names of the 37B reference match cell_id."
    )
  }

  if (column_overlap >= row_overlap) {
    position <- match(
      labels$cell_id,
      colnames(reference)
    )

    valid <- !is.na(position)

    labels <- labels[
      valid,
      ,
      drop = FALSE
    ]

    reference <- reference[
      ,
      position[valid],
      drop = FALSE
    ]

    orientation <-
      "genes_by_cells"
  } else {
    position <- match(
      labels$cell_id,
      rownames(reference)
    )

    valid <- !is.na(position)

    labels <- labels[
      valid,
      ,
      drop = FALSE
    ]

    reference <- Matrix::t(
      reference[
        position[valid],
        ,
        drop = FALSE
      ]
    )

    orientation <-
      "cells_by_genes_transposed_to_genes_by_cells"
  }

  colnames(reference) <-
    labels$cell_id

  if (is.null(rownames(reference))) {
    fail(
      "The 37B reference matrix has no gene row names."
    )
  }

  rownames(reference) <- toupper(
    trimws(
      rownames(reference)
    )
  )

  if (anyDuplicated(rownames(reference))) {
    fail(
      "Duplicate gene names in the 37B reference after conversion to uppercase."
    )
  }

  reference_qc <- data.frame(
    orientation_detected =
      orientation,
    reference_genes =
      nrow(reference),
    reference_cells =
      ncol(reference),
    label_rows =
      nrow(labels),
    unique_patients =
      length(
        unique(
          labels$Patient
        )
      ),
    unique_cell_types =
      length(
        unique(
          labels$cell_type
        )
      ),
    unique_cell_states =
      length(
        unique(
          labels$cell_state
        )
      ),
    stringsAsFactors = FALSE
  )

  write_tsv(
    reference_qc,
    REFERENCE_QC_FILE
  )

  list(
    reference = reference,
    labels = labels,
    qc = reference_qc
  )
}

aggregate_patient_group <- function(
    reference,
    patient,
    group_label
) {
  group_key <- paste(
    patient,
    group_label,
    sep = "|||GROUP|||"
  )

  group_factor <- factor(
    group_key
  )

  design <- Matrix::sparseMatrix(
    i = seq_along(group_factor),
    j = as.integer(group_factor),
    x = 1,
    dims = c(
      length(group_factor),
      nlevels(group_factor)
    ),
    dimnames = list(
      NULL,
      levels(group_factor)
    )
  )

  group_sum <- reference %*%
    design

  library_size <- Matrix::colSums(
    group_sum
  )

  if (
    any(
      !is.finite(library_size)
    ) ||
    any(
      library_size <= 0
    )
  ) {
    fail(
      "Invalid library size in reference pseudobulk."
    )
  }

  group_log_cpm <- log1p(
    sweep(
      as.matrix(group_sum),
      2L,
      library_size,
      "/"
    ) *
      10000
  )

  group_parts <- strsplit(
    colnames(group_log_cpm),
    "|||GROUP|||",
    fixed = TRUE
  )

  group_metadata <- data.frame(
    column_name =
      colnames(group_log_cpm),
    Patient = vapply(
      group_parts,
      `[`,
      character(1),
      1L
    ),
    group_label = vapply(
      group_parts,
      `[`,
      character(1),
      2L
    ),
    stringsAsFactors = FALSE
  )

  list(
    expression = group_log_cpm,
    metadata = group_metadata
  )
}

derive_paired_markers <- function(
    reference,
    labels,
    positive_label,
    negative_label,
    contrast_name,
    bulk_expression,
    excluded_genes = c(
      "LAMC2",
      "ITGB4"
    )
) {
  keep <- labels$contrast_group %in%
    c(
      positive_label,
      negative_label
    )

  contrast_reference <- reference[
    ,
    keep,
    drop = FALSE
  ]

  contrast_labels <- labels[
    keep,
    ,
    drop = FALSE
  ]

  aggregated <- aggregate_patient_group(
    reference =
      contrast_reference,
    patient =
      contrast_labels$Patient,
    group_label =
      contrast_labels$contrast_group
  )

  positive_patients <- aggregated$metadata$
    Patient[
      aggregated$metadata$
        group_label ==
        positive_label
    ]

  negative_patients <- aggregated$metadata$
    Patient[
      aggregated$metadata$
        group_label ==
        negative_label
    ]

  common_patients <- sort(
    intersect(
      positive_patients,
      negative_patients
    )
  )

  if (length(common_patients) < 3L) {
    fail(
      contrast_name,
      "Only ",
      length(common_patients),
      " paired patients; stable markers cannot be established."
    )
  }

  positive_columns <- match(
    paste(
      common_patients,
      positive_label,
      sep = "|||GROUP|||"
    ),
    colnames(
      aggregated$expression
    )
  )

  negative_columns <- match(
    paste(
      common_patients,
      negative_label,
      sep = "|||GROUP|||"
    ),
    colnames(
      aggregated$expression
    )
  )

  positive_matrix <- aggregated$expression[
    ,
    positive_columns,
    drop = FALSE
  ]

  negative_matrix <- aggregated$expression[
    ,
    negative_columns,
    drop = FALSE
  ]

  difference_matrix <- positive_matrix -
    negative_matrix

  median_difference <- apply(
    difference_matrix,
    1L,
    median,
    na.rm = TRUE
  )

  mean_difference <- rowMeans(
    difference_matrix,
    na.rm = TRUE
  )

  positive_consistency <- rowMeans(
    difference_matrix > 0,
    na.rm = TRUE
  )

  negative_consistency <- rowMeans(
    difference_matrix < 0,
    na.rm = TRUE
  )

  bulk_missing_fraction <- rowMeans(
    !is.finite(
      bulk_expression
    )
  )

  bulk_standard_deviation <- apply(
    bulk_expression,
    1L,
    sd,
    na.rm = TRUE
  )

  gene_names <- rownames(
    reference
  )

  marker_table <- data.frame(
    contrast =
      contrast_name,
    gene =
      gene_names,
    paired_patients =
      length(
        common_patients
      ),
    median_difference_positive_minus_negative =
      median_difference,
    mean_difference_positive_minus_negative =
      mean_difference,
    positive_consistency =
      positive_consistency,
    negative_consistency =
      negative_consistency,
    bulk_present =
      gene_names %in%
        rownames(
          bulk_expression
        ),
    bulk_missing_fraction =
      bulk_missing_fraction[
        gene_names
      ],
    bulk_standard_deviation =
      bulk_standard_deviation[
        gene_names
      ],
    excluded_primary_gene =
      gene_names %in%
        excluded_genes,
    stringsAsFactors = FALSE
  )

  marker_table$bulk_missing_fraction[
    !is.finite(
      marker_table$
        bulk_missing_fraction
    )
  ] <- 1

  marker_table$bulk_standard_deviation[
    !is.finite(
      marker_table$
        bulk_standard_deviation
    )
  ] <- 0

  available <- marker_table$
    bulk_present &
    marker_table$
      bulk_missing_fraction <=
        MAX_BULK_MISSING_FRACTION &
    marker_table$
      bulk_standard_deviation >
        1e-8 &
    !marker_table$
      excluded_primary_gene

  positive_candidates <- which(
    available &
      marker_table$
        median_difference_positive_minus_negative >
        0 &
      marker_table$
        positive_consistency >=
        0.55
  )

  negative_candidates <- which(
    available &
      marker_table$
        median_difference_positive_minus_negative <
        0 &
      marker_table$
        negative_consistency >=
        0.55
  )

  if (
    length(positive_candidates) <
      MIN_MARKERS_PER_DIRECTION
  ) {
    positive_candidates <- which(
      available &
        marker_table$
          median_difference_positive_minus_negative >
          0
    )
  }

  if (
    length(negative_candidates) <
      MIN_MARKERS_PER_DIRECTION
  ) {
    negative_candidates <- which(
      available &
        marker_table$
          median_difference_positive_minus_negative <
          0
    )
  }

  positive_order <- positive_candidates[
    order(
      -marker_table$
        median_difference_positive_minus_negative[
          positive_candidates
        ],
      -marker_table$
        positive_consistency[
          positive_candidates
        ]
    )
  ]

  negative_order <- negative_candidates[
    order(
      marker_table$
        median_difference_positive_minus_negative[
          negative_candidates
        ],
      -marker_table$
        negative_consistency[
          negative_candidates
        ]
    )
  ]

  positive_selected <- head(
    positive_order,
    MAX_MARKERS_PER_DIRECTION
  )

  negative_selected <- head(
    negative_order,
    MAX_MARKERS_PER_DIRECTION
  )

  marker_table$selected_direction <-
    "not_selected"

  marker_table$selected_direction[
    positive_selected
  ] <- "positive_marker"

  marker_table$selected_direction[
    negative_selected
  ] <- "negative_marker"

  marker_table <- marker_table[
    order(
      factor(
        marker_table$selected_direction,
        levels = c(
          "positive_marker",
          "negative_marker",
          "not_selected"
        )
      ),
      -abs(
        marker_table$
          median_difference_positive_minus_negative
      )
    ),
    ,
    drop = FALSE
  ]

  positive_genes <- marker_table$gene[
    marker_table$
      selected_direction ==
      "positive_marker"
  ]

  negative_genes <- marker_table$gene[
    marker_table$
      selected_direction ==
      "negative_marker"
  ]

  if (
    length(positive_genes) <
      MIN_MARKERS_PER_DIRECTION ||
    length(negative_genes) <
      MIN_MARKERS_PER_DIRECTION
  ) {
    fail(
      contrast_name,
      "Insufficient final markers. Positive=",
      length(positive_genes),
      ", negative=",
      length(negative_genes),
      "."
    )
  }

  list(
    marker_table =
      marker_table,
    positive_genes =
      positive_genes,
    negative_genes =
      negative_genes,
    paired_patients =
      common_patients
  )
}

rank01 <- function(x) {
  result <- rep(
    NA_real_,
    length(x)
  )

  valid <- is.finite(x)

  n_valid <- sum(
    valid
  )

  if (n_valid == 0L) {
    return(result)
  }

  if (n_valid == 1L) {
    result[valid] <- 0.5

    return(result)
  }

  ranks <- rank(
    x[valid],
    ties.method = "average"
  )

  result[valid] <- (
    ranks - 1
  ) / (
    n_valid - 1
  )

  result
}

mean_rank_for_genes <- function(
    rank_vector,
    gene_names,
    matrix_genes
) {
  position <- match(
    gene_names,
    matrix_genes
  )

  position <- position[
    !is.na(position)
  ]

  values <- rank_vector[
    position
  ]

  values <- values[
    is.finite(values)
  ]

  if (length(values) == 0L) {
    return(
      c(
        score = NA_real_,
        n_genes = 0
      )
    )
  }

  c(
    score = mean(values),
    n_genes = length(values)
  )
}

score_bulk_samples <- function(
    bulk_expression,
    state_positive,
    state_negative,
    epi_positive,
    epi_negative
) {
  fixed_gene_sets <- list(
    Basal = c(
      "KRT5",
      "KRT6A",
      "KRT14",
      "TP63",
      "CD44",
      "EGFR"
    ),
    Luminal = c(
      "FOXA1",
      "GATA3",
      "PPARG",
      "KRT20",
      "UPK2",
      "FGFR3"
    ),
    Adhesion = c(
      "LAMA3",
      "LAMB3",
      "ITGA6"
    ),
    Invasion = c(
      "MMP14",
      "MMP2",
      "SERPINE1"
    ),
    Cell_cycle = c(
      "MKI67",
      "TOP2A",
      "BUB1",
      "TTK",
      "KIF23"
    )
  )

  matrix_genes <- rownames(
    bulk_expression
  )

  results <- vector(
    "list",
    ncol(
      bulk_expression
    )
  )

  for (
    sample_index in seq_len(
      ncol(
        bulk_expression
      )
    )
  ) {
    if (
      sample_index == 1L ||
      sample_index %% 50L == 0L ||
      sample_index ==
        ncol(bulk_expression)
    ) {
      log_msg(
        "Rank projection sample ",
        sample_index,
        "/",
        ncol(bulk_expression)
      )
    }

    one_expression <- bulk_expression[
      ,
      sample_index
    ]

    one_rank <- rank01(
      one_expression
    )

    state_positive_result <-
      mean_rank_for_genes(
        rank_vector = one_rank,
        gene_names =
          state_positive,
        matrix_genes =
          matrix_genes
      )

    state_negative_result <-
      mean_rank_for_genes(
        rank_vector = one_rank,
        gene_names =
          state_negative,
        matrix_genes =
          matrix_genes
      )

    epi_positive_result <-
      mean_rank_for_genes(
        rank_vector = one_rank,
        gene_names =
          epi_positive,
        matrix_genes =
          matrix_genes
      )

    epi_negative_result <-
      mean_rank_for_genes(
        rank_vector = one_rank,
        gene_names =
          epi_negative,
        matrix_genes =
          matrix_genes
      )

    fixed_scores <- lapply(
      fixed_gene_sets,
      function(gene_set) {
        mean_rank_for_genes(
          rank_vector =
            one_rank,
          gene_names =
            gene_set,
          matrix_genes =
            matrix_genes
        )
      }
    )

    get_expression <- function(gene) {
      position <- match(
        gene,
        matrix_genes
      )

      if (is.na(position)) {
        return(
          NA_real_
        )
      }

      as.numeric(
        one_expression[
          position
        ]
      )
    }

    get_rank <- function(gene) {
      position <- match(
        gene,
        matrix_genes
      )

      if (is.na(position)) {
        return(
          NA_real_
        )
      }

      as.numeric(
        one_rank[
          position
        ]
      )
    }

    results[[sample_index]] <- data.frame(
      sample_id =
        colnames(bulk_expression)[
          sample_index
        ],
      patient_id =
        patient_id_from_barcode(
          colnames(bulk_expression)[
            sample_index
          ]
        ),
      state_high_marker_mean_rank =
        state_positive_result[
          "score"
        ],
      state_non_high_marker_mean_rank =
        state_negative_result[
          "score"
        ],
      LAMC2_state_rank_score =
        state_positive_result[
          "score"
        ] -
        state_negative_result[
          "score"
        ],
      n_state_high_markers_used =
        state_positive_result[
          "n_genes"
        ],
      n_state_non_high_markers_used =
        state_negative_result[
          "n_genes"
        ],
      epithelial_marker_mean_rank =
        epi_positive_result[
          "score"
        ],
      non_epithelial_marker_mean_rank =
        epi_negative_result[
          "score"
        ],
      epithelial_rank_score =
        epi_positive_result[
          "score"
        ] -
        epi_negative_result[
          "score"
        ],
      n_epithelial_markers_used =
        epi_positive_result[
          "n_genes"
        ],
      n_non_epithelial_markers_used =
        epi_negative_result[
          "n_genes"
        ],
      LAMC2_expression_processed =
        get_expression(
          "LAMC2"
        ),
      LAMC2_within_sample_rank =
        get_rank(
          "LAMC2"
        ),
      ITGB4_expression_processed =
        get_expression(
          "ITGB4"
        ),
      ITGB4_within_sample_rank =
        get_rank(
          "ITGB4"
        ),
      Basal_rank_score =
        fixed_scores$Basal[
          "score"
        ],
      Luminal_rank_score =
        fixed_scores$Luminal[
          "score"
        ],
      Adhesion_rank_score =
        fixed_scores$Adhesion[
          "score"
        ],
      Invasion_rank_score =
        fixed_scores$Invasion[
          "score"
        ],
      Cell_cycle_rank_score =
        fixed_scores$Cell_cycle[
          "score"
        ],
      stringsAsFactors = FALSE
    )
  }

  score_table <- do.call(
    rbind,
    results
  )

  numeric_columns <- setdiff(
    colnames(score_table),
    c(
      "sample_id",
      "patient_id"
    )
  )

  score_table[
    numeric_columns
  ] <- lapply(
    score_table[
      numeric_columns
    ],
    as.numeric
  )

  score_table
}

make_adjusted_scores <- function(score_table) {
  score_table$LAMC2_state_rank_score_z <-
    safe_z(
      score_table$
        LAMC2_state_rank_score
    )

  score_table$epithelial_rank_score_z <-
    safe_z(
      score_table$
        epithelial_rank_score
    )

  complete <- is.finite(
    score_table$
      LAMC2_state_rank_score
  ) &
    is.finite(
      score_table$
        epithelial_rank_score
    )

  adjusted <- rep(
    NA_real_,
    nrow(score_table)
  )

  if (sum(complete) >= 20L) {
    fit <- lm(
      LAMC2_state_rank_score ~
        epithelial_rank_score,
      data = score_table[
        complete,
        ,
        drop = FALSE
      ]
    )

    adjusted[
      complete
    ] <- residuals(
      fit
    )

    coefficient_table <- data.frame(
      term =
        names(
          coef(fit)
        ),
      estimate =
        as.numeric(
          coef(fit)
        ),
      stringsAsFactors = FALSE
    )

    write_tsv(
      coefficient_table,
      file.path(
        OUT,
        "37I_epithelial_adjustment_model.tsv"
      )
    )
  }

  score_table$
    LAMC2_state_epithelial_adjusted <-
    adjusted

  score_table$
    LAMC2_state_epithelial_adjusted_z <-
    safe_z(
      adjusted
    )

  valid_state <- is.finite(
    score_table$
      LAMC2_state_rank_score
  )

  state_percentile <- rep(
    NA_real_,
    nrow(score_table)
  )

  state_percentile[
    valid_state
  ] <- (
    rank(
      score_table$
        LAMC2_state_rank_score[
          valid_state
        ],
      ties.method = "average"
    ) -
      0.5
  ) /
    sum(
      valid_state
    )

  score_table$
    LAMC2_state_rank_percentile <-
    state_percentile

  valid_adjusted <- is.finite(
    score_table$
      LAMC2_state_epithelial_adjusted
  )

  adjusted_percentile <- rep(
    NA_real_,
    nrow(score_table)
  )

  adjusted_percentile[
    valid_adjusted
  ] <- (
    rank(
      score_table$
        LAMC2_state_epithelial_adjusted[
          valid_adjusted
        ],
      ties.method = "average"
    ) -
      0.5
  ) /
    sum(
      valid_adjusted
    )

  score_table$
    LAMC2_state_adjusted_percentile <-
    adjusted_percentile

  score_table
}

make_correlation_table <- function(
    score_table
) {
  target_scores <- c(
    "LAMC2_state_rank_score",
    "LAMC2_state_epithelial_adjusted"
  )

  features <- c(
    "LAMC2_expression_processed",
    "LAMC2_within_sample_rank",
    "ITGB4_expression_processed",
    "ITGB4_within_sample_rank",
    "epithelial_rank_score",
    "Basal_rank_score",
    "Luminal_rank_score",
    "Adhesion_rank_score",
    "Invasion_rank_score",
    "Cell_cycle_rank_score"
  )

  rows <- list()
  row_index <- 0L

  for (target in target_scores) {
    for (feature in features) {
      x <- score_table[[target]]
      y <- score_table[[feature]]

      valid <- is.finite(x) &
        is.finite(y)

      row_index <- row_index + 1L

      if (sum(valid) >= 5L) {
        test <- suppressWarnings(
          cor.test(
            x[valid],
            y[valid],
            method = "spearman",
            exact = FALSE
          )
        )

        rho <- unname(
          test$estimate
        )

        p_value <- test$p.value
      } else {
        rho <- NA_real_
        p_value <- NA_real_
      }

      rows[[row_index]] <- data.frame(
        target_score =
          target,
        feature =
          feature,
        n_complete =
          sum(valid),
        spearman_rho =
          rho,
        p_value =
          p_value,
        stringsAsFactors = FALSE
      )
    }
  }

  result <- do.call(
    rbind,
    rows
  )

  result$FDR_BH <- p.adjust(
    result$p_value,
    method = "BH"
  )

  result
}

merge_clinical <- function(
    score_table
) {
  clinical <- data.table::fread(
    CLINICAL_FILE,
    data.table = FALSE,
    check.names = FALSE,
    na.strings = c(
      "NA",
      "NaN",
      ""
    )
  )

  candidate_columns <- c(
    "patient_id",
    "sample_id",
    "submitter_id",
    "bcr_patient_barcode",
    "case_submitter_id"
  )

  id_column <- candidate_columns[
    candidate_columns %in%
      colnames(clinical)
  ]

  if (length(id_column) == 0L) {
    fail(
      "No recognizable patient ID column in the clinical file."
    )
  }

  clinical$patient_id_for_merge <-
    patient_id_from_barcode(
      clinical[[id_column[1L]]]
    )

  clinical <- clinical[
    !is.na(
      clinical$
        patient_id_for_merge
    ) &
      nzchar(
        clinical$
          patient_id_for_merge
      ),
    ,
    drop = FALSE
  ]

  clinical <- clinical[
    !duplicated(
      clinical$
        patient_id_for_merge
    ),
    ,
    drop = FALSE
  ]

  merged <- merge(
    score_table,
    clinical,
    by.x = "patient_id",
    by.y = "patient_id_for_merge",
    all = FALSE,
    sort = FALSE
  )

  merged <- merged[
    match(
      score_table$patient_id,
      merged$patient_id
    ),
    ,
    drop = FALSE
  ]

  merged <- merged[
    !is.na(
      merged$patient_id
    ),
    ,
    drop = FALSE
  ]

  list(
    merged = merged,
    clinical_id_column =
      id_column[1L],
    clinical_unique_patients =
      nrow(clinical)
  )
}

main <- function() {
  writeLines(
    "RUNNING",
    STATUS_FILE
  )

  log_msg(
    "============================================================"
  )

  log_msg(
    "37K1 P01: LAMC2 excluded, ITGB4 eligible"
  )

  log_msg(
    "Start: ",
    format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S"
    )
  )

  log_msg(
    "============================================================"
  )

  required_files <- c(
    EXPR_FILE,
    CLINICAL_FILE,
    REFERENCE_FILE,
    LABEL_FILE
  )

  missing_files <- required_files[
    !file.exists(
      required_files
    )
  ]

  if (length(missing_files) > 0L) {
    fail(
      "Missing required file: ",
      paste(
        missing_files,
        collapse = " | "
      )
    )
  }

  bulk_result <- read_processed_bulk()

  bulk_expression <-
    bulk_result$expression

  log_msg(
    "Reading the 37B single-cell reference."
  )

  reference <- readRDS(
    REFERENCE_FILE
  )

  labels <- readRDS(
    LABEL_FILE
  )

  labels <- as.data.frame(
    labels,
    stringsAsFactors = FALSE
  )

  oriented <- orient_reference(
    reference =
      reference,
    labels =
      labels
  )

  reference <- oriented$reference
  labels <- oriented$labels

  common_reference_bulk_genes <- intersect(
    rownames(reference),
    rownames(bulk_expression)
  )

  if (
    length(
      common_reference_bulk_genes
    ) < 3000L
  ) {
    fail(
      "The 37B reference and original bulk matrix share only ",
      length(
        common_reference_bulk_genes
      ),
      " genes."
    )
  }

  reference <- reference[
    common_reference_bulk_genes,
    ,
    drop = FALSE
  ]

  bulk_for_markers <- bulk_expression[
    common_reference_bulk_genes,
    ,
    drop = FALSE
  ]

  log_msg(
    "Number of shared genes: ",
    length(
      common_reference_bulk_genes
    )
  )

  state_labels <- labels[
    labels$cell_state %in%
      c(
        "Epi_LAMC2_high",
        "Epi_LAMC2_low_detected",
        "Epi_LAMC2_negative"
      ),
    ,
    drop = FALSE
  ]

  state_reference <- reference[
    ,
    match(
      state_labels$cell_id,
      colnames(reference)
    ),
    drop = FALSE
  ]

  state_labels$contrast_group <- ifelse(
    state_labels$cell_state ==
      "Epi_LAMC2_high",
    "Epi_LAMC2_high",
    "Epi_LAMC2_non_high"
  )

  state_markers <- derive_paired_markers(
    reference =
      state_reference,
    labels =
      state_labels,
    positive_label =
      "Epi_LAMC2_high",
    negative_label =
      "Epi_LAMC2_non_high",
    contrast_name =
      "LAMC2_high_vs_non_high_within_epithelial",
    bulk_expression =
      bulk_for_markers,
    excluded_genes =
      c(
        "LAMC2"
      )
  )

  write_tsv(
    state_markers$
      marker_table,
    STATE_MARKER_FILE
  )

  epithelial_labels <- labels

  epithelial_labels$contrast_group <-
    ifelse(
      epithelial_labels$cell_type ==
        "Epithelial",
      "Epithelial",
      "Non_epithelial"
    )

  epithelial_markers <-
    derive_paired_markers(
      reference =
        reference,
      labels =
        epithelial_labels,
      positive_label =
        "Epithelial",
      negative_label =
        "Non_epithelial",
      contrast_name =
        "epithelial_vs_non_epithelial",
      bulk_expression =
        bulk_for_markers,
      excluded_genes =
        c(
          "LAMC2"
        )
    )

  write_tsv(
    epithelial_markers$
      marker_table,
    EPI_MARKER_FILE
  )

  marker_qc <- data.frame(
    contrast = c(
      "LAMC2_high_vs_non_high_within_epithelial",
      "epithelial_vs_non_epithelial"
    ),
    paired_reference_patients = c(
      length(
        state_markers$
          paired_patients
      ),
      length(
        epithelial_markers$
          paired_patients
      )
    ),
    positive_markers = c(
      length(
        state_markers$
          positive_genes
      ),
      length(
        epithelial_markers$
          positive_genes
      )
    ),
    negative_markers = c(
      length(
        state_markers$
          negative_genes
      ),
      length(
        epithelial_markers$
          negative_genes
      )
    ),
    LAMC2_in_positive_markers = c(
      "LAMC2" %in%
        state_markers$
          positive_genes,
      "LAMC2" %in%
        epithelial_markers$
          positive_genes
    ),
    ITGB4_in_state_markers = c(
      "ITGB4" %in%
        c(
          state_markers$
            positive_genes,
          state_markers$
            negative_genes
        ),
      NA
    ),
    stringsAsFactors = FALSE
  )

  write_tsv(
    marker_qc,
    MARKER_QC_FILE
  )

  log_msg(
    "Starting within-sample rank projection for 408 samples."
  )

  score_table <- score_bulk_samples(
    bulk_expression =
      bulk_expression,
    state_positive =
      state_markers$
        positive_genes,
    state_negative =
      state_markers$
        negative_genes,
    epi_positive =
      epithelial_markers$
        positive_genes,
    epi_negative =
      epithelial_markers$
        negative_genes
  )

  score_table <- make_adjusted_scores(
    score_table
  )

  if (nrow(score_table) != EXPECTED_TUMORS) {
    fail(
      "Expected the final score table to contain ",
      EXPECTED_TUMORS,
      ", but found ",
      nrow(score_table),
      "."
    )
  }

  if (anyDuplicated(score_table$patient_id)) {
    fail(
      "Duplicate patients in the final score table."
    )
  }

  write_tsv(
    score_table,
    SCORE_FILE
  )

  saveRDS(
    score_table,
    SCORE_RDS,
    compress = "gzip"
  )

  correlations <- make_correlation_table(
    score_table
  )

  write_tsv(
    correlations,
    CORRELATION_FILE
  )

  clinical_result <- merge_clinical(
    score_table
  )

  write_tsv(
    clinical_result$merged,
    CLINICAL_SCORE_FILE
  )

  ready <- (
    nrow(score_table) ==
      EXPECTED_TUMORS &&
      length(
        state_markers$
          positive_genes
      ) >=
        MIN_MARKERS_PER_DIRECTION &&
      length(
        state_markers$
          negative_genes
      ) >=
        MIN_MARKERS_PER_DIRECTION &&
      !"LAMC2" %in%
        c(
          state_markers$
            positive_genes,
          state_markers$
            negative_genes
        ) &&
      nrow(
        clinical_result$merged
      ) >=
        400L
  )

  review <- c(
    "37K1 P01 state rank projection completed (LAMC2 excluded, ITGB4 eligible).",
    "",
    "Input matrix: ",
    EXPR_FILE,
    "",
    paste0(
      "Number of primary tumor samples: ",
      nrow(score_table)
    ),
    paste0(
      "Number of clinically matched patients: ",
      nrow(
        clinical_result$merged
      )
    ),
    paste0(
      "Clinical patient ID column: ",
      clinical_result$
        clinical_id_column
    ),
    "",
    "Expression matrix processing principles:",
    "- Neither 2^x nor 2^x-1 was applied.",
    "- Negative values were not shifted.",
    "- Negative values were not truncated to zero.",
    "- Continuous expression values were not rounded to counts.",
    "- Missing values were not globally replaced with zero.",
    "- Only finite expression values were ranked within each sample.",
    "",
    "Analysis name:",
    "- single-cell-guided within-sample rank projection",
    "- Single-cell reference-guided within-sample rank projection",
    "",
    "Interpretation:",
    "- LAMC2_state_rank_score is a relative state enrichment score.",
    "- LAMC2_state_epithelial_adjusted is a relative score adjusted for epithelial enrichment.",
    "- Neither score represents a cell proportion or a BayesPrism fraction.",
    "",
    paste0(
      "Number of LAMC2-high positive markers: ",
      length(
        state_markers$
          positive_genes
      )
    ),
    paste0(
      "Number of LAMC2-non-high negative markers: ",
      length(
        state_markers$
          negative_genes
      )
    ),
    paste0(
      "Number of paired reference patients for state markers: ",
      length(
        state_markers$
          paired_patients
      )
    ),
    paste0(
      "Number of epithelial positive markers: ",
      length(
        epithelial_markers$
          positive_genes
      )
    ),
    paste0(
      "Number of non-epithelial negative markers: ",
      length(
        epithelial_markers$
          negative_genes
      )
    ),
    "",
    paste0(
      "LAMC2 included in the primary state markers: ",
      "LAMC2" %in%
        c(
          state_markers$
            positive_genes,
          state_markers$
            negative_genes
        )
    ),
    paste0(
      "ITGB4 included in the primary state markers: ",
      "ITGB4" %in%
        c(
          state_markers$
            positive_genes,
          state_markers$
            negative_genes
        )
    ),
    "",
    paste0(
      "Final QC passed: ",
      ready
    ),
    "",
    "Main outputs:",
    SCORE_FILE,
    CLINICAL_SCORE_FILE,
    STATE_MARKER_FILE,
    EPI_MARKER_FILE,
    CORRELATION_FILE,
    "",
    "No data were downloaded in this step.",
    "BayesPrism was not run in this step.",
    "Processed expression values were not treated as raw counts."
  )

  writeLines(
    review,
    REVIEW_FILE,
    useBytes = TRUE
  )

  writeLines(
    capture.output(
      sessionInfo()
    ),
    file.path(
      OUT,
      "37I_sessionInfo.txt"
    )
  )

  if (ready) {
    writeLines(
      "SUCCESS_PROCESSED_BULK_RANK_PROJECTION",
      STATUS_FILE
    )
  } else {
    writeLines(
      "FAILED_FINAL_QC",
      STATUS_FILE
    )
  }

  log_msg(
    "============================================================"
  )

  log_msg(
    "37I completed. Final QC: ",
    ready
  )

  log_msg(
    "End: ",
    format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S"
    )
  )

  log_msg(
    "============================================================"
  )
}

tryCatch(
  main(),
  error = function(e) {
    error_text <- conditionMessage(e)

    writeLines(
      paste0(
        "FAILED: ",
        error_text
      ),
      STATUS_FILE
    )

    writeLines(
      c(
        "37I run failed.",
        paste0(
          "Error: ",
          error_text
        ),
        "",
        "No data were downloaded.",
        "The original expression matrix was not modified.",
        "BayesPrism was not run.",
        "See 37I_run_log.txt."
      ),
      REVIEW_FILE,
      useBytes = TRUE
    )

    log_msg(
      "[ERROR] ",
      error_text
    )

    stop(
      error_text,
      call. = FALSE
    )
  }
)
