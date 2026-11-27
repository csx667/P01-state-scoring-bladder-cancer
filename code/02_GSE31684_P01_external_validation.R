options(
  stringsAsFactors = FALSE,
  warn = 1,
  width = 180
)

BASE <- "/root/autodl-tmp/scrna_project"
RESULT_ROOT <- file.path(
  BASE,
  "result",
  "37_BayesPrism_rebuild_v5",
  "37K8_GEO_P01_external_validation",
  "01_GSE31684_FIRST_PART"
)
RUN_TAG <- format(Sys.time(), "%Y%m%d_%H%M%S")
OUT <- file.path(RESULT_ROOT, paste0("run_", RUN_TAG))

DIRS <- list(
  audit = file.path(OUT, "00_input_audit"),
  projection = file.path(OUT, "01_projection_main_median"),
  pathology = file.path(OUT, "02_pathology_main"),
  survival = file.path(OUT, "03_survival_main"),
  figures_main = file.path(OUT, "04_figures_main"),
  sensitivity = file.path(OUT, "05_sensitivity"),
  tables = file.path(OUT, "06_tables"),
  report = file.path(OUT, "07_report"),
  code = file.path(OUT, "08_code")
)

SENS_DIRS <- list(
  iqr_probe = file.path(DIRS$sensitivity, "01_IQR_highest_probe"),
  P00 = file.path(DIRS$sensitivity, "02_policy_P00"),
  P11 = file.path(DIRS$sensitivity, "03_policy_P11"),
  strict_mibc = file.path(DIRS$sensitivity, "04_strict_MIBC_pT2_pT4"),
  n_known = file.path(DIRS$sensitivity, "05_N_known"),
  pure_tcc = file.path(DIRS$sensitivity, "06_pure_TCC"),
  no_preop = file.path(DIRS$sensitivity, "07_no_preop_chemo"),
  nonph = file.path(DIRS$sensitivity, "08_nonPH_timevarying"),
  top10 = file.path(DIRS$sensitivity, "09_Top10_descriptive")
)

for (path in c(OUT, unlist(DIRS, use.names = FALSE), unlist(SENS_DIRS, use.names = FALSE))) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(path, "tables"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(path, "figures"), recursive = TRUE, showWarnings = FALSE)
}

dir.create(RESULT_ROOT, recursive = TRUE, showWarnings = FALSE)
writeLines(OUT, file.path(RESULT_ROOT, "LATEST_RUN_PATH.txt"), useBytes = TRUE)

STATUS_FILE <- file.path(OUT, "GSE31684_FIRST_PART_RUN_STATUS.txt")
REPORT_FILE <- file.path(OUT, "GSE31684_FIRST_PART_SEND_TO_CHATGPT.txt")
LOG_FILE <- file.path(OUT, "GSE31684_FIRST_PART_run_log.txt")

log_con <- file(LOG_FILE, open = "wt")

log_msg <- function(...) {
  text <- paste0(..., collapse = "")
  cat(text, "\n")
  cat(text, "\n", file = log_con)
  flush(log_con)
}

fail <- function(...) {
  text <- paste0(..., collapse = "")
  log_msg("[FATAL] ", text)
  stop(text, call. = FALSE)
}

write_tsv <- function(x, path) {
  if (is.null(x)) {
    x <- data.frame(message = "NULL_OBJECT", stringsAsFactors = FALSE)
  }
  write.table(
    x,
    path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = "NA",
    fileEncoding = "UTF-8"
  )
}

bind_rows_fill <- function(items) {
  items <- items[vapply(items, function(x) is.data.frame(x) && nrow(x) > 0L, logical(1))]
  if (length(items) == 0L) {
    return(data.frame())
  }
  all_names <- unique(unlist(lapply(items, names), use.names = FALSE))
  filled <- lapply(
    items,
    function(x) {
      missing <- setdiff(all_names, names(x))
      for (name in missing) {
        x[[name]] <- NA
      }
      x[, all_names, drop = FALSE]
    }
  )
  result <- do.call(rbind, filled)
  rownames(result) <- NULL
  result
}

safe_read_tsv <- function(path) {
  if (!file.exists(path)) {
    fail("TSV file not found: ", path)
  }
  read.delim(
    path,
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )
}

locate_data_file <- function(exact_basename, label = exact_basename) {
  path <- file.path(BASE, "data", exact_basename)
  if (!file.exists(path)) {
    fail(
      "Required file not found in the data directory: ",
      label,
      "; expected path=",
      path
    )
  }
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

read_gse31684_series_matrix <- function(path) {
  if (!file.exists(path)) {
    fail("GSE31684 Series Matrix not found: ", path)
  }

  connection <- if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    gzfile(path, open = "rt")
  } else {
    file(path, open = "rt")
  }

  on.exit(
    try(close(connection), silent = TRUE),
    add = TRUE
  )

  lines <- readLines(
    connection,
    warn = FALSE,
    encoding = "UTF-8"
  )

  begin_index <- which(trimws(lines) == "!series_matrix_table_begin")
  end_index <- which(trimws(lines) == "!series_matrix_table_end")

  if (length(begin_index) != 1L || length(end_index) != 1L) {
    fail(
      "Series Matrix table boundaries are not unique: begin=",
      length(begin_index),
      "; end=",
      length(end_index)
    )
  }

  if (end_index <= begin_index + 1L) {
    fail("The Series Matrix table is empty or its boundaries are out of order.")
  }

  table_lines <- lines[(begin_index + 1L):(end_index - 1L)]
  text_connection <- textConnection(table_lines)
  on.exit(
    try(close(text_connection), silent = TRUE),
    add = TRUE
  )

  table <- read.delim(
    text_connection,
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "\"",
    comment.char = "",
    na.strings = c("NA", "NaN", "null", "NULL", "")
  )

  if (ncol(table) < 2L || nrow(table) < 1000L) {
    fail(
      "Unexpected dimensions after reading the Series Matrix: ",
      nrow(table),
      "×",
      ncol(table)
    )
  }

  probe_column <- names(table)[1L]
  probes <- clean_character(table[[probe_column]])
  if (any(is.na(probes)) || anyDuplicated(probes)) {
    fail("Missing or duplicate probe IDs in the Series Matrix.")
  }

  expression_data <- table[, -1L, drop = FALSE]
  expression_data[] <- lapply(
    expression_data,
    function(x) suppressWarnings(as.numeric(x))
  )
  expression_matrix <- as.matrix(expression_data)
  storage.mode(expression_matrix) <- "double"
  rownames(expression_matrix) <- probes
  colnames(expression_matrix) <- toupper(clean_character(names(expression_data)))

  if (any(is.na(colnames(expression_matrix))) || anyDuplicated(colnames(expression_matrix))) {
    fail("Missing or duplicate GSM column names in the Series Matrix.")
  }

  metadata_lines <- lines[seq_len(begin_index - 1L)]
  platform_lines <- metadata_lines[grepl("^!Series_platform_id", metadata_lines)]
  sample_platform_lines <- metadata_lines[grepl("^!Sample_platform_id", metadata_lines)]

  audit <- data.frame(
    item = c(
      "source_file",
      "compressed_size_bytes",
      "total_text_lines",
      "table_begin_line",
      "table_end_line",
      "probe_rows",
      "sample_columns",
      "first_probe",
      "last_probe",
      "first_sample",
      "last_sample",
      "series_platform_line",
      "sample_platform_line"
    ),
    value = c(
      normalizePath(path, winslash = "/", mustWork = TRUE),
      as.character(file.info(path)$size),
      as.character(length(lines)),
      as.character(begin_index),
      as.character(end_index),
      as.character(nrow(expression_matrix)),
      as.character(ncol(expression_matrix)),
      rownames(expression_matrix)[1L],
      rownames(expression_matrix)[nrow(expression_matrix)],
      colnames(expression_matrix)[1L],
      colnames(expression_matrix)[ncol(expression_matrix)],
      if (length(platform_lines) > 0L) platform_lines[1L] else "NOT_FOUND",
      if (length(sample_platform_lines) > 0L) sample_platform_lines[1L] else "NOT_FOUND"
    ),
    stringsAsFactors = FALSE
  )

  list(
    matrix = expression_matrix,
    audit = audit,
    metadata_lines = metadata_lines
  )
}

normalize_name <- function(x) {
  gsub("[^a-z0-9]+", "", tolower(as.character(x)))
}

pick_column <- function(data, candidates, required = TRUE, label = "column") {
  normalized <- normalize_name(names(data))
  for (candidate in candidates) {
    hit <- which(normalized == normalize_name(candidate))
    if (length(hit) == 1L) {
      return(names(data)[hit])
    }
    if (length(hit) > 1L) {
      fail("Field match is not unique: ", label, "; candidates=", candidate, "; matches=", paste(names(data)[hit], collapse = "|"))
    }
  }
  if (required) {
    fail("Required field not found: ", label, "; candidates=", paste(candidates, collapse = "|"))
  }
  NA_character_
}

clean_character <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "N/A", "NULL", "null", ".")] <- NA_character_
  x
}

parse_binary_event <- function(x) {
  raw <- tolower(clean_character(x))
  result <- rep(NA_integer_, length(raw))
  numeric_value <- suppressWarnings(as.numeric(raw))
  result[is.finite(numeric_value) & numeric_value == 0] <- 0L
  result[is.finite(numeric_value) & numeric_value == 1] <- 1L
  result[grepl("dead|death|deceased|event|yes|positive|dod", raw)] <- 1L
  result[grepl("alive|censor|disease free", raw)] <- 0L
  result[raw %in% c("no", "negative", "none")] <- 0L
  result
}

parse_yes_no <- function(x) {
  raw <- tolower(clean_character(x))
  result <- rep(NA_character_, length(raw))
  result[raw %in% c("1", "yes", "y", "true", "positive", "received")] <- "Yes"
  result[raw %in% c("0", "no", "n", "false", "negative", "none", "not received")] <- "No"
  factor(result, levels = c("No", "Yes"))
}

standardize_sex <- function(x) {
  raw <- tolower(clean_character(x))
  raw_key <- gsub("[^a-z0-9]", "", raw)
  result <- rep(NA_character_, length(raw_key))

  female_values <- c("f", "female", "woman", "women")
  male_values <- c("m", "male", "man", "men")

  result[raw_key %in% female_values] <- "Female"
  result[raw_key %in% male_values] <- "Male"

  factor(result, levels = c("Female", "Male"))
}

standardize_t_stage4 <- function(x) {
  raw <- gsub("[^A-Z0-9]", "", toupper(clean_character(x)))
  result <- rep(NA_character_, length(raw))
  result[grepl("^PTA|^TA$|^PT1|^T1", raw)] <- "pTa/pT1"
  result[grepl("^PT2|^T2", raw)] <- "pT2"
  result[grepl("^PT3|^T3", raw)] <- "pT3"
  result[grepl("^PT4|^T4", raw)] <- "pT4"
  factor(result, levels = c("pTa/pT1", "pT2", "pT3", "pT4"), ordered = TRUE)
}

standardize_node <- function(x) {
  raw <- gsub("[^A-Z0-9+-]", "", toupper(clean_character(x)))
  result <- rep(NA_character_, length(raw))
  result[grepl("UNKNOWN|^NX$|^NAX$|^NA$", raw)] <- "Unknown"
  result[grepl("NEGATIVE|^N0$|^0$|^N-$", raw)] <- "Negative"
  result[grepl("POSITIVE|^N[1234]|^N\\+$|^[1234]$", raw)] <- "Positive"
  factor(result, levels = c("Negative", "Positive", "Unknown"))
}

z_score <- function(x) {
  x <- as.numeric(x)
  sd_value <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(sd_value) || sd_value <= 0) {
    return(rep(NA_real_, length(x)))
  }
  as.numeric((x - mean(x, na.rm = TRUE)) / sd_value)
}

make_quartile <- function(x) {
  result <- rep(NA_character_, length(x))
  valid <- which(is.finite(x))
  if (length(valid) == 0L) {
    return(factor(result, levels = paste0("Q", 1:4)))
  }
  ord <- order(x[valid], valid, na.last = NA)
  rank_index <- integer(length(valid))
  rank_index[ord] <- seq_along(valid)
  quartile <- pmin(4L, ceiling(rank_index * 4 / length(valid)))
  result[valid] <- paste0("Q", quartile)
  factor(result, levels = paste0("Q", 1:4))
}

make_top10 <- function(x) {
  result <- rep(NA_character_, length(x))
  valid <- which(is.finite(x))
  if (length(valid) == 0L) {
    return(factor(result, levels = c("Other90", "Top10")))
  }
  cutoff <- as.numeric(stats::quantile(x[valid], probs = 0.90, na.rm = TRUE, type = 7))
  result[valid] <- ifelse(x[valid] >= cutoff, "Top10", "Other90")
  factor(result, levels = c("Other90", "Top10"))
}

format_p <- function(p) {
  if (!is.finite(p)) {
    return("NA")
  }
  if (p < 0.001) {
    return(format(p, scientific = TRUE, digits = 2))
  }
  sprintf("%.3f", p)
}

SEARCH_ROOTS <- unique(c(
  file.path(BASE, "data"),
  file.path(BASE, "result"),
  "/data",
  file.path(OUT, "00_input_audit", "extracted_inputs")
))
SEARCH_ROOTS <- SEARCH_ROOTS[dir.exists(SEARCH_ROOTS)]

extract_known_archives <- function() {
  initial_roots <- unique(c(file.path(BASE, "data"), file.path(BASE, "result"), "/data"))
  initial_roots <- initial_roots[dir.exists(initial_roots)]
  zip_files <- unique(unlist(lapply(initial_roots, list.files, recursive = TRUE, full.names = TRUE, pattern = "\\.zip$", ignore.case = TRUE), use.names = FALSE))
  selected <- zip_files[grepl("AUTO_DL_external_validation_starter|37K7_CLEAN_REVIEW_BUNDLE", basename(zip_files), ignore.case = TRUE)]
  extraction_dir <- file.path(OUT, "00_input_audit", "extracted_inputs")
  dir.create(extraction_dir, recursive = TRUE, showWarnings = FALSE)
  rows <- list()
  for (i in seq_along(selected)) {
    archive <- selected[i]
    target <- file.path(extraction_dir, paste0("archive_", i))
    dir.create(target, recursive = TRUE, showWarnings = FALSE)
    status <- tryCatch(
      {
        utils::unzip(archive, exdir = target)
        "SUCCESS"
      },
      error = function(e) paste0("FAILED: ", conditionMessage(e))
    )
    rows[[i]] <- data.frame(archive = archive, extracted_to = target, status = status, stringsAsFactors = FALSE)
  }
  bind_rows_fill(rows)
}

archive_audit <- extract_known_archives()
write_tsv(archive_audit, file.path(DIRS$audit, "A00_extracted_known_archives.tsv"))
SEARCH_ROOTS <- unique(c(SEARCH_ROOTS, file.path(OUT, "00_input_audit", "extracted_inputs")))
SEARCH_ROOTS <- SEARCH_ROOTS[dir.exists(SEARCH_ROOTS)]

inventory_files <- function() {
  files <- unique(unlist(lapply(SEARCH_ROOTS, list.files, recursive = TRUE, full.names = TRUE, include.dirs = FALSE), use.names = FALSE))
  info <- file.info(files)
  data.frame(
    path = files,
    basename = basename(files),
    size_bytes = as.numeric(info$size),
    modified = as.character(info$mtime),
    stringsAsFactors = FALSE
  )
}

FILE_INVENTORY <- inventory_files()
write_tsv(FILE_INVENTORY, file.path(DIRS$audit, "A01_all_candidate_file_inventory.tsv"))

choose_file <- function(exact_basename, preferred_path_pattern = NULL, required = TRUE, label = exact_basename) {
  candidates <- FILE_INVENTORY[FILE_INVENTORY$basename == exact_basename, , drop = FALSE]
  if (nrow(candidates) == 0L) {
    if (required) {
      fail("Required file not found: ", label, "; exact filename=", exact_basename)
    }
    return(NA_character_)
  }
  score <- rep(0, nrow(candidates))
  if (!is.null(preferred_path_pattern)) {
    score <- score + 1000 * grepl(preferred_path_pattern, candidates$path, ignore.case = TRUE)
  }
  score <- score + log10(pmax(candidates$size_bytes, 1))
  candidates$selection_score <- score
  candidates <- candidates[order(-candidates$selection_score, candidates$path), , drop = FALSE]
  write_tsv(candidates, file.path(DIRS$audit, paste0("A_candidate_", gsub("[^A-Za-z0-9]+", "_", label), ".tsv")))
  candidates$path[1L]
}

find_files_pattern <- function(pattern) {
  FILE_INVENTORY$path[grepl(pattern, FILE_INVENTORY$path, ignore.case = TRUE)]
}

load_marker_definitions <- function() {
  clean_file <- choose_file(
    "37K7_CLEAN_POLICY_MARKERS.tsv",
    preferred_path_pattern = "37K7",
    required = FALSE,
    label = "37K7 clean policy markers"
  )
  if (!is.na(clean_file)) {
    table <- safe_read_tsv(clean_file)
    return(list(table = table, source = clean_file))
  }

  rows <- list()
  policies <- c("P00", "P01", "P11")
  for (policy in policies) {
    candidates <- FILE_INVENTORY$path[
      FILE_INVENTORY$basename == "37I_LAMC2_high_vs_non_high_marker_table.tsv" &
        grepl(paste0("37K1_", policy, "_"), FILE_INVENTORY$path, ignore.case = TRUE)
    ]
    if (length(candidates) != 1L) {
      fail("Unable to uniquely locate the ", policy, " marker table; number of candidates=", length(candidates), "; ", paste(candidates, collapse = " | "))
    }
    one <- safe_read_tsv(candidates[1L])
    if (!"policy" %in% names(one)) {
      one$policy <- policy
    }
    rows[[policy]] <- one
  }
  list(table = bind_rows_fill(rows), source = paste(vapply(rows, function(x) "source_marker_table", character(1)), collapse = "|"))
}

validate_marker_definitions <- function(marker_table) {
  policy_col <- pick_column(marker_table, c("policy"), TRUE, "marker policy")
  gene_col <- pick_column(marker_table, c("gene", "gene_symbol", "symbol"), TRUE, "marker gene")
  direction_col <- pick_column(marker_table, c("selected_direction", "direction"), TRUE, "marker direction")
  clean <- data.frame(
    policy = toupper(clean_character(marker_table[[policy_col]])),
    gene = toupper(clean_character(marker_table[[gene_col]])),
    direction = tolower(clean_character(marker_table[[direction_col]])),
    stringsAsFactors = FALSE
  )
  clean <- clean[
    clean$policy %in% c("P00", "P01", "P11") &
      clean$direction %in% c("positive_marker", "negative_marker") &
      !is.na(clean$gene),
    ,
    drop = FALSE
  ]
  clean <- unique(clean)

  summary_rows <- list()
  for (policy in c("P00", "P01", "P11")) {
    one <- clean[clean$policy == policy, , drop = FALSE]
    pos <- one$gene[one$direction == "positive_marker"]
    neg <- one$gene[one$direction == "negative_marker"]
    summary_rows[[policy]] <- data.frame(
      policy = policy,
      positive_markers = length(unique(pos)),
      negative_markers = length(unique(neg)),
      LAMC2_direction = ifelse("LAMC2" %in% pos, "positive_marker", ifelse("LAMC2" %in% neg, "negative_marker", "not_selected")),
      ITGB4_direction = ifelse("ITGB4" %in% pos, "positive_marker", ifelse("ITGB4" %in% neg, "negative_marker", "not_selected")),
      stringsAsFactors = FALSE
    )
  }
  summary <- bind_rows_fill(summary_rows)
  write_tsv(summary, file.path(DIRS$audit, "A02_marker_policy_validation.tsv"))

  if (any(summary$positive_markers != 120L) || any(summary$negative_markers != 120L)) {
    fail("P00/P01/P11 marker counts failed the 120 positive + 120 negative check.")
  }
  p00 <- summary[summary$policy == "P00", , drop = FALSE]
  p01 <- summary[summary$policy == "P01", , drop = FALSE]
  p11 <- summary[summary$policy == "P11", , drop = FALSE]
  if (p00$LAMC2_direction != "not_selected" || p00$ITGB4_direction != "not_selected") {
    fail("P00 policy check failed: both LAMC2 and ITGB4 must be excluded.")
  }
  if (p01$LAMC2_direction != "not_selected" || p01$ITGB4_direction != "positive_marker") {
    fail("P01 policy check failed: LAMC2 must be excluded and ITGB4 must be a positive marker.")
  }
  if (p11$LAMC2_direction != "positive_marker" || p11$ITGB4_direction != "positive_marker") {
    fail("P11 policy check failed: both LAMC2 and ITGB4 must be positive markers.")
  }
  clean
}

find_expression_matrix <- function(root_object) {
  preferred <- c("expression_matrix", "expr", "expr_matrix", "expression")
  for (name in preferred) {
    if (is.list(root_object) && name %in% names(root_object)) {
      candidate <- root_object[[name]]
      if (is.matrix(candidate) || is.data.frame(candidate)) {
        matrix_candidate <- as.matrix(candidate)
        if (nrow(matrix_candidate) > 1000L && ncol(matrix_candidate) >= 50L) {
          return(list(matrix = matrix_candidate, component = name))
        }
      }
    }
  }
  if (is.list(root_object)) {
    for (name in names(root_object)) {
      candidate <- root_object[[name]]
      if (is.matrix(candidate) || is.data.frame(candidate)) {
        matrix_candidate <- as.matrix(candidate)
        if (nrow(matrix_candidate) > 1000L && ncol(matrix_candidate) >= 50L) {
          return(list(matrix = matrix_candidate, component = name))
        }
      }
    }
  }
  fail("No complete expression matrix with >1000 rows and >=50 columns found in the D98 object.")
}

find_mapping_in_object <- function(root_object) {
  if (!is.list(root_object)) {
    return(NULL)
  }
  for (name in names(root_object)) {
    candidate <- root_object[[name]]
    if (!is.data.frame(candidate) || nrow(candidate) < 1000L) {
      next
    }
    normalized <- normalize_name(names(candidate))
    has_probe <- any(normalized %in% normalize_name(c("probe", "probe_id", "probeset_id", "ID", "ID_REF"))) || !is.null(rownames(candidate))
    has_gene <- any(normalized %in% normalize_name(c("gene", "gene_symbol", "symbol", "Gene Symbol")))
    if (has_probe && has_gene) {
      return(list(table = candidate, component = name))
    }
  }
  NULL
}

read_annotation_file <- function(path) {
  connection <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, open = "rt") else file(path, open = "rt")
  preview <- readLines(connection, n = 60L, warn = FALSE)
  try(close(connection), silent = TRUE)
  if (length(preview) == 0L) {
    stop("empty annotation file")
  }
  header_index <- which(
    grepl("probe|probeset|id_ref|gene.symbol|gene_symbol|symbol", preview, ignore.case = TRUE) &
      (grepl("\\t", preview) | grepl(",", preview))
  )
  if (length(header_index) == 0L) {
    stop("no recognizable annotation header")
  }
  header_index <- header_index[1L]
  separator <- if (lengths(regmatches(preview[header_index], gregexpr("\\t", preview[header_index]))) >= lengths(regmatches(preview[header_index], gregexpr(",", preview[header_index])))) "\t" else ","
  if (separator == "\t") {
    read.delim(path, header = TRUE, skip = header_index - 1L, check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
  } else {
    read.csv(path, header = TRUE, skip = header_index - 1L, check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
  }
}

extract_probe_gene_mapping <- function(annotation_table) {
  probe_col <- pick_column(annotation_table, c("probe_id", "probeset_id", "probe set id", "ID_REF", "ID"), required = FALSE, label = "probe ID")
  gene_col <- pick_column(annotation_table, c("gene_symbol", "Gene Symbol", "symbol", "gene"), required = TRUE, label = "gene symbol")
  probe <- if (is.na(probe_col)) rownames(annotation_table) else annotation_table[[probe_col]]
  gene_raw <- annotation_table[[gene_col]]

  parse_one_gene <- function(value) {
    if (is.na(value) || !nzchar(trimws(as.character(value)))) {
      return(NA_character_)
    }
    value <- toupper(trimws(as.character(value)))
    value <- gsub('"', "", value, fixed = TRUE)
    pieces <- unlist(strsplit(value, "\\s*(///|//|;|,|\\|)\\s*", perl = TRUE), use.names = FALSE)
    pieces <- unique(trimws(pieces))
    pieces <- pieces[grepl("^[A-Z0-9][A-Z0-9.-]*$", pieces)]
    pieces <- pieces[!pieces %in% c("NA", "N/A", "---", "NULL")]
    if (length(pieces) != 1L) {
      return(NA_character_)
    }
    pieces[1L]
  }

  gene <- vapply(gene_raw, parse_one_gene, character(1))
  mapping <- data.frame(
    probe_id = clean_character(probe),
    gene = gene,
    stringsAsFactors = FALSE
  )
  mapping <- mapping[!is.na(mapping$probe_id) & !is.na(mapping$gene), , drop = FALSE]
  mapping <- unique(mapping)
  duplicate_probe <- duplicated(mapping$probe_id) | duplicated(mapping$probe_id, fromLast = TRUE)
  if (any(duplicate_probe)) {
    duplicate_ids <- unique(mapping$probe_id[duplicate_probe])
    keep_ids <- duplicate_ids[vapply(
      duplicate_ids,
      function(id) length(unique(mapping$gene[mapping$probe_id == id])) == 1L,
      logical(1)
    )]
    mapping <- mapping[!mapping$probe_id %in% duplicate_ids | mapping$probe_id %in% keep_ids, , drop = FALSE]
    mapping <- unique(mapping)
  }
  mapping
}

load_probe_mapping <- function(root_object, expression_probes) {
  object_mapping <- find_mapping_in_object(root_object)
  candidates_log <- list()
  if (!is.null(object_mapping)) {
    mapping <- tryCatch(extract_probe_gene_mapping(object_mapping$table), error = function(e) NULL)
    if (!is.null(mapping) && sum(expression_probes %in% mapping$probe_id) >= 10000L) {
      return(list(mapping = mapping, source = paste0("D98 component: ", object_mapping$component)))
    }
  }

  annotation_candidates <- FILE_INVENTORY$path[
    grepl("GPL570|hgu133plus2|probe.*gene|probe.*symbol", FILE_INVENTORY$path, ignore.case = TRUE) &
      grepl("\\.(tsv|csv|txt|gz)$", FILE_INVENTORY$path, ignore.case = TRUE) &
      !grepl("series_matrix|clinical|family.soft", FILE_INVENTORY$path, ignore.case = TRUE)
  ]
  annotation_candidates <- annotation_candidates[order(-file.info(annotation_candidates)$size)]

  for (path in annotation_candidates) {
    result <- tryCatch(
      {
        table <- read_annotation_file(path)
        mapping <- extract_probe_gene_mapping(table)
        overlap <- sum(expression_probes %in% mapping$probe_id)
        candidates_log[[length(candidates_log) + 1L]] <- data.frame(path = path, rows = nrow(mapping), overlap = overlap, status = "PARSED", stringsAsFactors = FALSE)
        if (overlap >= 10000L) list(mapping = mapping, source = path) else NULL
      },
      error = function(e) {
        candidates_log[[length(candidates_log) + 1L]] <<- data.frame(path = path, rows = NA_integer_, overlap = NA_integer_, status = paste0("FAILED: ", conditionMessage(e)), stringsAsFactors = FALSE)
        NULL
      }
    )
    if (!is.null(result)) {
      write_tsv(bind_rows_fill(candidates_log), file.path(DIRS$audit, "A03_annotation_candidate_audit.tsv"))
      return(result)
    }
  }

  if (requireNamespace("AnnotationDbi", quietly = TRUE) && requireNamespace("hgu133plus2.db", quietly = TRUE)) {
    mapping_table <- AnnotationDbi::select(
      hgu133plus2.db::hgu133plus2.db,
      keys = unique(expression_probes),
      keytype = "PROBEID",
      columns = c("SYMBOL")
    )
    names(mapping_table)[names(mapping_table) == "PROBEID"] <- "probe_id"
    names(mapping_table)[names(mapping_table) == "SYMBOL"] <- "gene"
    mapping_table$probe_id <- clean_character(mapping_table$probe_id)
    mapping_table$gene <- toupper(clean_character(mapping_table$gene))
    mapping_table <- mapping_table[
      !is.na(mapping_table$probe_id) & !is.na(mapping_table$gene),
      c("probe_id", "gene"),
      drop = FALSE
    ]
    mapping_table <- unique(mapping_table)
    gene_count_by_probe <- tapply(
      mapping_table$gene,
      mapping_table$probe_id,
      function(x) length(unique(x))
    )
    unambiguous_probes <- names(gene_count_by_probe)[gene_count_by_probe == 1L]
    mapping_table <- mapping_table[
      mapping_table$probe_id %in% unambiguous_probes,
      ,
      drop = FALSE
    ]
    mapping_table <- mapping_table[!duplicated(mapping_table$probe_id), , drop = FALSE]
    write_tsv(bind_rows_fill(candidates_log), file.path(DIRS$audit, "A03_annotation_candidate_audit.tsv"))
    return(list(mapping = mapping_table, source = "scrna_full installed hgu133plus2.db"))
  }

  write_tsv(bind_rows_fill(candidates_log), file.path(DIRS$audit, "A03_annotation_candidate_audit.tsv"))
  fail("No usable full GPL570 probe annotation found. Ensure a GPL570 probe-to-gene TSV is present in /data or hgu133plus2.db is installed in scrna_full.")
}

aggregate_probe_median <- function(expr, mapping) {
  match_index <- match(rownames(expr), mapping$probe_id)
  gene <- mapping$gene[match_index]
  valid <- !is.na(gene) & nzchar(gene)
  expr_use <- expr[valid, , drop = FALSE]
  gene_use <- gene[valid]
  groups <- split(seq_along(gene_use), gene_use)
  aggregated <- vapply(
    groups,
    function(index) {
      if (length(index) == 1L) {
        as.numeric(expr_use[index, , drop = TRUE])
      } else {
        apply(expr_use[index, , drop = FALSE], 2L, stats::median, na.rm = TRUE)
      }
    },
    numeric(ncol(expr_use))
  )
  aggregated <- t(aggregated)
  colnames(aggregated) <- colnames(expr_use)
  rownames(aggregated) <- names(groups)
  aggregated
}

aggregate_probe_iqr <- function(expr, mapping) {
  match_index <- match(rownames(expr), mapping$probe_id)
  gene <- mapping$gene[match_index]
  valid <- !is.na(gene) & nzchar(gene)
  expr_use <- expr[valid, , drop = FALSE]
  gene_use <- gene[valid]
  probe_use <- rownames(expr_use)
  probe_iqr <- apply(expr_use, 1L, stats::IQR, na.rm = TRUE)
  groups <- split(seq_along(gene_use), gene_use)
  chosen <- vapply(
    groups,
    function(index) {
      index <- index[order(-probe_iqr[index], probe_use[index])]
      index[1L]
    },
    integer(1)
  )
  aggregated <- expr_use[chosen, , drop = FALSE]
  rownames(aggregated) <- names(chosen)
  list(matrix = aggregated, selected = data.frame(gene = names(chosen), probe_id = probe_use[chosen], probe_iqr = probe_iqr[chosen], stringsAsFactors = FALSE))
}

compute_policy_scores <- function(gene_expr, marker_clean, aggregation_label) {
  rank_matrix <- apply(gene_expr, 2L, rank, na.last = "keep", ties.method = "average")
  if (!is.matrix(rank_matrix)) {
    rank_matrix <- as.matrix(rank_matrix)
  }
  rownames(rank_matrix) <- rownames(gene_expr)
  colnames(rank_matrix) <- colnames(gene_expr)

  score_rows <- list()
  coverage_rows <- list()
  for (policy in c("P00", "P01", "P11")) {
    one <- marker_clean[marker_clean$policy == policy, , drop = FALSE]
    positive_all <- unique(one$gene[one$direction == "positive_marker"])
    negative_all <- unique(one$gene[one$direction == "negative_marker"])
    positive <- intersect(positive_all, rownames(rank_matrix))
    negative <- intersect(negative_all, rownames(rank_matrix))
    if (length(positive) < 25L || length(negative) < 25L) {
      fail(policy, " has insufficient marker coverage in GSE31684: positive=", length(positive), "; negative=", length(negative))
    }
    raw <- colMeans(rank_matrix[positive, , drop = FALSE], na.rm = TRUE) - colMeans(rank_matrix[negative, , drop = FALSE], na.rm = TRUE)
    score_rows[[policy]] <- data.frame(
      gsm_id = names(raw),
      policy = policy,
      aggregation = aggregation_label,
      score_raw = as.numeric(raw),
      score_z = z_score(raw),
      stringsAsFactors = FALSE
    )
    coverage_rows[[policy]] <- data.frame(
      policy = policy,
      aggregation = aggregation_label,
      positive_total = length(positive_all),
      positive_available = length(positive),
      positive_coverage = length(positive) / length(positive_all),
      negative_total = length(negative_all),
      negative_available = length(negative),
      negative_coverage = length(negative) / length(negative_all),
      missing_positive = paste(setdiff(positive_all, positive), collapse = "|"),
      missing_negative = paste(setdiff(negative_all, negative), collapse = "|"),
      stringsAsFactors = FALSE
    )
  }
  list(scores = bind_rows_fill(score_rows), coverage = bind_rows_fill(coverage_rows))
}

scores_to_wide <- function(score_long, suffix = "") {
  gsms <- unique(score_long$gsm_id)
  result <- data.frame(gsm_id = gsms, stringsAsFactors = FALSE)
  for (policy in c("P00", "P01", "P11")) {
    one <- score_long[score_long$policy == policy, c("gsm_id", "score_raw", "score_z"), drop = FALSE]
    names(one)[2:3] <- paste0(policy, c("_raw", "_z"), suffix)
    result <- merge(result, one, by = "gsm_id", all.x = TRUE, sort = FALSE)
  }
  result
}

hc3_vcov <- function(fit) {
  X <- stats::model.matrix(fit)
  residual <- stats::residuals(fit)
  leverage <- stats::hatvalues(fit)
  omega <- residual^2 / pmax((1 - leverage)^2, 1e-12)
  bread <- tryCatch(solve(crossprod(X)), error = function(e) qr.solve(crossprod(X)))
  meat <- crossprod(X * sqrt(omega))
  bread %*% meat %*% bread
}

robust_lm_contrasts <- function(data, score_col, group_col, group_levels, family_label) {
  required <- c(score_col, group_col, "age_z", "sex")
  d <- data[complete.cases(data[, required, drop = FALSE]), required, drop = FALSE]
  names(d)[names(d) == score_col] <- "score"
  names(d)[names(d) == group_col] <- "group"
  d$group <- factor(as.character(d$group), levels = group_levels)
  d$sex <- factor(as.character(d$sex), levels = c("Female", "Male"))
  d <- droplevels(d)

  group_levels_present <- levels(d$group)
  sex_levels_present <- levels(d$sex)
  base_status <- data.frame(
    family = family_label,
    score = score_col,
    group_variable = group_col,
    n = nrow(d),
    group_levels_present = paste(group_levels_present, collapse = "|"),
    sex_levels_present = paste(sex_levels_present, collapse = "|"),
    stringsAsFactors = FALSE
  )

  if (nrow(d) < 20L || nlevels(d$group) < 2L) {
    overall <- base_status
    overall$df <- NA_integer_
    overall$robust_Wald_chisq <- NA_real_
    overall$p_robust_Wald <- NA_real_
    overall$fit_status <- "SKIPPED_INSUFFICIENT_N_OR_GROUP_LEVELS"
    overall$message <- "Adjusted pathology model requires at least 20 complete cases and at least two pathology-group levels."
    return(list(contrasts = data.frame(), overall = overall))
  }

  if (nlevels(d$sex) < 2L) {
    overall <- base_status
    overall$df <- NA_integer_
    overall$robust_Wald_chisq <- NA_real_
    overall$p_robust_Wald <- NA_real_
    overall$fit_status <- "SKIPPED_SEX_ONE_LEVEL"
    overall$message <- paste0(
      "Age/sex-adjusted pathology model was not fitted because sex had fewer than two levels after complete-case filtering. Levels: ",
      paste(sex_levels_present, collapse = "|")
    )
    return(list(contrasts = data.frame(), overall = overall))
  }

  fit <- tryCatch(
    stats::lm(score ~ group + age_z + sex, data = d),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    overall <- base_status
    overall$df <- NA_integer_
    overall$robust_Wald_chisq <- NA_real_
    overall$p_robust_Wald <- NA_real_
    overall$fit_status <- "FAILED_LM"
    overall$message <- conditionMessage(fit)
    return(list(contrasts = data.frame(), overall = overall))
  }

  vcov_hc3 <- hc3_vcov(fit)
  beta <- stats::coef(fit)
  reference_sex <- levels(d$sex)[1L]
  levels_present <- levels(d$group)
  design_rows <- lapply(
    levels_present,
    function(level) {
      new_data <- data.frame(
        group = factor(level, levels = levels_present),
        age_z = 0,
        sex = factor(reference_sex, levels = levels(d$sex))
      )
      as.numeric(
        stats::model.matrix(
          stats::delete.response(stats::terms(fit)),
          new_data,
          contrasts.arg = fit$contrasts
        )
      )
    }
  )
  names(design_rows) <- levels_present

  rows <- list()
  counter <- 0L
  if (length(levels_present) >= 2L) {
    pairs <- combn(levels_present, 2L, simplify = FALSE)
    for (pair in pairs) {
      counter <- counter + 1L
      contrast <- design_rows[[pair[2L]]] - design_rows[[pair[1L]]]
      estimate <- sum(contrast * beta)
      variance <- as.numeric(t(contrast) %*% vcov_hc3 %*% contrast)
      se <- sqrt(max(variance, 0))
      z <- estimate / se
      p <- 2 * stats::pnorm(abs(z), lower.tail = FALSE)
      rows[[counter]] <- data.frame(
        family = family_label,
        score = score_col,
        group_variable = group_col,
        comparison = paste0(pair[2L], " vs ", pair[1L]),
        n = nrow(d),
        adjusted_difference_SD = estimate,
        CI_low = estimate - 1.96 * se,
        CI_high = estimate + 1.96 * se,
        p_raw = p,
        fit_status = "SUCCESS",
        sex_levels_present = paste(sex_levels_present, collapse = "|"),
        stringsAsFactors = FALSE
      )
    }
  }
  contrasts <- bind_rows_fill(rows)
  if (nrow(contrasts) > 0L) {
    contrasts$p_holm <- stats::p.adjust(contrasts$p_raw, method = "holm")
    contrasts$p_BH <- stats::p.adjust(contrasts$p_raw, method = "BH")
  }

  group_columns <- grep("^group", names(beta))
  overall <- base_status
  overall$df <- NA_integer_
  overall$robust_Wald_chisq <- NA_real_
  overall$p_robust_Wald <- NA_real_
  overall$fit_status <- "SUCCESS"
  overall$message <- ""
  if (length(group_columns) > 0L) {
    b <- beta[group_columns]
    V <- vcov_hc3[group_columns, group_columns, drop = FALSE]
    statistic <- tryCatch(
      as.numeric(t(b) %*% solve(V) %*% b),
      error = function(e) NA_real_
    )
    df <- length(b)
    overall$df <- df
    overall$robust_Wald_chisq <- statistic
    overall$p_robust_Wald <- ifelse(
      is.finite(statistic),
      stats::pchisq(statistic, df = df, lower.tail = FALSE),
      NA_real_
    )
  }
  list(contrasts = contrasts, overall = overall)
}


rank_biserial <- function(x_high, x_low) {
  x_high <- x_high[is.finite(x_high)]
  x_low <- x_low[is.finite(x_low)]
  if (length(x_high) == 0L || length(x_low) == 0L) {
    return(NA_real_)
  }
  differences <- outer(x_high, x_low, "-")
  mean(differences > 0) - mean(differences < 0)
}

pathology_analysis <- function(data, score_col, group_col, group_levels, ordered_group, family_label) {
  required <- c(score_col, group_col)
  d <- data[complete.cases(data[, required, drop = FALSE]), required, drop = FALSE]
  names(d) <- c("score", "group")
  d$group <- factor(as.character(d$group), levels = group_levels, ordered = ordered_group)
  d <- droplevels(d)
  summary <- bind_rows_fill(lapply(
    levels(d$group),
    function(level) {
      values <- d$score[d$group == level]
      data.frame(
        family = family_label,
        score = score_col,
        group_variable = group_col,
        group = level,
        n = sum(is.finite(values)),
        mean = mean(values, na.rm = TRUE),
        sd = stats::sd(values, na.rm = TRUE),
        median = stats::median(values, na.rm = TRUE),
        Q1 = as.numeric(stats::quantile(values, 0.25, na.rm = TRUE)),
        Q3 = as.numeric(stats::quantile(values, 0.75, na.rm = TRUE)),
        stringsAsFactors = FALSE
      )
    }
  ))
  kw <- tryCatch(stats::kruskal.test(score ~ group, data = d), error = function(e) NULL)
  trend <- if (ordered_group) tryCatch(stats::cor.test(d$score, as.numeric(d$group), method = "spearman", exact = FALSE), error = function(e) NULL) else NULL
  overall <- data.frame(
    family = family_label,
    score = score_col,
    group_variable = group_col,
    n = nrow(d),
    groups = nlevels(d$group),
    kruskal_chisq = if (!is.null(kw)) as.numeric(kw$statistic) else NA_real_,
    kruskal_p = if (!is.null(kw)) kw$p.value else NA_real_,
    spearman_rho = if (!is.null(trend)) as.numeric(trend$estimate) else NA_real_,
    spearman_p = if (!is.null(trend)) trend$p.value else NA_real_,
    stringsAsFactors = FALSE
  )

  pair_rows <- list()
  pairs <- if (nlevels(d$group) >= 2L) combn(levels(d$group), 2L, simplify = FALSE) else list()
  for (i in seq_along(pairs)) {
    pair <- pairs[[i]]
    low <- d$score[d$group == pair[1L]]
    high <- d$score[d$group == pair[2L]]
    test <- tryCatch(stats::wilcox.test(high, low, exact = FALSE, conf.int = TRUE), error = function(e) NULL)
    pair_rows[[i]] <- data.frame(
      family = family_label,
      score = score_col,
      group_variable = group_col,
      comparison = paste0(pair[2L], " vs ", pair[1L]),
      n_reference = sum(is.finite(low)),
      n_comparison = sum(is.finite(high)),
      median_reference = stats::median(low, na.rm = TRUE),
      median_comparison = stats::median(high, na.rm = TRUE),
      median_difference = stats::median(high, na.rm = TRUE) - stats::median(low, na.rm = TRUE),
      Hodges_Lehmann_shift = if (!is.null(test) && length(test$estimate) > 0L) as.numeric(test$estimate[1L]) else NA_real_,
      CI_low = if (!is.null(test) && length(test$conf.int) == 2L) as.numeric(test$conf.int[1L]) else NA_real_,
      CI_high = if (!is.null(test) && length(test$conf.int) == 2L) as.numeric(test$conf.int[2L]) else NA_real_,
      rank_biserial = rank_biserial(high, low),
      p_raw = if (!is.null(test)) test$p.value else NA_real_,
      stringsAsFactors = FALSE
    )
  }
  pairs_table <- bind_rows_fill(pair_rows)
  if (nrow(pairs_table) > 0L) {
    pairs_table$p_holm <- stats::p.adjust(pairs_table$p_raw, method = "holm")
    pairs_table$p_BH <- stats::p.adjust(pairs_table$p_raw, method = "BH")
  }
  list(summary = summary, overall = overall, pairs = pairs_table)
}

safe_ggsave <- function(plot, path_base, width, height) {
  pdf_device <- if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf
  ggplot2::ggsave(paste0(path_base, ".pdf"), plot = plot, width = width, height = height, units = "in", device = pdf_device)
  ggplot2::ggsave(paste0(path_base, ".png"), plot = plot, width = width, height = height, units = "in", dpi = 600, bg = "white")
}

base_theme <- function() {
  ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(color = "#E8E8E8", linewidth = 0.3),
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      plot.subtitle = ggplot2::element_text(size = 10, color = "#444444"),
      axis.title = ggplot2::element_text(face = "bold"),
      legend.title = ggplot2::element_text(face = "bold"),
      strip.text = ggplot2::element_text(face = "bold")
    )
}

plot_pathology <- function(data, score_col, group_col, group_levels, title, subtitle, path_base) {
  d <- data[complete.cases(data[, c(score_col, group_col), drop = FALSE]), c(score_col, group_col), drop = FALSE]
  names(d) <- c("score", "group")
  d$group <- factor(as.character(d$group), levels = group_levels)
  palette <- c("#4477AA", "#66CCAA", "#EEAA55", "#AA3377", "#CC6677")
  p <- ggplot2::ggplot(d, ggplot2::aes(x = group, y = score, fill = group, color = group)) +
    ggplot2::geom_violin(alpha = 0.28, trim = FALSE, linewidth = 0.35) +
    ggplot2::geom_boxplot(width = 0.22, outlier.shape = NA, alpha = 0.72, linewidth = 0.45) +
    ggplot2::geom_jitter(width = 0.12, height = 0, size = 1.45, alpha = 0.68) +
    ggplot2::stat_summary(fun = median, geom = "point", shape = 23, size = 3.1, fill = "white", color = "#222222") +
    ggplot2::scale_fill_manual(values = palette[seq_along(group_levels)], drop = FALSE) +
    ggplot2::scale_color_manual(values = palette[seq_along(group_levels)], drop = FALSE) +
    ggplot2::labs(title = title, subtitle = subtitle, x = NULL, y = "P01 projection score (z)") +
    base_theme() +
    ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = 18, hjust = 1))
  safe_ggsave(p, path_base, 7.2, 5.6)
  p
}

plot_policy_correlations <- function(score_wide, path_base) {
  policies <- c("P00_z", "P01_z", "P11_z")
  matrix_data <- as.matrix(score_wide[, policies, drop = FALSE])
  correlation <- stats::cor(matrix_data, use = "pairwise.complete.obs", method = "spearman")
  grid <- expand.grid(x = policies, y = policies, stringsAsFactors = FALSE)
  grid$rho <- mapply(function(x, y) correlation[x, y], grid$x, grid$y)
  grid$x <- factor(grid$x, levels = policies)
  grid$y <- factor(grid$y, levels = rev(policies))
  p <- ggplot2::ggplot(grid, ggplot2::aes(x = x, y = y, fill = rho)) +
    ggplot2::geom_tile(color = "white", linewidth = 1.2) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.3f", rho)), fontface = "bold", size = 5) +
    ggplot2::scale_fill_gradient2(low = "#4477AA", mid = "white", high = "#AA3377", midpoint = 0, limits = c(-1, 1)) +
    ggplot2::coord_equal() +
    ggplot2::labs(title = "GSE31684 policy-score concordance", x = NULL, y = NULL, fill = "Spearman ρ") +
    base_theme() +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
  safe_ggsave(p, path_base, 6.2, 5.5)
  correlation
}

plot_marker_coverage <- function(coverage, path_base) {
  long <- rbind(
    data.frame(policy = coverage$policy, direction = "Positive", coverage = coverage$positive_coverage),
    data.frame(policy = coverage$policy, direction = "Negative", coverage = coverage$negative_coverage)
  )
  long$policy <- factor(long$policy, levels = c("P00", "P01", "P11"))
  p <- ggplot2::ggplot(long, ggplot2::aes(x = policy, y = coverage, fill = direction)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.72), width = 0.64, color = "white") +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f%%", 100 * coverage)), position = ggplot2::position_dodge(width = 0.72), vjust = -0.35, size = 3.7) +
    ggplot2::scale_fill_manual(values = c(Positive = "#CC6677", Negative = "#4477AA")) +
    ggplot2::scale_y_continuous(limits = c(0, 1.08), expand = ggplot2::expansion(mult = c(0, 0.02))) +
    ggplot2::labs(title = "Marker coverage in GSE31684", x = NULL, y = "Coverage", fill = "Marker direction") +
    base_theme()
  safe_ggsave(p, path_base, 6.8, 5.0)
  p
}

fit_cox_model <- function(data, endpoint, score_col, covariates, model_id, cohort_label, policy_label, aggregation_label) {
  time_col <- paste0(tolower(endpoint), "_time")
  event_col <- paste0(tolower(endpoint), "_event")
  required <- c(time_col, event_col, score_col, covariates)
  d <- data[complete.cases(data[, required, drop = FALSE]), required, drop = FALSE]
  d <- droplevels(d)
  if (nrow(d) < 20L || sum(d[[event_col]] == 1L) < 8L) {
    return(list(summary = data.frame(
      endpoint = endpoint, cohort = cohort_label, policy = policy_label, aggregation = aggregation_label,
      model_id = model_id, n = nrow(d), events = sum(d[[event_col]] == 1L), score = score_col,
      HR = NA_real_, CI_low = NA_real_, CI_high = NA_real_, p_wald = NA_real_, concordance = NA_real_, AIC = NA_real_,
      fit_status = "SKIPPED_SPARSE", message = "Insufficient n or events", stringsAsFactors = FALSE
    ), terms = data.frame(), ph = data.frame(), fit = NULL, data = d))
  }
  rhs <- paste(c(score_col, covariates), collapse = " + ")
  formula <- stats::as.formula(paste0("survival::Surv(", time_col, ", ", event_col, ") ~ ", rhs))
  warnings <- character()
  fit <- tryCatch(
    withCallingHandlers(
      survival::coxph(formula, data = d, ties = "efron", x = TRUE, y = TRUE, model = TRUE),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    return(list(summary = data.frame(
      endpoint = endpoint, cohort = cohort_label, policy = policy_label, aggregation = aggregation_label,
      model_id = model_id, n = nrow(d), events = sum(d[[event_col]] == 1L), score = score_col,
      HR = NA_real_, CI_low = NA_real_, CI_high = NA_real_, p_wald = NA_real_, concordance = NA_real_, AIC = NA_real_,
      fit_status = "FAILED", message = conditionMessage(fit), stringsAsFactors = FALSE
    ), terms = data.frame(), ph = data.frame(), fit = NULL, data = d))
  }
  fit_summary <- summary(fit)
  coefficient_table <- fit_summary$coefficients
  confidence_table <- fit_summary$conf.int
  if (!score_col %in% rownames(coefficient_table)) {
    return(list(summary = data.frame(
      endpoint = endpoint, cohort = cohort_label, policy = policy_label, aggregation = aggregation_label,
      model_id = model_id, n = nrow(d), events = sum(d[[event_col]] == 1L), score = score_col,
      HR = NA_real_, CI_low = NA_real_, CI_high = NA_real_, p_wald = NA_real_, concordance = fit_summary$concordance[1L], AIC = stats::AIC(fit),
      fit_status = "FAILED_SCORE_TERM_MISSING", message = paste(warnings, collapse = " | "), stringsAsFactors = FALSE
    ), terms = data.frame(), ph = data.frame(), fit = fit, data = d))
  }
  summary_row <- data.frame(
    endpoint = endpoint,
    cohort = cohort_label,
    policy = policy_label,
    aggregation = aggregation_label,
    model_id = model_id,
    n = nrow(d),
    events = sum(d[[event_col]] == 1L),
    score = score_col,
    HR = as.numeric(confidence_table[score_col, "exp(coef)"]),
    CI_low = as.numeric(confidence_table[score_col, "lower .95"]),
    CI_high = as.numeric(confidence_table[score_col, "upper .95"]),
    p_wald = as.numeric(coefficient_table[score_col, "Pr(>|z|)"]),
    concordance = as.numeric(fit_summary$concordance[1L]),
    AIC = stats::AIC(fit),
    fit_status = ifelse(length(warnings) > 0L, "SUCCESS_WITH_WARNING", "SUCCESS"),
    message = paste(unique(warnings), collapse = " | "),
    stringsAsFactors = FALSE
  )
  term_rows <- lapply(
    rownames(coefficient_table),
    function(term) data.frame(
      endpoint = endpoint, cohort = cohort_label, policy = policy_label, aggregation = aggregation_label,
      model_id = model_id, n = nrow(d), events = sum(d[[event_col]] == 1L), term = term,
      coefficient = as.numeric(coefficient_table[term, "coef"]),
      HR = as.numeric(confidence_table[term, "exp(coef)"]),
      CI_low = as.numeric(confidence_table[term, "lower .95"]),
      CI_high = as.numeric(confidence_table[term, "upper .95"]),
      p_wald = as.numeric(coefficient_table[term, "Pr(>|z|)"]),
      stringsAsFactors = FALSE
    )
  )
  ph_table <- tryCatch(
    {
      zph <- survival::cox.zph(fit, transform = "km")
      table <- as.data.frame(zph$table)
      table$term <- rownames(table)
      rownames(table) <- NULL
      names(table)[names(table) == "p"] <- "PH_p"
      table$endpoint <- endpoint
      table$cohort <- cohort_label
      table$policy <- policy_label
      table$aggregation <- aggregation_label
      table$model_id <- model_id
      table
    },
    error = function(e) data.frame(
      endpoint = endpoint, cohort = cohort_label, policy = policy_label, aggregation = aggregation_label,
      model_id = model_id, term = "PH_TEST_FAILED", chisq = NA_real_, df = NA_real_, PH_p = NA_real_, message = conditionMessage(e),
      stringsAsFactors = FALSE
    )
  )
  list(summary = summary_row, terms = bind_rows_fill(term_rows), ph = ph_table, fit = fit, data = d)
}

bootstrap_cox <- function(data, endpoint, score_col, covariates, cohort_label, policy_label, aggregation_label, B = 1000L, seed = 20260803L) {
  time_col <- paste0(tolower(endpoint), "_time")
  event_col <- paste0(tolower(endpoint), "_event")
  required <- c(time_col, event_col, score_col, covariates)
  d <- data[complete.cases(data[, required, drop = FALSE]), required, drop = FALSE]
  d <- droplevels(d)
  rhs <- paste(c(score_col, covariates), collapse = " + ")
  formula <- stats::as.formula(paste0("survival::Surv(", time_col, ", ", event_col, ") ~ ", rhs))
  set.seed(seed)
  values <- rep(NA_real_, B)
  for (i in seq_len(B)) {
    index <- sample.int(nrow(d), nrow(d), replace = TRUE)
    boot <- droplevels(d[index, , drop = FALSE])
    if (sum(boot[[event_col]] == 1L) < 8L) {
      next
    }
    fit <- tryCatch(survival::coxph(formula, data = boot, ties = "efron"), error = function(e) NULL)
    if (!is.null(fit) && score_col %in% names(stats::coef(fit))) {
      values[i] <- exp(stats::coef(fit)[score_col])
    }
  }
  use <- values[is.finite(values)]
  data.frame(
    endpoint = endpoint,
    cohort = cohort_label,
    policy = policy_label,
    aggregation = aggregation_label,
    requested = B,
    successful = length(use),
    median_HR = ifelse(length(use) > 0L, stats::median(use), NA_real_),
    percentile_2_5 = ifelse(length(use) > 0L, as.numeric(stats::quantile(use, 0.025)), NA_real_),
    percentile_97_5 = ifelse(length(use) > 0L, as.numeric(stats::quantile(use, 0.975)), NA_real_),
    proportion_HR_above_1 = ifelse(length(use) > 0L, mean(use > 1), NA_real_),
    stringsAsFactors = FALSE
  )
}

landmark_cox <- function(data, endpoint, score_col, covariates, landmarks = c(1, 3)) {
  time_col <- paste0(tolower(endpoint), "_time")
  event_col <- paste0(tolower(endpoint), "_event")
  rows <- list()
  for (i in seq_along(landmarks)) {
    landmark <- landmarks[i]
    d <- data[is.finite(data[[time_col]]) & data[[time_col]] > landmark, , drop = FALSE]
    d$landmark_time <- d[[time_col]] - landmark
    required <- c("landmark_time", event_col, score_col, covariates)
    d <- droplevels(d[complete.cases(d[, required, drop = FALSE]), required, drop = FALSE])
    formula <- stats::as.formula(paste0("survival::Surv(landmark_time, ", event_col, ") ~ ", paste(c(score_col, covariates), collapse = " + ")))
    fit <- tryCatch(survival::coxph(formula, data = d, ties = "efron"), error = function(e) NULL)
    if (is.null(fit) || !score_col %in% names(stats::coef(fit))) {
      rows[[i]] <- data.frame(endpoint = endpoint, landmark_month = landmark, n = nrow(d), events = sum(d[[event_col]] == 1L), HR = NA_real_, CI_low = NA_real_, CI_high = NA_real_, p = NA_real_, status = "FAILED", stringsAsFactors = FALSE)
    } else {
      s <- summary(fit)
      rows[[i]] <- data.frame(
        endpoint = endpoint, landmark_month = landmark, n = nrow(d), events = sum(d[[event_col]] == 1L),
        HR = s$conf.int[score_col, "exp(coef)"], CI_low = s$conf.int[score_col, "lower .95"], CI_high = s$conf.int[score_col, "upper .95"],
        p = s$coefficients[score_col, "Pr(>|z|)"], status = "SUCCESS", stringsAsFactors = FALSE
      )
    }
  }
  bind_rows_fill(rows)
}

fit_time_varying <- function(data, endpoint, score_col) {
  time_col <- paste0(tolower(endpoint), "_time")
  event_col <- paste0(tolower(endpoint), "_event")
  required <- c(time_col, event_col, score_col, "age_z", "sex", "t_stage4_model", "node_status")
  d <- droplevels(data[complete.cases(data[, required, drop = FALSE]), required, drop = FALSE])
  formula <- stats::as.formula(paste0(
    "survival::Surv(", time_col, ", ", event_col, ") ~ ", score_col,
    " + tt(", score_col, ") + age_z + sex + strata(t_stage4_model) + node_status"
  ))
  fit <- tryCatch(
    survival::coxph(
      formula,
      data = d,
      ties = "efron",
      tt = function(x, t, ...) x * log(pmax(t, 0.1) / 12),
      x = TRUE
    ),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    return(list(table = data.frame(endpoint = endpoint, status = "FAILED", message = conditionMessage(fit), stringsAsFactors = FALSE), curve = data.frame(), fit = NULL))
  }
  beta <- stats::coef(fit)
  variance <- stats::vcov(fit)
  base_name <- score_col
  tt_name <- paste0("tt(", score_col, ")")
  if (!all(c(base_name, tt_name) %in% names(beta))) {
    return(list(table = data.frame(endpoint = endpoint, status = "FAILED_TERMS", message = paste(names(beta), collapse = "|"), stringsAsFactors = FALSE), curve = data.frame(), fit = fit))
  }
  times <- c(3, 6, 12, 24, 36, 60)
  curve_rows <- lapply(
    times,
    function(time) {
      contrast <- rep(0, length(beta))
      names(contrast) <- names(beta)
      contrast[base_name] <- 1
      contrast[tt_name] <- log(time / 12)
      estimate <- sum(contrast * beta)
      se <- sqrt(as.numeric(t(contrast) %*% variance %*% contrast))
      data.frame(endpoint = endpoint, time_month = time, HR = exp(estimate), CI_low = exp(estimate - 1.96 * se), CI_high = exp(estimate + 1.96 * se), stringsAsFactors = FALSE)
    }
  )
  curve <- bind_rows_fill(curve_rows)
  coefficient_table <- summary(fit)$coefficients
  table <- data.frame(
    endpoint = endpoint,
    n = nrow(d),
    events = sum(d[[event_col]] == 1L),
    baseline_logHR = beta[base_name],
    time_interaction = beta[tt_name],
    interaction_p = coefficient_table[tt_name, "Pr(>|z|)"],
    status = "SUCCESS",
    message = "",
    stringsAsFactors = FALSE
  )
  list(table = table, curve = curve, fit = fit)
}

fit_piecewise <- function(data, endpoint, score_col, cut_month = 12) {
  time_col <- paste0(tolower(endpoint), "_time")
  event_col <- paste0(tolower(endpoint), "_event")
  required <- c("patient_id", time_col, event_col, score_col, "age_z", "sex", "t_stage4_model", "node_status")
  d <- droplevels(data[complete.cases(data[, required, drop = FALSE]), required, drop = FALSE])

  if (nrow(d) < 20L || sum(d[[event_col]] == 1L) < 8L) {
    return(data.frame(
      endpoint = endpoint,
      status = "SKIPPED_SPARSE",
      message = "Insufficient complete cases or events",
      stringsAsFactors = FALSE
    ))
  }

  d$split_time <- d[[time_col]]
  d$split_event <- d[[event_col]]
  split_formula <- stats::as.formula(
    "Surv(split_time, split_event) ~ .",
    env = asNamespace("survival")
  )
  split_data <- tryCatch(
    survival::survSplit(
      split_formula,
      data = d,
      cut = cut_month,
      start = "tstart",
      end = "tstop",
      event = "split_event",
      episode = "period"
    ),
    error = function(e) e
  )
  if (inherits(split_data, "error")) {
    return(data.frame(endpoint = endpoint, status = "FAILED_SPLIT", message = conditionMessage(split_data), stringsAsFactors = FALSE))
  }

  split_data$score_early <- split_data[[score_col]] * as.numeric(split_data$period == 1L)
  split_data$score_late <- split_data[[score_col]] * as.numeric(split_data$period == 2L)
  formula <- stats::as.formula(
    "survival::Surv(tstart, tstop, split_event) ~ score_early + score_late + age_z + sex + strata(t_stage4_model) + node_status + cluster(patient_id)"
  )
  fit <- tryCatch(survival::coxph(formula, data = split_data, ties = "efron"), error = function(e) e)
  if (inherits(fit, "error")) {
    return(data.frame(endpoint = endpoint, status = "FAILED_MODEL", message = conditionMessage(fit), stringsAsFactors = FALSE))
  }

  s <- summary(fit)
  beta <- stats::coef(fit)
  variance <- stats::vcov(fit)
  difference <- beta["score_early"] - beta["score_late"]
  contrast <- rep(0, length(beta))
  names(contrast) <- names(beta)
  contrast["score_early"] <- 1
  contrast["score_late"] <- -1
  se_difference <- sqrt(as.numeric(t(contrast) %*% variance %*% contrast))

  data.frame(
    endpoint = endpoint,
    cut_month = cut_month,
    n = nrow(d),
    events = sum(d[[event_col]] == 1L),
    early_HR = s$conf.int["score_early", "exp(coef)"],
    early_CI_low = s$conf.int["score_early", "lower .95"],
    early_CI_high = s$conf.int["score_early", "upper .95"],
    early_p = s$coefficients["score_early", "Pr(>|z|)"],
    late_HR = s$conf.int["score_late", "exp(coef)"],
    late_CI_low = s$conf.int["score_late", "lower .95"],
    late_CI_high = s$conf.int["score_late", "upper .95"],
    late_p = s$coefficients["score_late", "Pr(>|z|)"],
    early_vs_late_logHR_difference = difference,
    early_vs_late_p = 2 * stats::pnorm(abs(difference / se_difference), lower.tail = FALSE),
    status = "SUCCESS",
    message = "",
    stringsAsFactors = FALSE
  )
}

km_rmst <- function(time, event, tau) {
  valid <- is.finite(time) & !is.na(event)
  time <- time[valid]
  event <- event[valid]
  if (length(time) < 2L) {
    return(NA_real_)
  }
  fit <- survival::survfit(survival::Surv(time, event) ~ 1)
  keep <- fit$time < tau
  event_times <- fit$time[keep]
  survival_levels <- fit$surv[keep]
  boundaries <- c(0, event_times, tau)
  levels <- c(1, survival_levels)
  sum(diff(boundaries) * levels)
}

rmst_pseudo_regression <- function(data, endpoint, score_col, tau) {
  time_col <- paste0(tolower(endpoint), "_time")
  event_col <- paste0(tolower(endpoint), "_event")
  required <- c(time_col, event_col, score_col, "age_z", "sex", "t_stage4_model", "node_status")
  d <- droplevels(data[complete.cases(data[, required, drop = FALSE]), required, drop = FALSE])
  n <- nrow(d)
  if (n < 20L) {
    return(data.frame(endpoint = endpoint, tau_month = tau, n = n, events = sum(d[[event_col]] == 1L), RMST_change_months_per_SD = NA_real_, CI_low = NA_real_, CI_high = NA_real_, p = NA_real_, status = "SKIPPED", stringsAsFactors = FALSE))
  }
  theta <- km_rmst(d[[time_col]], d[[event_col]], tau)
  leave_one_out <- vapply(
    seq_len(n),
    function(i) km_rmst(d[[time_col]][-i], d[[event_col]][-i], tau),
    numeric(1)
  )
  pseudo <- n * theta - (n - 1) * leave_one_out
  model_data <- data.frame(
    pseudo = pseudo,
    score = d[[score_col]],
    age_z = d$age_z,
    sex = d$sex,
    t_stage4_model = d$t_stage4_model,
    node_status = d$node_status
  )
  fit <- stats::lm(pseudo ~ score + age_z + sex + t_stage4_model + node_status, data = model_data)
  V <- hc3_vcov(fit)
  beta <- stats::coef(fit)
  if (!"score" %in% names(beta)) {
    return(data.frame(endpoint = endpoint, tau_month = tau, n = n, events = sum(d[[event_col]] == 1L), RMST_change_months_per_SD = NA_real_, CI_low = NA_real_, CI_high = NA_real_, p = NA_real_, status = "FAILED_TERM", stringsAsFactors = FALSE))
  }
  estimate <- beta["score"]
  se <- sqrt(V["score", "score"])
  data.frame(
    endpoint = endpoint,
    tau_month = tau,
    n = n,
    events = sum(d[[event_col]] == 1L),
    RMST_change_months_per_SD = estimate,
    CI_low = estimate - 1.96 * se,
    CI_high = estimate + 1.96 * se,
    p = 2 * stats::pnorm(abs(estimate / se), lower.tail = FALSE),
    status = "SUCCESS",
    stringsAsFactors = FALSE
  )
}

km_plot <- function(data, endpoint, group_col, title, path_base) {
  time_col <- paste0(tolower(endpoint), "_time")
  event_col <- paste0(tolower(endpoint), "_event")
  required <- c(time_col, event_col, group_col)
  d <- droplevels(data[complete.cases(data[, required, drop = FALSE]), required, drop = FALSE])
  formula <- stats::as.formula(paste0("survival::Surv(", time_col, ", ", event_col, ") ~ ", group_col))
  fit <- survival::survfit(formula, data = d)
  summary_fit <- summary(fit)
  if (length(summary_fit$time) == 0L) {
    return(data.frame())
  }
  plot_data <- data.frame(
    time = summary_fit$time,
    survival = summary_fit$surv,
    lower = summary_fit$lower,
    upper = summary_fit$upper,
    group = sub(paste0("^", group_col, "="), "", summary_fit$strata),
    stringsAsFactors = FALSE
  )
  start <- data.frame(time = 0, survival = 1, lower = 1, upper = 1, group = unique(plot_data$group), stringsAsFactors = FALSE)
  plot_data <- rbind(start, plot_data)
  logrank <- tryCatch(survival::survdiff(formula, data = d), error = function(e) NULL)
  logrank_p <- if (!is.null(logrank)) stats::pchisq(logrank$chisq, df = length(logrank$n) - 1L, lower.tail = FALSE) else NA_real_
  palette <- c("#4477AA", "#66CCAA", "#EEAA55", "#AA3377", "#CC6677")
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = time, y = survival, color = group)) +
    ggplot2::geom_step(linewidth = 0.9) +
    ggplot2::scale_color_manual(values = palette[seq_along(unique(plot_data$group))]) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(title = title, subtitle = paste0("Log-rank P = ", format_p(logrank_p), "; n = ", nrow(d)), x = "Time (months)", y = "Survival probability", color = NULL) +
    base_theme()
  safe_ggsave(p, path_base, 7.2, 5.6)

  risk_times <- c(0, 12, 24, 36, 60)
  risk_summary <- summary(fit, times = risk_times, extend = TRUE)
  risk_table <- data.frame(
    endpoint = endpoint,
    group = sub(paste0("^", group_col, "="), "", risk_summary$strata),
    time_month = risk_summary$time,
    n_risk = risk_summary$n.risk,
    n_event = risk_summary$n.event,
    stringsAsFactors = FALSE
  )
  write_tsv(risk_table, paste0(path_base, "_risk_table.tsv"))
  data.frame(endpoint = endpoint, group_variable = group_col, n = nrow(d), logrank_p = logrank_p, stringsAsFactors = FALSE)
}

plot_forest <- function(table, label_col, title, path_base, color = "#CC6677", x_label = "Hazard ratio") {
  d <- table[is.finite(table$HR) & is.finite(table$CI_low) & is.finite(table$CI_high), , drop = FALSE]
  if (nrow(d) == 0L) {
    return(NULL)
  }
  d$label <- as.character(d[[label_col]])
  d$label <- factor(d$label, levels = rev(d$label))
  p <- ggplot2::ggplot(d, ggplot2::aes(x = HR, y = label)) +
    ggplot2::geom_vline(xintercept = 1, linetype = 2, color = "#777777") +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = CI_low, xmax = CI_high), height = 0.18, color = color, linewidth = 0.65) +
    ggplot2::geom_point(size = 3.0, shape = 21, fill = color, color = "white", stroke = 0.5) +
    ggplot2::scale_x_log10() +
    ggplot2::labs(title = title, x = x_label, y = NULL) +
    base_theme()
  safe_ggsave(p, path_base, 8.0, max(4.5, 0.55 * nrow(d) + 2.0))
  p
}

plot_time_varying <- function(curve, title, path_base) {
  if (nrow(curve) == 0L) {
    return(NULL)
  }
  p <- ggplot2::ggplot(curve, ggplot2::aes(x = time_month, y = HR)) +
    ggplot2::geom_hline(yintercept = 1, linetype = 2, color = "#777777") +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = CI_low, ymax = CI_high), fill = "#CC6677", alpha = 0.22) +
    ggplot2::geom_line(color = "#AA3377", linewidth = 1.0) +
    ggplot2::geom_point(color = "#AA3377", size = 2.6) +
    ggplot2::scale_y_log10() +
    ggplot2::labs(title = title, x = "Time (months)", y = "Time-specific HR per 1 SD P01") +
    base_theme()
  safe_ggsave(p, path_base, 7.2, 5.2)
  p
}

main <- function() {
  writeLines("RUNNING", STATUS_FILE, useBytes = TRUE)
  log_msg("============================================================")
  log_msg("GSE31684 P01 external validation - first part V2.4.2 case-normalization fix")
  log_msg("Start: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  log_msg("Output: ", OUT)
  log_msg("============================================================")

  required_packages <- c("survival", "ggplot2")
  optional_annotation_packages <- c("AnnotationDbi", "hgu133plus2.db")
  all_audited_packages <- c(required_packages, optional_annotation_packages)
  package_status <- data.frame(
    package = all_audited_packages,
    role = c(
      rep("required_analysis", length(required_packages)),
      rep("GPL570_annotation_fallback", length(optional_annotation_packages))
    ),
    available = vapply(
      all_audited_packages,
      requireNamespace,
      logical(1),
      quietly = TRUE
    ),
    stringsAsFactors = FALSE
  )
  write_tsv(package_status, file.path(DIRS$audit, "A04_required_package_status.tsv"))
  required_status <- package_status[
    package_status$role == "required_analysis",
    ,
    drop = FALSE
  ]
  if (!all(required_status$available)) {
    fail(
      "Required R packages missing from scrna_full: ",
      paste(required_status$package[!required_status$available], collapse = ", "),
      ". Packages will not be installed automatically."
    )
  }

  log_msg("[STEP 1/10] Locate locked inputs in project data directory")
  series_file <- locate_data_file(
    "GSE31684_series_matrix (4).txt.gz",
    "Complete GSE31684 Series Matrix"
  )
  d01_file <- locate_data_file(
    "D01_GSE31684_all93_integrated.tsv",
    "GSE31684 D01 GSM-to-patient mapping"
  )
  f02_file <- locate_data_file(
    "F02_all93_analysis_data.tsv",
    "Final GSE31684 F02 clinical matrix"
  )
  marker_loaded <- load_marker_definitions()
  marker_clean <- validate_marker_definitions(marker_loaded$table)
  write_tsv(marker_clean, file.path(DIRS$audit, "A05_locked_policy_markers.tsv"))

  input_manifest <- data.frame(
    role = c(
      "GSE31684 GPL570 Series Matrix expression",
      "D01 GSM-patient map",
      "F02 final clinical model matrix",
      "P00/P01/P11 markers"
    ),
    path = c(series_file, d01_file, f02_file, marker_loaded$source),
    stringsAsFactors = FALSE
  )
  write_tsv(input_manifest, file.path(DIRS$audit, "A06_locked_input_manifest.tsv"))

  log_msg("[STEP 2/10] Read GSE31684 Series Matrix and locked clinical matrices")
  expression_result <- read_gse31684_series_matrix(series_file)
  expr <- expression_result$matrix
  write_tsv(
    expression_result$audit,
    file.path(DIRS$audit, "A07_GSE31684_series_matrix_audit.tsv")
  )
  writeLines(
    expression_result$metadata_lines,
    file.path(DIRS$audit, "A07_GSE31684_series_matrix_metadata.txt"),
    useBytes = TRUE
  )
  if (is.null(rownames(expr)) || is.null(colnames(expr))) {
    fail("The complete Series Matrix lacks probe row names or GSM column names.")
  }
  if (nrow(expr) < 50000L || ncol(expr) != 93L) {
    fail(
      "Unexpected GSE31684 Series Matrix dimensions (expected approximately 54,675 probes x 93 samples): actual=",
      nrow(expr),
      "×",
      ncol(expr)
    )
  }
  d01 <- safe_read_tsv(d01_file)
  f02 <- safe_read_tsv(f02_file)

  d01_gsm_col <- pick_column(d01, c("gsm_id", "gsm", "sample_id"), TRUE, "D01 gsm_id")
  d01_patient_col <- pick_column(d01, c("patient_id", "patient", "patientid"), TRUE, "D01 patient_id")
  f02_patient_col <- pick_column(f02, c("patient_id", "patient", "patientid"), TRUE, "F02 patient_id")
  d01[[d01_gsm_col]] <- toupper(clean_character(d01[[d01_gsm_col]]))
  d01[[d01_patient_col]] <- clean_character(d01[[d01_patient_col]])
  f02[[f02_patient_col]] <- clean_character(f02[[f02_patient_col]])
  if (anyDuplicated(d01[[d01_gsm_col]]) || anyDuplicated(d01[[d01_patient_col]]) || anyDuplicated(f02[[f02_patient_col]])) {
    fail("Duplicate GSM/patient_id in D01 or F02.")
  }
  if (nrow(d01) != 93L || nrow(f02) != 93L) {
    fail("D01/F02 do not contain 93 patients: D01=", nrow(d01), "; F02=", nrow(f02))
  }
  gsm_overlap <- intersect(colnames(expr), d01[[d01_gsm_col]])
  if (length(gsm_overlap) != 93L) {
    fail("Expression matrix and D01 GSM overlap is not 93 samples: ", length(gsm_overlap))
  }
  expr <- expr[, d01[[d01_gsm_col]], drop = FALSE]

  log_msg("[STEP 3/10] Load GPL570 mapping and aggregate probes")
  mapping_result <- load_probe_mapping(NULL, rownames(expr))
  mapping <- mapping_result$mapping
  mapping_overlap <- sum(rownames(expr) %in% mapping$probe_id)
  mapping_summary <- data.frame(
    mapping_source = mapping_result$source,
    mapping_rows = nrow(mapping),
    expression_probes = nrow(expr),
    matched_expression_probes = mapping_overlap,
    matched_fraction = mapping_overlap / nrow(expr),
    stringsAsFactors = FALSE
  )
  write_tsv(mapping_summary, file.path(DIRS$audit, "A08_probe_mapping_summary.tsv"))
  write_tsv(mapping, file.path(DIRS$audit, "A09_probe_to_gene_mapping_used.tsv"))
  if (mapping_overlap < 10000L) {
    fail("GPL570 mapping coverage is too low: ", mapping_overlap, " probes")
  }

  gene_median <- aggregate_probe_median(expr, mapping)
  iqr_result <- aggregate_probe_iqr(expr, mapping)
  gene_iqr <- iqr_result$matrix
  write_tsv(iqr_result$selected, file.path(SENS_DIRS$iqr_probe, "tables", "S01_IQR_highest_probe_selection.tsv"))
  saveRDS(gene_median, file.path(DIRS$projection, "GSE31684_gene_by_sample_expression_median_probe_aggregation.rds"), compress = TRUE)
  saveRDS(gene_iqr, file.path(SENS_DIRS$iqr_probe, "GSE31684_gene_by_sample_expression_IQR_probe.rds"), compress = TRUE)

  log_msg("[STEP 4/10] Compute locked P00/P01/P11 rank-projection scores")
  score_median <- compute_policy_scores(gene_median, marker_clean, "samplewise_probe_median")
  score_iqr <- compute_policy_scores(gene_iqr, marker_clean, "highest_IQR_probe")
  write_tsv(score_median$coverage, file.path(DIRS$projection, "P00_P01_P11_marker_coverage_main_median.tsv"))
  write_tsv(score_iqr$coverage, file.path(SENS_DIRS$iqr_probe, "tables", "P00_P01_P11_marker_coverage_IQR.tsv"))
  write_tsv(score_median$scores, file.path(DIRS$projection, "P00_P01_P11_patient_scores_long_main_median.tsv"))
  write_tsv(score_iqr$scores, file.path(SENS_DIRS$iqr_probe, "tables", "P00_P01_P11_patient_scores_long_IQR.tsv"))
  score_wide <- scores_to_wide(score_median$scores)
  score_iqr_wide <- scores_to_wide(score_iqr$scores, suffix = "_IQR")
  score_wide <- merge(score_wide, score_iqr_wide, by = "gsm_id", all.x = TRUE, sort = FALSE)
  score_wide$P01_quartile <- make_quartile(score_wide$P01_z)
  score_wide$P01_top10 <- make_top10(score_wide$P01_z)
  policy_correlation <- plot_policy_correlations(score_wide, file.path(DIRS$figures_main, "Figure01_policy_score_correlations"))
  correlation_table <- data.frame(policy = rownames(policy_correlation), as.data.frame(policy_correlation, check.names = FALSE), stringsAsFactors = FALSE)
  write_tsv(correlation_table, file.path(DIRS$projection, "P00_P01_P11_Spearman_correlation_matrix.tsv"))
  plot_marker_coverage(score_median$coverage, file.path(DIRS$figures_main, "Figure02_marker_coverage"))

  log_msg("[STEP 5/10] Merge scores with final 93-patient clinical matrix")
  mapping_table <- d01[, c(d01_gsm_col, d01_patient_col), drop = FALSE]
  names(mapping_table) <- c("gsm_id", "patient_id")
  clinical <- merge(f02, mapping_table, by.x = f02_patient_col, by.y = "patient_id", all.x = TRUE, sort = FALSE)
  names(clinical)[names(clinical) == f02_patient_col] <- "patient_id"
  clinical <- merge(clinical, score_wide, by = "gsm_id", all.x = TRUE, sort = FALSE)
  d01_extra <- d01
  names(d01_extra)[names(d01_extra) == d01_gsm_col] <- "gsm_id"
  names(d01_extra)[names(d01_extra) == d01_patient_col] <- "patient_id"
  clinical <- merge(clinical, d01_extra, by = c("gsm_id", "patient_id"), all.x = TRUE, sort = FALSE, suffixes = c("", "_D01"))

  age_col <- pick_column(clinical, c("age_D01", "age", "age_at_index"), TRUE, "age")
  sex_col <- pick_column(clinical, c("sex_D01", "gender_D01", "sex", "gender"), TRUE, "sex")
  t_col <- pick_column(
    clinical,
    c("pathologic_t_D01", "pathologic_t", "pathologic_T", "t_stage4"),
    TRUE,
    "pathologic T"
  )
  node_col <- pick_column(
    clinical,
    c("node_status_D01", "pathologic_n_D01", "pathologic_n", "pathologic_N", "node_status"),
    TRUE,
    "node status"
  )
  time_col <- pick_column(clinical, c("survival_month", "survival_months", "os_time", "survival_month_D01"), TRUE, "survival month")
  os_event_col <- pick_column(clinical, c("os_event", "OS_event", "overall_survival_event", "os_event_D01"), TRUE, "OS event")
  dss_event_col <- pick_column(clinical, c("dss_event", "DSS_event", "css_event", "disease_specific_event", "dss_event_D01"), TRUE, "DSS event")
  pure_tcc_col <- pick_column(clinical, c("pure_tcc_D01", "pure_tcc", "histology_pure_tcc"), FALSE, "pure TCC")
  preop_col <- pick_column(
    clinical,
    c(
      "pre_rc_chemo_D01",
      "pre_rc_chemo",
      "preop_chemo_D01",
      "preop_chemo",
      "preoperative_chemo",
      "neoadjuvant_chemo",
      "nac"
    ),
    FALSE,
    "preoperative chemotherapy"
  )

  clinical$age <- suppressWarnings(as.numeric(clinical[[age_col]]))
  clinical$age_z <- z_score(clinical$age)
  clinical$sex <- standardize_sex(clinical[[sex_col]])
  sex_raw_value <- clean_character(clinical[[sex_col]])
  sex_raw_mapping <- data.frame(
    raw_value = ifelse(is.na(sex_raw_value), "<NA>", sex_raw_value),
    standardized_sex = ifelse(is.na(clinical$sex), "<NA>", as.character(clinical$sex)),
    stringsAsFactors = FALSE
  )
  sex_raw_mapping <- stats::aggregate(
    rep(1L, nrow(sex_raw_mapping)),
    by = sex_raw_mapping,
    FUN = sum
  )
  names(sex_raw_mapping)[names(sex_raw_mapping) == "x"] <- "n"

  sex_audit <- data.frame(
    source_column = sex_col,
    level = c("Female", "Male", "Missing_or_unrecognized"),
    n = c(
      sum(clinical$sex == "Female", na.rm = TRUE),
      sum(clinical$sex == "Male", na.rm = TRUE),
      sum(is.na(clinical$sex))
    ),
    stringsAsFactors = FALSE
  )
  age_sex_source_audit <- data.frame(
    variable = c("age", "sex"),
    source_column = c(age_col, sex_col),
    nonmissing_n = c(sum(is.finite(clinical$age)), sum(!is.na(clinical$sex))),
    distinct_nonmissing_levels = c(length(unique(clinical$age[is.finite(clinical$age)])), nlevels(droplevels(clinical$sex))),
    stringsAsFactors = FALSE
  )
  write_tsv(sex_audit, file.path(DIRS$audit, "A10A_sex_source_and_level_audit.tsv"))
  write_tsv(sex_raw_mapping, file.path(DIRS$audit, "A10A2_sex_raw_to_standardized_mapping.tsv"))
  write_tsv(age_sex_source_audit, file.path(DIRS$audit, "A10B_age_sex_source_audit.tsv"))

  female_all93 <- sum(clinical$sex == "Female", na.rm = TRUE)
  male_all93 <- sum(clinical$sex == "Male", na.rm = TRUE)
  missing_all93 <- sum(is.na(clinical$sex))

  if (
    female_all93 != 25L ||
    male_all93 != 68L ||
    missing_all93 != 0L
  ) {
    fail(
      "Sex QC failed for all 93 GSE31684 patients. Expected Female=25, Male=68, Missing=0; actual Female=",
      female_all93,
      ", Male=", male_all93,
      ", Missing=", missing_all93,
      ". See A10A2_sex_raw_to_standardized_mapping.tsv."
    )
  }
  raw_t_value <- clean_character(clinical[[t_col]])
  raw_node_value <- clean_character(clinical[[node_col]])
  raw_t_node_audit <- rbind(
    data.frame(
      variable = "pathologic_T_raw",
      source_column = t_col,
      raw_value = names(table(raw_t_value, useNA = "ifany")),
      n = as.integer(table(raw_t_value, useNA = "ifany")),
      stringsAsFactors = FALSE
    ),
    data.frame(
      variable = "node_status_raw",
      source_column = node_col,
      raw_value = names(table(raw_node_value, useNA = "ifany")),
      n = as.integer(table(raw_node_value, useNA = "ifany")),
      stringsAsFactors = FALSE
    )
  )
  raw_t_node_audit$raw_value[is.na(raw_t_node_audit$raw_value) | !nzchar(raw_t_node_audit$raw_value)] <- "<NA>"
  write_tsv(
    raw_t_node_audit,
    file.path(DIRS$audit, "A10D0_raw_T_and_N_value_counts.tsv")
  )

  clinical$t_stage4 <- standardize_t_stage4(clinical[[t_col]])
  clinical$node_status <- standardize_node(clinical[[node_col]])

  t_counts <- c(
    "pTa/pT1" = sum(clinical$t_stage4 == "pTa/pT1", na.rm = TRUE),
    "pT2" = sum(clinical$t_stage4 == "pT2", na.rm = TRUE),
    "pT3" = sum(clinical$t_stage4 == "pT3", na.rm = TRUE),
    "pT4" = sum(clinical$t_stage4 == "pT4", na.rm = TRUE)
  )
  node_counts <- c(
    "Negative" = sum(clinical$node_status == "Negative", na.rm = TRUE),
    "Positive" = sum(clinical$node_status == "Positive", na.rm = TRUE),
    "Unknown" = sum(clinical$node_status == "Unknown", na.rm = TRUE)
  )

  clinical_source_audit <- data.frame(
    variable = c(
      rep("pathologic_T", 5L),
      rep("node_status", 4L)
    ),
    source_column = c(
      rep(t_col, 5L),
      rep(node_col, 4L)
    ),
    standardized_level = c(
      "pTa/pT1", "pT2", "pT3", "pT4", "Missing_or_unrecognized",
      "Negative", "Positive", "Unknown", "Missing_or_unrecognized"
    ),
    n = c(
      unname(t_counts),
      sum(is.na(clinical$t_stage4)),
      unname(node_counts),
      sum(is.na(clinical$node_status))
    ),
    stringsAsFactors = FALSE
  )
  write_tsv(
    clinical_source_audit,
    file.path(DIRS$audit, "A10D_T_and_N_source_level_audit.tsv")
  )

  if (!identical(as.integer(t_counts), c(15L, 17L, 42L, 19L))) {
    fail(
      "GSE31684 T-stage QC failed. Expected pTa/pT1=15, pT2=17, pT3=42, pT4=19; actual=",
      paste(names(t_counts), t_counts, sep = ":", collapse = ", "),
      ". See A10D_T_and_N_source_level_audit.tsv."
    )
  }
  if (!identical(as.integer(node_counts), c(49L, 28L, 16L))) {
    fail(
      "GSE31684 N-status QC failed. Expected Negative=49, Positive=28, Unknown=16; actual=",
      paste(names(node_counts), node_counts, sep = ":", collapse = ", "),
      ". See A10D_T_and_N_source_level_audit.tsv."
    )
  }

  log_msg(
    "Clinical source QC passed: T source=", t_col,
    " [15/17/42/19]; N source=", node_col,
    " [49/28/16]."
  )

  clinical$os_time <- suppressWarnings(as.numeric(clinical[[time_col]]))
  clinical$dss_time <- clinical$os_time
  clinical$os_event <- parse_binary_event(clinical[[os_event_col]])
  clinical$dss_event <- parse_binary_event(clinical[[dss_event_col]])
  clinical$t_stage4_model <- factor(as.character(clinical$t_stage4), levels = c("pTa/pT1", "pT2", "pT3", "pT4"), ordered = FALSE)
  clinical$strict_MIBC <- clinical$t_stage4 %in% c("pT2", "pT3", "pT4")

  female_strict_mibc <- sum(clinical$strict_MIBC & clinical$sex == "Female", na.rm = TRUE)
  male_strict_mibc <- sum(clinical$strict_MIBC & clinical$sex == "Male", na.rm = TRUE)
  strict_mibc_n <- sum(clinical$strict_MIBC, na.rm = TRUE)

  strict_mibc_sex_qc <- data.frame(
    cohort = "strict_pT2_pT4",
    n = strict_mibc_n,
    Female = female_strict_mibc,
    Male = male_strict_mibc,
    stringsAsFactors = FALSE
  )
  write_tsv(
    strict_mibc_sex_qc,
    file.path(DIRS$audit, "A10C_strict_MIBC_sex_QC.tsv")
  )

  if (
    strict_mibc_n != 78L ||
    female_strict_mibc != 21L ||
    male_strict_mibc != 57L
  ) {
    fail(
      "Sex QC failed for the strict pT2-pT4 cohort. Expected n=78, Female=21, Male=57; actual n=",
      strict_mibc_n,
      ", Female=", female_strict_mibc,
      ", Male=", male_strict_mibc,
      "."
    )
  }

  log_msg(
    "Sex QC passed: all93 Female=25, Male=68; strict pT2-pT4 Female=21, Male=57."
  )

  clinical$t_binary <- factor(ifelse(is.na(clinical$t_stage4), NA_character_, ifelse(clinical$t_stage4 %in% c("pT3", "pT4"), "High_T_pT3_pT4", "Low_T_pTa_pT1_pT2")), levels = c("Low_T_pTa_pT1_pT2", "High_T_pT3_pT4"))
  clinical$node_known <- clinical$node_status %in% c("Negative", "Positive")
  clinical$tn_burden <- NA_character_
  clinical$tn_burden[clinical$node_status == "Negative" & clinical$t_binary == "Low_T_pTa_pT1_pT2"] <- "G1_LowT_Negative"
  clinical$tn_burden[clinical$node_status == "Negative" & clinical$t_binary == "High_T_pT3_pT4"] <- "G2_HighT_Negative"
  clinical$tn_burden[clinical$node_status == "Positive"] <- "G3_Node_Positive"
  clinical$tn_burden <- factor(clinical$tn_burden, levels = c("G1_LowT_Negative", "G2_HighT_Negative", "G3_Node_Positive"), ordered = TRUE)
  clinical$pure_tcc_flag <- if (!is.na(pure_tcc_col)) parse_yes_no(clinical[[pure_tcc_col]]) else factor(rep(NA_character_, nrow(clinical)), levels = c("No", "Yes"))
  clinical$preop_chemo_flag <- if (!is.na(preop_col)) parse_yes_no(clinical[[preop_col]]) else factor(rep(NA_character_, nrow(clinical)), levels = c("No", "Yes"))

  treatment_histology_audit <- data.frame(
    variable = c(
      "pure_tcc_Yes", "pure_tcc_No", "pure_tcc_missing",
      "preop_chemo_No", "preop_chemo_Yes", "preop_chemo_missing"
    ),
    source_column = c(
      rep(ifelse(is.na(pure_tcc_col), "NOT_FOUND", pure_tcc_col), 3L),
      rep(ifelse(is.na(preop_col), "NOT_FOUND", preop_col), 3L)
    ),
    n = c(
      sum(clinical$pure_tcc_flag == "Yes", na.rm = TRUE),
      sum(clinical$pure_tcc_flag == "No", na.rm = TRUE),
      sum(is.na(clinical$pure_tcc_flag)),
      sum(clinical$preop_chemo_flag == "No", na.rm = TRUE),
      sum(clinical$preop_chemo_flag == "Yes", na.rm = TRUE),
      sum(is.na(clinical$preop_chemo_flag))
    ),
    stringsAsFactors = FALSE
  )
  write_tsv(
    treatment_histology_audit,
    file.path(DIRS$audit, "A10E_histology_and_preop_source_level_audit.tsv")
  )

  if (
    sum(clinical$pure_tcc_flag == "Yes", na.rm = TRUE) != 86L ||
    sum(clinical$pure_tcc_flag == "No", na.rm = TRUE) != 7L ||
    sum(is.na(clinical$pure_tcc_flag)) != 0L
  ) {
    fail("Pure TCC QC failed. Expected Yes=86, No=7, Missing=0.")
  }
  if (
    sum(clinical$preop_chemo_flag == "No", na.rm = TRUE) != 90L ||
    sum(clinical$preop_chemo_flag == "Yes", na.rm = TRUE) != 3L ||
    sum(is.na(clinical$preop_chemo_flag)) != 0L
  ) {
    fail(
      "Preoperative chemotherapy QC failed. Expected No=90, Yes=3, Missing=0; source column=",
      ifelse(is.na(preop_col), "NOT_FOUND", preop_col),
      ". See A10E_histology_and_preop_source_level_audit.tsv."
    )
  }

  log_msg(
    "Histology/preop QC passed: pure TCC 86/7; preop No/Yes 90/3; preop source=",
    preop_col,
    "."
  )

  clinical$P01_quartile <- make_quartile(clinical$P01_z)
  clinical$P01_top10 <- make_top10(clinical$P01_z)

  if (nrow(clinical) != 93L || anyDuplicated(clinical$patient_id) || anyDuplicated(clinical$gsm_id)) {
    fail("The final clinical matrix does not contain 93 unique patients/GSMs.")
  }
  if (any(!is.finite(clinical$P01_z)) || sum(is.finite(clinical$os_time) & !is.na(clinical$os_event)) != 93L) {
    fail("Unexpected missing P01 or OS values among the 93 patients.")
  }
  write_tsv(clinical, file.path(DIRS$tables, "GSE31684_FINAL_93_PATIENT_MATRIX_WITH_P00_P01_P11.tsv"))
  saveRDS(clinical, file.path(DIRS$tables, "GSE31684_FINAL_93_PATIENT_MATRIX_WITH_P00_P01_P11.rds"), compress = TRUE)

  cohort_summary <- data.frame(
    item = c("all93", "strict_pT2_pT4", "N_known", "pure_TCC_yes", "preop_chemo_no", "OS_events", "DSS_events"),
    n = c(
      nrow(clinical), sum(clinical$strict_MIBC, na.rm = TRUE), sum(clinical$node_known, na.rm = TRUE),
      sum(clinical$pure_tcc_flag == "Yes", na.rm = TRUE), sum(clinical$preop_chemo_flag == "No", na.rm = TRUE),
      sum(clinical$os_event == 1L, na.rm = TRUE), sum(clinical$dss_event == 1L, na.rm = TRUE)
    ),
    stringsAsFactors = FALSE
  )
  write_tsv(cohort_summary, file.path(DIRS$audit, "A10_final_cohort_summary.tsv"))

  log_msg("[STEP 6/10] Run P01 pathology associations and figures")
  pathology_specs <- list(
    T = list(group = "t_stage4", levels = c("pTa/pT1", "pT2", "pT3", "pT4"), ordered = TRUE, title = "P01 across pathological T categories"),
    N = list(group = "node_status", levels = c("Negative", "Positive"), ordered = FALSE, title = "P01 by lymph-node status"),
    TN = list(group = "tn_burden", levels = c("G1_LowT_Negative", "G2_HighT_Negative", "G3_Node_Positive"), ordered = TRUE, title = "P01 across combined TN burden")
  )
  pathology_all <- list(summary = list(), overall = list(), pairs = list(), adjusted = list(), adjusted_overall = list())
  for (name in names(pathology_specs)) {
    spec <- pathology_specs[[name]]
    d_use <- if (name == "N") clinical[clinical$node_known, , drop = FALSE] else if (name == "TN") clinical[!is.na(clinical$tn_burden), , drop = FALSE] else clinical
    result <- pathology_analysis(d_use, "P01_z", spec$group, spec$levels, spec$ordered, paste0("P01_", name))
    adjusted <- robust_lm_contrasts(d_use, "P01_z", spec$group, spec$levels, paste0("P01_", name))
    pathology_all$summary[[name]] <- result$summary
    pathology_all$overall[[name]] <- result$overall
    pathology_all$pairs[[name]] <- result$pairs
    pathology_all$adjusted[[name]] <- adjusted$contrasts
    pathology_all$adjusted_overall[[name]] <- adjusted$overall
    subtitle <- paste0("Kruskal-Wallis P = ", format_p(result$overall$kruskal_p[1L]), if (spec$ordered) paste0("; trend ρ = ", sprintf("%.2f", result$overall$spearman_rho[1L]), ", P = ", format_p(result$overall$spearman_p[1L])) else "")
    plot_pathology(d_use, "P01_z", spec$group, spec$levels, spec$title, subtitle, file.path(DIRS$figures_main, paste0("Figure_pathology_", name, "_P01")))
  }
  pathology_summary <- bind_rows_fill(pathology_all$summary)
  pathology_overall <- bind_rows_fill(pathology_all$overall)
  pathology_pairs <- bind_rows_fill(pathology_all$pairs)
  pathology_adjusted <- bind_rows_fill(pathology_all$adjusted)
  pathology_adjusted_overall <- bind_rows_fill(pathology_all$adjusted_overall)
  write_tsv(pathology_summary, file.path(DIRS$pathology, "P01_pathology_group_summaries.tsv"))
  write_tsv(pathology_overall, file.path(DIRS$pathology, "P01_pathology_overall_tests.tsv"))
  write_tsv(pathology_pairs, file.path(DIRS$pathology, "P01_pathology_all_pairwise_tests.tsv"))
  write_tsv(pathology_adjusted, file.path(DIRS$pathology, "P01_pathology_age_sex_adjusted_contrasts_HC3.tsv"))
  write_tsv(pathology_adjusted_overall, file.path(DIRS$pathology, "P01_pathology_age_sex_adjusted_overall_HC3.tsv"))

  quartile_counts <- as.data.frame(table(clinical$P01_quartile, clinical$t_stage4, useNA = "no"), stringsAsFactors = FALSE)
  names(quartile_counts) <- c("P01_quartile", "T_stage", "n")
  quartile_counts$within_quartile_fraction <- ave(quartile_counts$n, quartile_counts$P01_quartile, FUN = function(x) x / sum(x))
  p_quartile <- ggplot2::ggplot(quartile_counts, ggplot2::aes(x = P01_quartile, y = within_quartile_fraction, fill = T_stage)) +
    ggplot2::geom_col(color = "white", linewidth = 0.4) +
    ggplot2::scale_fill_manual(values = c("pTa/pT1" = "#4477AA", "pT2" = "#66CCAA", "pT3" = "#EEAA55", "pT4" = "#AA3377")) +
    ggplot2::labs(title = "Pathological T composition across P01 quartiles", x = "P01 quartile", y = "Within-quartile fraction", fill = "Pathological T") +
    base_theme()
  safe_ggsave(p_quartile, file.path(DIRS$figures_main, "Figure_P01_quartile_T_composition"), 7.2, 5.3)

  log_msg("[STEP 7/10] Run main OS/DSS Cox models, PH checks, bootstrap and KM")
  main_covariates <- list(
    M1 = character(),
    M2 = c("age_z", "sex"),
    M3 = c("age_z", "sex", "t_stage4_model", "node_status")
  )
  main_results <- list()
  main_terms <- list()
  main_ph <- list()
  main_fits <- list()
  counter <- 0L
  for (endpoint in c("OS", "DSS")) {
    for (model_id in names(main_covariates)) {
      counter <- counter + 1L
      result <- fit_cox_model(clinical, endpoint, "P01_z", main_covariates[[model_id]], model_id, "All93", "P01", "median_probe")
      main_results[[counter]] <- result$summary
      main_terms[[counter]] <- result$terms
      main_ph[[counter]] <- result$ph
      main_fits[[paste(endpoint, model_id, sep = "_")]] <- result$fit
    }
  }
  main_cox <- bind_rows_fill(main_results)
  main_cox$label <- paste(main_cox$endpoint, main_cox$model_id, sep = " — ")
  write_tsv(main_cox, file.path(DIRS$survival, "P01_main_OS_DSS_Cox_models.tsv"))
  write_tsv(bind_rows_fill(main_terms), file.path(DIRS$survival, "P01_main_OS_DSS_Cox_all_terms.tsv"))
  ph_table <- bind_rows_fill(main_ph)
  write_tsv(ph_table, file.path(DIRS$survival, "P01_main_OS_DSS_PH_tests.tsv"))
  plot_forest(main_cox, "label", "GSE31684 P01 continuous survival models", file.path(DIRS$figures_main, "Figure_survival_P01_main_models"), color = "#CC6677")

  bootstrap_table <- bind_rows_fill(lapply(
    c("OS", "DSS"),
    function(endpoint) bootstrap_cox(clinical, endpoint, "P01_z", main_covariates$M3, "All93", "P01", "median_probe", B = 1000L, seed = 20260803L + ifelse(endpoint == "OS", 1L, 2L))
  ))
  write_tsv(bootstrap_table, file.path(DIRS$survival, "P01_main_M3_bootstrap_1000.tsv"))

  km_rows <- list()
  for (endpoint in c("OS", "DSS")) {
    km_rows[[paste0(endpoint, "_quartile")]] <- km_plot(clinical, endpoint, "P01_quartile", paste0(endpoint, " by P01 quartile"), file.path(DIRS$figures_main, paste0("Figure_KM_", endpoint, "_P01_quartiles")))
    km_rows[[paste0(endpoint, "_top10")]] <- km_plot(clinical, endpoint, "P01_top10", paste0(endpoint, " by prespecified P01 Top10"), file.path(SENS_DIRS$top10, "figures", paste0("KM_", endpoint, "_P01_Top10")))
  }
  write_tsv(bind_rows_fill(km_rows), file.path(DIRS$survival, "P01_KM_logrank_summary.tsv"))

  log_msg("[STEP 8/10] Run policy, aggregation and cohort sensitivities")
  sensitivity_results <- list()
  sensitivity_index <- 0L

  for (policy in c("P00", "P11")) {
    score_col <- paste0(policy, "_z")
    folder <- SENS_DIRS[[policy]]
    policy_rows <- list()
    for (endpoint in c("OS", "DSS")) {
      fit <- fit_cox_model(clinical, endpoint, score_col, main_covariates$M3, "M3", "All93", policy, "median_probe")
      policy_rows[[endpoint]] <- fit$summary
      sensitivity_index <- sensitivity_index + 1L
      sensitivity_results[[sensitivity_index]] <- fit$summary
    }
    policy_table <- bind_rows_fill(policy_rows)
    write_tsv(policy_table, file.path(folder, "tables", paste0(policy, "_M3_OS_DSS.tsv")))
    for (name in names(pathology_specs)) {
      spec <- pathology_specs[[name]]
      d_use <- if (name == "N") clinical[clinical$node_known, , drop = FALSE] else if (name == "TN") clinical[!is.na(clinical$tn_burden), , drop = FALSE] else clinical
      result <- pathology_analysis(d_use, score_col, spec$group, spec$levels, spec$ordered, paste0(policy, "_", name))
      write_tsv(result$overall, file.path(folder, "tables", paste0(policy, "_pathology_", name, "_overall.tsv")))
    }
  }

  iqr_rows <- list()
  for (endpoint in c("OS", "DSS")) {
    fit <- fit_cox_model(clinical, endpoint, "P01_z_IQR", main_covariates$M3, "M3", "All93", "P01", "highest_IQR_probe")
    iqr_rows[[endpoint]] <- fit$summary
    sensitivity_index <- sensitivity_index + 1L
    sensitivity_results[[sensitivity_index]] <- fit$summary
  }
  write_tsv(bind_rows_fill(iqr_rows), file.path(SENS_DIRS$iqr_probe, "tables", "P01_IQR_probe_M3_OS_DSS.tsv"))
  iqr_correlation <- suppressWarnings(stats::cor(clinical$P01_z, clinical$P01_z_IQR, method = "spearman", use = "complete.obs"))
  write_tsv(data.frame(comparison = "P01 median-probe aggregation vs highest-IQR probe", spearman_rho = iqr_correlation, stringsAsFactors = FALSE), file.path(SENS_DIRS$iqr_probe, "tables", "P01_aggregation_concordance.tsv"))

  cohort_definitions <- list(
    strict_mibc = clinical$strict_MIBC,
    n_known = clinical$node_known,
    pure_tcc = if (sum(clinical$pure_tcc_flag == "Yes", na.rm = TRUE) >= 20L) clinical$pure_tcc_flag == "Yes" else rep(FALSE, nrow(clinical)),
    no_preop = if (sum(clinical$preop_chemo_flag == "No", na.rm = TRUE) >= 20L) clinical$preop_chemo_flag == "No" else rep(FALSE, nrow(clinical))
  )
  cohort_labels <- c(strict_mibc = "Strict pT2-pT4", n_known = "N known", pure_tcc = "Pure TCC", no_preop = "No preoperative chemotherapy")
  cohort_folders <- list(strict_mibc = SENS_DIRS$strict_mibc, n_known = SENS_DIRS$n_known, pure_tcc = SENS_DIRS$pure_tcc, no_preop = SENS_DIRS$no_preop)
  for (cohort_name in names(cohort_definitions)) {
    index <- cohort_definitions[[cohort_name]]
    index[is.na(index)] <- FALSE
    d <- clinical[index, , drop = FALSE]
    rows <- list()
    for (endpoint in c("OS", "DSS")) {
      fit <- fit_cox_model(d, endpoint, "P01_z", main_covariates$M3, "M3", cohort_labels[[cohort_name]], "P01", "median_probe")
      rows[[endpoint]] <- fit$summary
      sensitivity_index <- sensitivity_index + 1L
      sensitivity_results[[sensitivity_index]] <- fit$summary
    }
    write_tsv(bind_rows_fill(rows), file.path(cohort_folders[[cohort_name]], "tables", paste0("P01_M3_", cohort_name, "_OS_DSS.tsv")))
  }
  sensitivity_cox <- bind_rows_fill(sensitivity_results)
  sensitivity_cox$label <- paste(sensitivity_cox$endpoint, sensitivity_cox$cohort, sensitivity_cox$policy, sensitivity_cox$aggregation, sep = " — ")
  write_tsv(sensitivity_cox, file.path(DIRS$tables, "P01_policy_aggregation_cohort_sensitivity_Cox.tsv"))
  plot_forest(sensitivity_cox, "label", "GSE31684 P01 sensitivity analyses", file.path(DIRS$figures_main, "Figure_survival_sensitivity_forest"), color = "#AA3377")

  policy_forest <- sensitivity_cox[sensitivity_cox$cohort == "All93" & sensitivity_cox$aggregation == "median_probe" & sensitivity_cox$policy %in% c("P00", "P11"), , drop = FALSE]
  p01_m3 <- main_cox[main_cox$model_id == "M3", , drop = FALSE]
  policy_combined <- rbind(p01_m3[, intersect(names(p01_m3), names(policy_forest)), drop = FALSE], policy_forest[, intersect(names(p01_m3), names(policy_forest)), drop = FALSE])
  policy_combined$label <- paste(policy_combined$endpoint, policy_combined$policy, sep = " — ")
  plot_forest(policy_combined, "label", "P00/P01/P11 policy sensitivity", file.path(DIRS$figures_main, "Figure_policy_survival_forest"), color = "#AA3377")

  log_msg("[STEP 9/10] Evaluate PH trigger and run prespecified non-PH sensitivities only when triggered")
  landmark_rows <- list()
  timevarying_tables <- list()
  timevarying_curves <- list()
  piecewise_rows <- list()
  rmst_rows <- list()
  ph_trigger_rows <- list()

  for (endpoint in c("OS", "DSS")) {
    endpoint_ph <- ph_table[
      ph_table$endpoint == endpoint &
        ph_table$model_id == "M3" &
        ph_table$policy == "P01" &
        ph_table$term %in% c("P01_z", "GLOBAL"),
      ,
      drop = FALSE
    ]
    trigger <- nrow(endpoint_ph) > 0L && any(is.finite(endpoint_ph$PH_p) & endpoint_ph$PH_p < 0.05)
    trigger_terms <- if (trigger) {
      paste(endpoint_ph$term[is.finite(endpoint_ph$PH_p) & endpoint_ph$PH_p < 0.05], collapse = "|")
    } else {
      "NONE"
    }
    ph_trigger_rows[[endpoint]] <- data.frame(
      endpoint = endpoint,
      M3_PH_triggered = trigger,
      triggering_terms = trigger_terms,
      minimum_M3_PH_p = if (nrow(endpoint_ph) > 0L && any(is.finite(endpoint_ph$PH_p))) min(endpoint_ph$PH_p, na.rm = TRUE) else NA_real_,
      action = ifelse(trigger, "RUN_NONPH_SENSITIVITIES", "NOT_TRIGGERED_PH_NOT_VIOLATED"),
      stringsAsFactors = FALSE
    )

    if (trigger) {
      landmark_rows[[endpoint]] <- landmark_cox(
        clinical,
        endpoint,
        "P01_z",
        c("age_z", "sex", "t_stage4_model", "node_status"),
        c(1, 3)
      )
      timevarying <- fit_time_varying(clinical, endpoint, "P01_z")
      timevarying_tables[[endpoint]] <- timevarying$table
      timevarying_curves[[endpoint]] <- timevarying$curve
      piecewise_rows[[endpoint]] <- fit_piecewise(clinical, endpoint, "P01_z", 12)
      rmst_rows[[paste0(endpoint, "_36")]] <- rmst_pseudo_regression(clinical, endpoint, "P01_z", 36)
      rmst_rows[[paste0(endpoint, "_60")]] <- rmst_pseudo_regression(clinical, endpoint, "P01_z", 60)
      if (is.data.frame(timevarying$curve) && nrow(timevarying$curve) > 0L) {
        plot_time_varying(
          timevarying$curve,
          paste0(endpoint, " time-varying P01 hazard ratio"),
          file.path(SENS_DIRS$nonph, "figures", paste0("Time_varying_HR_", endpoint))
        )
      }
    } else {
      landmark_rows[[endpoint]] <- data.frame(
        endpoint = endpoint,
        landmark_month = c(1, 3),
        n = NA_integer_,
        events = NA_integer_,
        HR = NA_real_,
        CI_low = NA_real_,
        CI_high = NA_real_,
        p = NA_real_,
        status = "NOT_TRIGGERED_PH_NOT_VIOLATED",
        stringsAsFactors = FALSE
      )
      timevarying_tables[[endpoint]] <- data.frame(
        endpoint = endpoint,
        n = NA_integer_,
        events = NA_integer_,
        baseline_logHR = NA_real_,
        time_interaction = NA_real_,
        interaction_p = NA_real_,
        status = "NOT_TRIGGERED_PH_NOT_VIOLATED",
        message = "M3 P01/global PH tests did not cross the prespecified 0.05 trigger.",
        stringsAsFactors = FALSE
      )
      timevarying_curves[[endpoint]] <- data.frame()
      piecewise_rows[[endpoint]] <- data.frame(
        endpoint = endpoint,
        status = "NOT_TRIGGERED_PH_NOT_VIOLATED",
        message = "M3 P01/global PH tests did not cross the prespecified 0.05 trigger.",
        stringsAsFactors = FALSE
      )
      for (tau in c(36, 60)) {
        rmst_rows[[paste0(endpoint, "_", tau)]] <- data.frame(
          endpoint = endpoint,
          tau_month = tau,
          n = NA_integer_,
          events = NA_integer_,
          RMST_change_months_per_SD = NA_real_,
          CI_low = NA_real_,
          CI_high = NA_real_,
          p = NA_real_,
          status = "NOT_TRIGGERED_PH_NOT_VIOLATED",
          stringsAsFactors = FALSE
        )
      }
    }
  }

  ph_trigger_table <- bind_rows_fill(ph_trigger_rows)
  landmark_table <- bind_rows_fill(landmark_rows)
  timevarying_table <- bind_rows_fill(timevarying_tables)
  timevarying_curve <- bind_rows_fill(timevarying_curves)
  piecewise_table <- bind_rows_fill(piecewise_rows)
  rmst_table <- bind_rows_fill(rmst_rows)
  write_tsv(ph_trigger_table, file.path(SENS_DIRS$nonph, "tables", "P01_M3_PH_trigger_decision.tsv"))
  write_tsv(landmark_table, file.path(SENS_DIRS$nonph, "tables", "P01_landmark_1_3_months.tsv"))
  write_tsv(timevarying_table, file.path(SENS_DIRS$nonph, "tables", "P01_time_varying_Cox.tsv"))
  write_tsv(timevarying_curve, file.path(SENS_DIRS$nonph, "tables", "P01_time_specific_HR.tsv"))
  write_tsv(piecewise_table, file.path(SENS_DIRS$nonph, "tables", "P01_fixed_12_month_piecewise_Cox.tsv"))
  write_tsv(rmst_table, file.path(SENS_DIRS$nonph, "tables", "P01_pseudo_RMST_regression_36_60_months.tsv"))

  rmst_plot_data <- rmst_table[is.finite(rmst_table$RMST_change_months_per_SD), , drop = FALSE]
  if (nrow(rmst_plot_data) > 0L) {
    rmst_plot_data$HR <- rmst_plot_data$RMST_change_months_per_SD
    rmst_plot_data$label <- paste(rmst_plot_data$endpoint, paste0(rmst_plot_data$tau_month, " months"), sep = " — ")
    p_rmst <- ggplot2::ggplot(rmst_plot_data, ggplot2::aes(x = RMST_change_months_per_SD, y = factor(label, levels = rev(label)))) +
      ggplot2::geom_vline(xintercept = 0, linetype = 2, color = "#777777") +
      ggplot2::geom_errorbarh(ggplot2::aes(xmin = CI_low, xmax = CI_high), height = 0.18, color = "#4477AA", linewidth = 0.7) +
      ggplot2::geom_point(shape = 21, size = 3, fill = "#4477AA", color = "white") +
      ggplot2::labs(title = "Adjusted pseudo-observation RMST regression", x = "RMST change (months) per 1 SD P01", y = NULL) +
      base_theme()
    safe_ggsave(p_rmst, file.path(SENS_DIRS$nonph, "figures", "RMST_36_60_months_forest"), 7.5, 4.8)
  }

  log_msg("[STEP 10/10] Final QC, report and session information")
  ph_main_m3 <- ph_table[ph_table$model_id == "M3" & ph_table$policy == "P01", , drop = FALSE]
  qc <- data.frame(
    check = c(
      "GSE31684_series_matrix_54675x93_approximately",
      "D01_93_unique",
      "F02_93_unique",
      "expression_clinical_93_match",
      "marker_policy_120_positive_120_negative",
      "P01_complete_93",
      "OS_complete_93",
      "T_levels_15_17_42_19",
      "N_levels_49_28_16",
      "TN_burden_77_complete",
      "preop_chemo_90_no_3_yes",
      "main_pathology_tables_valid",
      "main_M3_OS_DSS_success",
      "bootstrap_1000_successful_majority",
      "sensitivity_M3_models_have_nonzero_n",
      "sensitivity_subfolders_created",
      "figures_PDF_and_PNG_created"
    ),
    passed = c(
      nrow(expr) >= 50000L && ncol(expr) == 93L,
      nrow(d01) == 93L && !anyDuplicated(d01[[d01_patient_col]]),
      nrow(f02) == 93L && !anyDuplicated(f02[[f02_patient_col]]),
      nrow(clinical) == 93L && length(gsm_overlap) == 93L,
      all(score_median$coverage$positive_total == 120L & score_median$coverage$negative_total == 120L),
      sum(is.finite(clinical$P01_z)) == 93L,
      sum(is.finite(clinical$os_time) & !is.na(clinical$os_event)) == 93L,
      identical(as.integer(table(clinical$t_stage4)), c(15L, 17L, 42L, 19L)),
      identical(as.integer(table(clinical$node_status)), c(49L, 28L, 16L)),
      sum(!is.na(clinical$tn_burden)) == 77L,
      sum(clinical$preop_chemo_flag == "No", na.rm = TRUE) == 90L && sum(clinical$preop_chemo_flag == "Yes", na.rm = TRUE) == 3L,
      nrow(pathology_overall) >= 3L && all(pathology_overall$n == c(93L, 77L, 77L)) && all(is.finite(pathology_overall$kruskal_p)),
      all(main_cox$fit_status[main_cox$model_id == "M3"] %in% c("SUCCESS", "SUCCESS_WITH_WARNING")) && all(main_cox$n[main_cox$model_id == "M3"] == 93L),
      all(bootstrap_table$requested == 1000L) && all(bootstrap_table$successful >= 800L),
      nrow(sensitivity_cox) > 0L && all(sensitivity_cox$n > 0L),
      all(dir.exists(unlist(SENS_DIRS, use.names = FALSE))),
      length(list.files(DIRS$figures_main, pattern = "\\.(pdf|png)$", ignore.case = TRUE)) >= 8L
    ),
    stringsAsFactors = FALSE
  )
  write_tsv(qc, file.path(OUT, "GSE31684_FIRST_PART_QC.tsv"))
  writeLines(capture.output(sessionInfo()), file.path(OUT, "GSE31684_FIRST_PART_sessionInfo.txt"), useBytes = TRUE)

  main_os_m3 <- main_cox[main_cox$endpoint == "OS" & main_cox$model_id == "M3", , drop = FALSE]
  main_dss_m3 <- main_cox[main_cox$endpoint == "DSS" & main_cox$model_id == "M3", , drop = FALSE]
  p01_t_overall <- pathology_overall[pathology_overall$group_variable == "t_stage4", , drop = FALSE]
  p01_n_overall <- pathology_overall[pathology_overall$group_variable == "node_status", , drop = FALSE]
  p01_tn_overall <- pathology_overall[pathology_overall$group_variable == "tn_burden", , drop = FALSE]
  extreme_adjusted <- pathology_adjusted[pathology_adjusted$comparison %in% c("pT4 vs pTa/pT1", "Positive vs Negative", "G3_Node_Positive vs G1_LowT_Negative"), , drop = FALSE]

  review <- c(
    "GSE31684 P01 external validation, part 1, V2.4 (clinical source correction) completed.",
    "",
    paste0("Run directory: ", OUT),
    paste0("GSE31684 Series Matrix: ", series_file),
    paste0("D01: ", d01_file),
    paste0("F02: ", f02_file),
    paste0("Marker source: ", marker_loaded$source),
    paste0("GPL570 mapping source: ", mapping_result$source),
    paste0("T-stage source: ", t_col, "; N-status source: ", node_col, "; preoperative chemotherapy source: ", preop_col, "."),
    "",
    "1. Cohort and projection",
    paste0("Full cohort: ", nrow(clinical), " patients; strict pT2-pT4: ", sum(clinical$strict_MIBC), "."),
    paste0("OS events: ", sum(clinical$os_event), "; DSS events: ", sum(clinical$dss_event), "."),
    paste0("P01 positive marker coverage: ", score_median$coverage$positive_available[score_median$coverage$policy == "P01"], "/120; negative marker coverage: ", score_median$coverage$negative_available[score_median$coverage$policy == "P01"], "/120."),
    paste0("Spearman rho for P01 median probe aggregation versus highest-IQR probe sensitivity analysis: ", signif(iqr_correlation, 4L), "."),
    "",
    "2. Primary pathological association results",
    paste0("Four-level T-stage Kruskal-Wallis P: ", if (nrow(p01_t_overall) > 0L) format_p(p01_t_overall$kruskal_p[1L]) else "NA", "; ordinal trend P: ", if (nrow(p01_t_overall) > 0L) format_p(p01_t_overall$spearman_p[1L]) else "NA", "."),
    paste0("N-status Kruskal-Wallis P: ", if (nrow(p01_n_overall) > 0L) format_p(p01_n_overall$kruskal_p[1L]) else "NA", "."),
    paste0("TN burden Kruskal-Wallis P: ", if (nrow(p01_tn_overall) > 0L) format_p(p01_tn_overall$kruskal_p[1L]) else "NA", "; ordinal trend P: ", if (nrow(p01_tn_overall) > 0L) format_p(p01_tn_overall$spearman_p[1L]) else "NA", "."),
    "Age- and sex-adjusted results for prespecified extreme-group comparisons:",
    if (nrow(extreme_adjusted) > 0L) paste(capture.output(print(extreme_adjusted, row.names = FALSE)), collapse = "\n") else "[No estimable results]",
    "",
    "3. Primary survival results",
    paste0("OS M3: HR=", if (nrow(main_os_m3) > 0L) signif(main_os_m3$HR[1L], 4L) else "NA", ", 95%CI ", if (nrow(main_os_m3) > 0L) paste0(signif(main_os_m3$CI_low[1L], 4L), "-", signif(main_os_m3$CI_high[1L], 4L)) else "NA", ", P=", if (nrow(main_os_m3) > 0L) format_p(main_os_m3$p_wald[1L]) else "NA", "."),
    paste0("DSS M3: HR=", if (nrow(main_dss_m3) > 0L) signif(main_dss_m3$HR[1L], 4L) else "NA", ", 95%CI ", if (nrow(main_dss_m3) > 0L) paste0(signif(main_dss_m3$CI_low[1L], 4L), "-", signif(main_dss_m3$CI_high[1L], 4L)) else "NA", ", P=", if (nrow(main_dss_m3) > 0L) format_p(main_dss_m3$p_wald[1L]) else "NA", "."),
    "Fixed M3 adjustment: age, sex, four-level T stage and three-category N status. OS is the primary cross-cohort endpoint; DSS is also reported.",
    "",
    "4. Sensitivity analysis directories",
    "01_IQR_highest_probe: highest-IQR probe per gene;",
    "02_policy_P00: exclude LAMC2 and ITGB4;",
    "03_policy_P11: allow LAMC2 and ITGB4;",
    "04_strict_MIBC_pT2_pT4: strict MIBC cohort of 78 patients;",
    "05_N_known: exclude unknown N status;",
    "06_pure_TCC: pure TCC (when the field is available and sample size is sufficient);",
    "07_no_preop_chemo: exclude preoperative chemotherapy (when the field is available and sample size is sufficient);",
    "08_nonPH_timevarying: time-varying effects, a 12-month split, 1- and 3-month landmarks, and 36- and 60-month pseudo-RMST are triggered only when the M3 P01 or global PH test has P<0.05;",
    "09_Top10_descriptive: top-10 KM analysis for descriptive purposes only.",
    "",
    "5. Methodological boundaries",
    "P01 is a single-cell reference-guided state rank projection score, not a cell proportion.",
    "The primary analysis did not select probes, markers, cutoffs or models based on pathological or survival outcomes.",
    "For genes with multiple probes, the primary analysis uses the within-sample median; the highest-IQR probe is used only for sensitivity analysis.",
    "All 93 patients form a radical cystectomy cohort and should not all be labeled MIBC; the 78 pT2-pT4 patients form the strict MIBC sensitivity cohort.",
    "",
    paste0("All final QC checks passed: ", all(qc$passed)),
    "",
    "Upload only the final ZIP generated in this run directory."
  )
  writeLines(review, REPORT_FILE, useBytes = TRUE)

  final_status <- if (all(qc$passed)) "SUCCESS_GSE31684_P01_FIRST_PART" else "COMPLETED_WITH_QC_WARNINGS_GSE31684_P01_FIRST_PART"
  writeLines(final_status, STATUS_FILE, useBytes = TRUE)
  log_msg(final_status)
  log_msg("End: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
}

tryCatch(
  main(),
  error = function(e) {
    error_text <- conditionMessage(e)
    writeLines(paste0("FAILED_GSE31684_P01_FIRST_PART: ", error_text), STATUS_FILE, useBytes = TRUE)
    calls <- capture.output(sys.calls())
    writeLines(
      c(
        "GSE31684 P01 part 1 run failed.",
        paste0("Error: ", error_text),
        "",
        "Call stack: ",
        calls,
        "",
        "No previous results were deleted or modified.",
        "No R packages were automatically installed or upgraded.",
        "Neither R nor the terminal was closed.",
        paste0("Log: ", LOG_FILE)
      ),
      REPORT_FILE,
      useBytes = TRUE
    )
    log_msg("[ERROR] ", error_text)
    stop(error_text, call. = FALSE)
  }
)




