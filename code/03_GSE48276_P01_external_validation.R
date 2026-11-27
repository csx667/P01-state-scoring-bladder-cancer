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
  "02_GSE48276_SECOND_PART"
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
  m0_only = file.path(DIRS$sensitivity, "04_M0_only"),
  strict_mibc = file.path(DIRS$sensitivity, "05_strict_pT2_pT4"),
  firth = file.path(DIRS$sensitivity, "06_Firth_optional"),
  nonph = file.path(DIRS$sensitivity, "07_nonPH_timevarying"),
  quartile = file.path(DIRS$sensitivity, "08_quartile_descriptive")
)

for (path in c(OUT, unlist(DIRS, use.names = FALSE), unlist(SENS_DIRS, use.names = FALSE))) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(path, "tables"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(path, "figures"), recursive = TRUE, showWarnings = FALSE)
}

dir.create(RESULT_ROOT, recursive = TRUE, showWarnings = FALSE)
writeLines(OUT, file.path(RESULT_ROOT, "LATEST_RUN_PATH.txt"), useBytes = TRUE)

STATUS_FILE <- file.path(OUT, "GSE48276_SECOND_PART_RUN_STATUS.txt")
REPORT_FILE <- file.path(OUT, "GSE48276_SECOND_PART_SEND_TO_CHATGPT.txt")
LOG_FILE <- file.path(OUT, "GSE48276_SECOND_PART_run_log.txt")

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

read_gse48276_series_matrix <- function(path) {
  # Dedicated GSE48276 Series Matrix reader; expected fixed QC is 29,377 probes x 116 samples.
  if (!file.exists(path)) {
    fail("GSE48276 Series Matrix not found: ", path)
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
      fail(policy, " has insufficient marker coverage in GSE48276: positive=", length(positive), "; negative=", length(negative))
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
    ggplot2::labs(title = "GSE48276 policy-score concordance", x = NULL, y = NULL, fill = "Spearman ρ") +
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
    ggplot2::labs(title = "Marker coverage in GSE48276", x = NULL, y = "Coverage", fill = "Marker direction") +
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


read_gpl14951_soft_mapping <- function(path) {
  if (!file.exists(path)) {
    fail("GSE48276 Family SOFT not found: ", path)
  }
  connection <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, open = "rt") else file(path, open = "rt")
  on.exit(try(close(connection), silent = TRUE), add = TRUE)
  lines <- readLines(connection, warn = FALSE, encoding = "UTF-8")
  begin <- which(trimws(lines) == "!platform_table_begin")
  end <- which(trimws(lines) == "!platform_table_end")
  if (length(begin) < 1L || length(end) < 1L) {
    fail("GPL14951 platform table boundaries not found in Family SOFT.")
  }
  begin <- begin[1L]
  end <- end[end > begin][1L]
  if (!is.finite(end) || end <= begin + 1L) {
    fail("Invalid Family SOFT platform table boundaries.")
  }
  table_lines <- lines[(begin + 1L):(end - 1L)]
  tc <- textConnection(table_lines)
  on.exit(try(close(tc), silent = TRUE), add = TRUE)
  annotation <- read.delim(
    tc,
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = "",
    fill = TRUE
  )
  if (nrow(annotation) < 20000L) {
    fail("Unexpected number of GPL14951 platform annotation rows: ", nrow(annotation))
  }
  mapping <- extract_probe_gene_mapping(annotation)
  audit <- data.frame(
    item = c(
      "source_file", "soft_total_lines", "platform_table_begin", "platform_table_end",
      "annotation_rows", "annotation_columns", "mapping_rows", "unique_probes", "unique_genes"
    ),
    value = c(
      normalizePath(path, winslash = "/", mustWork = TRUE),
      length(lines), begin, end, nrow(annotation), ncol(annotation), nrow(mapping),
      length(unique(mapping$probe_id)), length(unique(mapping$gene))
    ),
    stringsAsFactors = FALSE
  )
  list(mapping = mapping, annotation = annotation, audit = audit)
}

standardize_t_stage4_gse48276 <- function(x) {
  raw <- gsub("[^A-Z0-9]", "", toupper(clean_character(x)))
  result <- rep(NA_character_, length(raw))
  result[grepl("^PT1|^T1", raw)] <- "pT1"
  result[grepl("^PT2|^T2", raw)] <- "pT2"
  result[grepl("^PT3|^T3", raw)] <- "pT3"
  result[grepl("^PT4|^T4", raw)] <- "pT4"
  factor(result, levels = c("pT1", "pT2", "pT3", "pT4"), ordered = TRUE)
}

standardize_node_stage3_gse48276 <- function(x) {
  raw <- gsub("[^A-Z0-9]", "", toupper(clean_character(x)))
  result <- rep(NA_character_, length(raw))
  result[grepl("^PN0|^N0|^0$", raw)] <- "N0"
  result[grepl("^PN1|^N1|^1$", raw)] <- "N1"
  result[grepl("^PN2|^N2|^2$", raw)] <- "N2"
  factor(result, levels = c("N0", "N1", "N2"), ordered = TRUE)
}

standardize_m_stage_gse48276 <- function(x) {
  raw <- gsub("[^A-Z0-9]", "", toupper(clean_character(x)))
  result <- rep(NA_character_, length(raw))
  result[grepl("^PM0|^M0|^0$", raw)] <- "M0"
  result[grepl("^PM1|^M1|^1$", raw)] <- "M1"
  factor(result, levels = c("M0", "M1"))
}

raw_value_count <- function(x, variable, source_column) {
  raw <- clean_character(x)
  raw[is.na(raw)] <- "<NA>"
  table <- as.data.frame(table(raw, useNA = "ifany"), stringsAsFactors = FALSE)
  names(table) <- c("raw_value", "n")
  table$variable <- variable
  table$source_column <- source_column
  table[, c("variable", "source_column", "raw_value", "n"), drop = FALSE]
}

fit_firth_cox <- function(data, endpoint, score_col, covariates, model_id, cohort_label, policy_label, aggregation_label) {
  if (!requireNamespace("coxphf", quietly = TRUE)) {
    return(data.frame(
      endpoint = endpoint, cohort = cohort_label, policy = policy_label, aggregation = aggregation_label,
      model_id = model_id, n = NA_integer_, events = NA_integer_, score = score_col,
      HR = NA_real_, CI_low = NA_real_, CI_high = NA_real_, p = NA_real_,
      fit_status = "COXPHF_NOT_INSTALLED", message = "Optional package coxphf is not installed; No package was installed automatically.",
      stringsAsFactors = FALSE
    ))
  }
  time_col <- paste0(tolower(endpoint), "_time")
  event_col <- paste0(tolower(endpoint), "_event")
  required <- c(time_col, event_col, score_col, covariates)
  d <- droplevels(data[complete.cases(data[, required, drop = FALSE]), required, drop = FALSE])
  if (nrow(d) < 20L || sum(d[[event_col]] == 1L) < 8L) {
    return(data.frame(
      endpoint = endpoint, cohort = cohort_label, policy = policy_label, aggregation = aggregation_label,
      model_id = model_id, n = nrow(d), events = sum(d[[event_col]] == 1L), score = score_col,
      HR = NA_real_, CI_low = NA_real_, CI_high = NA_real_, p = NA_real_,
      fit_status = "SKIPPED_SPARSE", message = "Insufficient n or events", stringsAsFactors = FALSE
    ))
  }
  formula <- stats::as.formula(paste0(
    "survival::Surv(", time_col, ", ", event_col, ") ~ ",
    paste(c(score_col, covariates), collapse = " + ")
  ))
  fit <- tryCatch(coxphf::coxphf(formula, data = d), error = function(e) e)
  if (inherits(fit, "error")) {
    return(data.frame(
      endpoint = endpoint, cohort = cohort_label, policy = policy_label, aggregation = aggregation_label,
      model_id = model_id, n = nrow(d), events = sum(d[[event_col]] == 1L), score = score_col,
      HR = NA_real_, CI_low = NA_real_, CI_high = NA_real_, p = NA_real_,
      fit_status = "FAILED", message = conditionMessage(fit), stringsAsFactors = FALSE
    ))
  }
  coefficient_names <- names(fit$coefficients)
  if (!score_col %in% coefficient_names) {
    return(data.frame(
      endpoint = endpoint, cohort = cohort_label, policy = policy_label, aggregation = aggregation_label,
      model_id = model_id, n = nrow(d), events = sum(d[[event_col]] == 1L), score = score_col,
      HR = NA_real_, CI_low = NA_real_, CI_high = NA_real_, p = NA_real_,
      fit_status = "FAILED_SCORE_TERM_MISSING", message = paste(coefficient_names, collapse = "|"), stringsAsFactors = FALSE
    ))
  }
  index <- match(score_col, coefficient_names)
  data.frame(
    endpoint = endpoint, cohort = cohort_label, policy = policy_label, aggregation = aggregation_label,
    model_id = model_id, n = nrow(d), events = sum(d[[event_col]] == 1L), score = score_col,
    HR = exp(fit$coefficients[index]),
    CI_low = fit$ci.lower[index],
    CI_high = fit$ci.upper[index],
    p = fit$prob[index],
    fit_status = "SUCCESS", message = "", stringsAsFactors = FALSE
  )
}

likelihood_ratio_add_score <- function(data, endpoint, score_col, covariates, model_id) {
  time_col <- paste0(tolower(endpoint), "_time")
  event_col <- paste0(tolower(endpoint), "_event")
  required <- c(time_col, event_col, score_col, covariates)
  d <- droplevels(data[complete.cases(data[, required, drop = FALSE]), required, drop = FALSE])
  if (nrow(d) < 20L || sum(d[[event_col]] == 1L) < 8L) {
    return(data.frame(endpoint = endpoint, model_id = model_id, n = nrow(d), events = sum(d[[event_col]] == 1L), LR_chisq = NA_real_, df = NA_integer_, p = NA_real_, status = "SKIPPED_SPARSE", stringsAsFactors = FALSE))
  }
  base_rhs <- if (length(covariates) > 0L) paste(covariates, collapse = " + ") else "1"
  full_rhs <- paste(c(score_col, covariates), collapse = " + ")
  base_formula <- stats::as.formula(paste0("survival::Surv(", time_col, ", ", event_col, ") ~ ", base_rhs))
  full_formula <- stats::as.formula(paste0("survival::Surv(", time_col, ", ", event_col, ") ~ ", full_rhs))
  base_fit <- tryCatch(survival::coxph(base_formula, data = d, ties = "efron"), error = function(e) NULL)
  full_fit <- tryCatch(survival::coxph(full_formula, data = d, ties = "efron"), error = function(e) NULL)
  if (is.null(base_fit) || is.null(full_fit)) {
    return(data.frame(endpoint = endpoint, model_id = model_id, n = nrow(d), events = sum(d[[event_col]] == 1L), LR_chisq = NA_real_, df = NA_integer_, p = NA_real_, status = "FAILED", stringsAsFactors = FALSE))
  }
  statistic <- 2 * (as.numeric(logLik(full_fit)) - as.numeric(logLik(base_fit)))
  df <- attr(logLik(full_fit), "df") - attr(logLik(base_fit), "df")
  data.frame(
    endpoint = endpoint, model_id = model_id, n = nrow(d), events = sum(d[[event_col]] == 1L),
    LR_chisq = statistic, df = df, p = stats::pchisq(statistic, df = df, lower.tail = FALSE),
    base_AIC = stats::AIC(base_fit), full_AIC = stats::AIC(full_fit), AIC_change = stats::AIC(full_fit) - stats::AIC(base_fit),
    status = "SUCCESS", stringsAsFactors = FALSE
  )
}

fit_time_varying_demographic <- function(data, endpoint, score_col) {
  time_col <- paste0(tolower(endpoint), "_time")
  event_col <- paste0(tolower(endpoint), "_event")
  required <- c(time_col, event_col, score_col, "age_z", "sex")
  d <- droplevels(data[complete.cases(data[, required, drop = FALSE]), required, drop = FALSE])
  formula <- stats::as.formula(paste0(
    "survival::Surv(", time_col, ", ", event_col, ") ~ ", score_col,
    " + tt(", score_col, ") + age_z + sex"
  ))
  fit <- tryCatch(
    survival::coxph(
      formula, data = d, ties = "efron",
      tt = function(x, t, ...) x * log(pmax(t, 0.1) / 12), x = TRUE
    ),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    return(list(table = data.frame(endpoint = endpoint, status = "FAILED", message = conditionMessage(fit), stringsAsFactors = FALSE), curve = data.frame()))
  }
  beta <- stats::coef(fit)
  variance <- stats::vcov(fit)
  tt_name <- paste0("tt(", score_col, ")")
  if (!all(c(score_col, tt_name) %in% names(beta))) {
    return(list(table = data.frame(endpoint = endpoint, status = "FAILED_TERMS", message = paste(names(beta), collapse = "|"), stringsAsFactors = FALSE), curve = data.frame()))
  }
  times <- c(3, 6, 12, 24, 36, 60)
  curve <- bind_rows_fill(lapply(times, function(time) {
    contrast <- rep(0, length(beta)); names(contrast) <- names(beta)
    contrast[score_col] <- 1
    contrast[tt_name] <- log(time / 12)
    estimate <- sum(contrast * beta)
    se <- sqrt(as.numeric(t(contrast) %*% variance %*% contrast))
    data.frame(endpoint = endpoint, time_month = time, HR = exp(estimate), CI_low = exp(estimate - 1.96 * se), CI_high = exp(estimate + 1.96 * se), stringsAsFactors = FALSE)
  }))
  coefficient_table <- summary(fit)$coefficients
  table <- data.frame(
    endpoint = endpoint, n = nrow(d), events = sum(d[[event_col]] == 1L),
    baseline_logHR = beta[score_col], time_interaction = beta[tt_name],
    interaction_p = coefficient_table[tt_name, "Pr(>|z|)"], status = "SUCCESS", message = "",
    stringsAsFactors = FALSE
  )
  list(table = table, curve = curve)
}

plot_difference_forest <- function(table, title, path_base) {
  d <- table[is.finite(table$adjusted_difference_SD) & is.finite(table$CI_low) & is.finite(table$CI_high), , drop = FALSE]
  if (nrow(d) == 0L) return(NULL)
  d$label <- paste(d$family, d$policy, sep = " — ")
  d$label <- factor(d$label, levels = rev(d$label))
  p <- ggplot2::ggplot(d, ggplot2::aes(x = adjusted_difference_SD, y = label)) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, color = "#777777") +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = CI_low, xmax = CI_high), height = 0.18, color = "#AA3377", linewidth = 0.7) +
    ggplot2::geom_point(shape = 21, size = 3, fill = "#AA3377", color = "white") +
    ggplot2::labs(title = title, x = "Adjusted score difference (SD)", y = NULL) +
    base_theme()
  safe_ggsave(p, path_base, 8.2, max(4.5, 0.5 * nrow(d) + 2.0))
  p
}

main <- function() {
  writeLines("RUNNING", STATUS_FILE, useBytes = TRUE)
  log_msg("============================================================")
  log_msg("GSE48276 P01 external validation - second part V1.1 reader-function fix")
  log_msg("Start: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  log_msg("Output: ", OUT)
  log_msg("============================================================")

  required_packages <- c("survival", "ggplot2")
  optional_packages <- c("coxphf")
  package_status <- data.frame(
    package = c(required_packages, optional_packages),
    role = c(rep("required_analysis", length(required_packages)), rep("optional_firth_sensitivity", length(optional_packages))),
    available = vapply(c(required_packages, optional_packages), requireNamespace, logical(1), quietly = TRUE),
    stringsAsFactors = FALSE
  )
  write_tsv(package_status, file.path(DIRS$audit, "A04_required_package_status.tsv"))
  if (!all(package_status$available[package_status$role == "required_analysis"])) {
    fail("Required R packages missing from scrna_full: ", paste(package_status$package[package_status$role == "required_analysis" & !package_status$available], collapse = ", "), ". Packages will not be installed automatically.")
  }

  log_msg("[STEP 1/10] Locate locked GSE48276 inputs and marker policies")
  series_file <- locate_data_file("GSE48276_series_matrix (1).txt.gz", "Complete GSE48276 Series Matrix")
  soft_file <- locate_data_file("GSE48276_family.soft.gz", "GSE48276 Family SOFT and GPL14951 annotation")
  b01_file <- locate_data_file("B01_MDA_FFPE1_57_analysis_dataset.tsv", "Final GSE48276 MDA_FFPE1-57 clinical matrix")
  marker_loaded <- load_marker_definitions()
  marker_clean <- validate_marker_definitions(marker_loaded$table)
  write_tsv(marker_clean, file.path(DIRS$audit, "A05_locked_policy_markers.tsv"))
  input_manifest <- data.frame(
    role = c("GSE48276 Series Matrix expression", "GSE48276 Family SOFT/GPL14951 annotation", "B01 locked MDA_FFPE1-57 clinical matrix", "P00/P01/P11 markers"),
    path = c(series_file, soft_file, b01_file, marker_loaded$source), stringsAsFactors = FALSE
  )
  write_tsv(input_manifest, file.path(DIRS$audit, "A06_locked_input_manifest.tsv"))

  log_msg("[STEP 2/10] Read 29,377-probe x 116-sample Series Matrix and B01 clinical matrix")
  expression_result <- read_gse48276_series_matrix(series_file)
  expr <- expression_result$matrix
  write_tsv(expression_result$audit, file.path(DIRS$audit, "A07_GSE48276_series_matrix_audit.tsv"))
  writeLines(expression_result$metadata_lines, file.path(DIRS$audit, "A07_GSE48276_series_matrix_metadata.txt"), useBytes = TRUE)
  if (nrow(expr) != 29377L || ncol(expr) != 116L) {
    fail("GSE48276 Series Matrix dimensions failed fixed QC. Expected 29377x116; actual=", nrow(expr), "x", ncol(expr))
  }
  clinical <- safe_read_tsv(b01_file)
  if (nrow(clinical) != 57L) fail("B01 does not contain 57 patients: ", nrow(clinical))

  log_msg("[STEP 3/10] Parse GPL14951 annotation from Family SOFT and aggregate probes")
  soft_result <- read_gpl14951_soft_mapping(soft_file)
  mapping <- soft_result$mapping
  mapping_overlap <- sum(rownames(expr) %in% mapping$probe_id)
  write_tsv(soft_result$audit, file.path(DIRS$audit, "A08_GPL14951_soft_annotation_audit.tsv"))
  write_tsv(data.frame(column = names(soft_result$annotation), stringsAsFactors = FALSE), file.path(DIRS$audit, "A08A_GPL14951_annotation_columns.tsv"))
  write_tsv(mapping, file.path(DIRS$audit, "A09_probe_to_gene_mapping_used.tsv"))
  write_tsv(data.frame(mapping_rows = nrow(mapping), matched_expression_probes = mapping_overlap, matched_fraction = mapping_overlap / nrow(expr), stringsAsFactors = FALSE), file.path(DIRS$audit, "A09A_probe_mapping_summary.tsv"))
  if (mapping_overlap < 15000L) fail("GPL14951 probe mapping coverage is too low: ", mapping_overlap)
  lamc2_probes <- mapping[mapping$gene == "LAMC2", , drop = FALSE]
  write_tsv(lamc2_probes, file.path(DIRS$audit, "A09B_LAMC2_probe_mapping.tsv"))
  expected_lamc2 <- c("ILMN_1701424", "ILMN_1653824", "ILMN_1706519")
  if (!all(expected_lamc2 %in% lamc2_probes$probe_id)) {
    fail("GPL14951 annotation did not contain all three prespecified LAMC2 probes. Actual=", paste(lamc2_probes$probe_id, collapse = "|"))
  }
  gene_median <- aggregate_probe_median(expr, mapping)
  iqr_result <- aggregate_probe_iqr(expr, mapping)
  gene_iqr <- iqr_result$matrix
  write_tsv(iqr_result$selected, file.path(SENS_DIRS$iqr_probe, "tables", "S01_IQR_highest_probe_selection.tsv"))
  saveRDS(gene_median, file.path(DIRS$projection, "GSE48276_gene_by_sample_expression_median_probe_aggregation.rds"), compress = TRUE)
  saveRDS(gene_iqr, file.path(SENS_DIRS$iqr_probe, "GSE48276_gene_by_sample_expression_IQR_probe.rds"), compress = TRUE)

  log_msg("[STEP 4/10] Compute P00/P01/P11 rank projections across all 116 samples")
  score_median <- compute_policy_scores(gene_median, marker_clean, "samplewise_probe_median")
  score_iqr <- compute_policy_scores(gene_iqr, marker_clean, "highest_IQR_probe")
  write_tsv(score_median$coverage, file.path(DIRS$projection, "P00_P01_P11_marker_coverage_main_median.tsv"))
  write_tsv(score_iqr$coverage, file.path(SENS_DIRS$iqr_probe, "tables", "P00_P01_P11_marker_coverage_IQR.tsv"))
  write_tsv(score_median$scores, file.path(DIRS$projection, "P00_P01_P11_all116_scores_long_main_median.tsv"))
  write_tsv(score_iqr$scores, file.path(SENS_DIRS$iqr_probe, "tables", "P00_P01_P11_all116_scores_long_IQR.tsv"))
  score_wide_all <- scores_to_wide(score_median$scores)
  score_iqr_all <- scores_to_wide(score_iqr$scores, suffix = "_IQR")
  score_wide_all <- merge(score_wide_all, score_iqr_all, by = "gsm_id", all.x = TRUE, sort = FALSE)
  policy_correlation_all116 <- plot_policy_correlations(score_wide_all, file.path(DIRS$figures_main, "Figure01_policy_score_correlations_all116"))
  write_tsv(data.frame(policy = rownames(policy_correlation_all116), as.data.frame(policy_correlation_all116, check.names = FALSE), stringsAsFactors = FALSE), file.path(DIRS$projection, "P00_P01_P11_Spearman_correlation_all116.tsv"))
  plot_marker_coverage(score_median$coverage, file.path(DIRS$figures_main, "Figure02_marker_coverage"))

  log_msg("[STEP 5/10] Merge P00/P01/P11 scores into the locked MDA_FFPE1-57 clinical matrix")
  gsm_col <- pick_column(clinical, c("gsm_id", "gsm", "sample_id"), TRUE, "gsm_id")
  patient_col <- pick_column(clinical, c("ffpe_number", "patient_id", "sample_name"), TRUE, "FFPE/patient identifier")
  age_col <- pick_column(clinical, c("age"), TRUE, "age")
  sex_col <- pick_column(clinical, c("sex", "gender"), TRUE, "sex")
  time_col <- pick_column(clinical, c("survival_month", "survival_months", "survival_mo"), TRUE, "survival month")
  os_event_col <- pick_column(clinical, c("os_event", "OS_event"), TRUE, "OS event")
  dss_event_col <- pick_column(clinical, c("dss_event", "DSS_event"), TRUE, "DSS event")
  t_col <- pick_column(clinical, c("pathologic_T", "pathologic_t", "pT"), TRUE, "pathologic T")
  n_col <- pick_column(clinical, c("pathologic_N", "pathologic_n", "pN"), TRUE, "pathologic N")
  m_col <- pick_column(clinical, c("pathologic_M", "pathologic_m", "pM", "metastasis_group"), TRUE, "pathologic M")

  clinical[[gsm_col]] <- toupper(clean_character(clinical[[gsm_col]]))
  if (anyDuplicated(clinical[[gsm_col]])) fail("Duplicate B01 gsm_id values.")
  if (sum(clinical[[gsm_col]] %in% colnames(expr)) != 57L) fail("B01 does not match 57 samples in the 116-sample expression matrix.")
  names(clinical)[names(clinical) == gsm_col] <- "gsm_id"
  names(clinical)[names(clinical) == patient_col] <- "patient_id"

  score_merge <- score_wide_all
  for (policy in c("P00", "P01", "P11")) {
    names(score_merge)[names(score_merge) == paste0(policy, "_z")] <- paste0(policy, "_z_all116")
    names(score_merge)[names(score_merge) == paste0(policy, "_z_IQR")] <- paste0(policy, "_z_IQR_all116")
  }
  clinical <- merge(clinical, score_merge, by = "gsm_id", all.x = TRUE, sort = FALSE)
  if (nrow(clinical) != 57L || anyDuplicated(clinical$gsm_id)) fail("The merged score table does not contain 57 unique GSMs.")

  for (policy in c("P00", "P01", "P11")) {
    clinical[[paste0(policy, "_z")]] <- z_score(clinical[[paste0(policy, "_raw")]])
    clinical[[paste0(policy, "_z_IQR")]] <- z_score(clinical[[paste0(policy, "_raw_IQR")]])
  }
  clinical$age <- suppressWarnings(as.numeric(clinical[[age_col]]))
  clinical$age_z <- z_score(clinical$age)
  clinical$sex <- standardize_sex(clinical[[sex_col]])
  clinical$t_stage4 <- standardize_t_stage4_gse48276(clinical[[t_col]])
  clinical$node_stage3 <- standardize_node_stage3_gse48276(clinical[[n_col]])
  clinical$m_stage <- standardize_m_stage_gse48276(clinical[[m_col]])
  clinical$survival_month <- suppressWarnings(as.numeric(clinical[[time_col]]))
  clinical$os_time <- clinical$survival_month
  clinical$dss_time <- clinical$survival_month
  clinical$os_event <- parse_binary_event(clinical[[os_event_col]])
  clinical$dss_event <- parse_binary_event(clinical[[dss_event_col]])

  clinical$t_pathology3 <- factor(
    ifelse(clinical$t_stage4 %in% c("pT1", "pT2"), "pT1/pT2", as.character(clinical$t_stage4)),
    levels = c("pT1/pT2", "pT3", "pT4"), ordered = TRUE
  )
  clinical$t_group <- factor(
    ifelse(clinical$t_stage4 %in% c("pT1", "pT2"), "pT1/pT2", "pT3/pT4"),
    levels = c("pT1/pT2", "pT3/pT4")
  )
  clinical$node_group <- factor(
    ifelse(clinical$node_stage3 == "N0", "N0", "N_positive"),
    levels = c("N0", "N_positive")
  )
  clinical$tn_burden <- NA_character_
  clinical$tn_burden[clinical$t_group == "pT1/pT2" & clinical$node_stage3 == "N0"] <- "G1_LowT_N0"
  clinical$tn_burden[clinical$t_group == "pT3/pT4" & clinical$node_stage3 == "N0"] <- "G2_HighT_N0"
  clinical$tn_burden[clinical$node_group == "N_positive"] <- "G3_Node_Positive"
  clinical$tn_burden <- factor(clinical$tn_burden, levels = c("G1_LowT_N0", "G2_HighT_N0", "G3_Node_Positive"), ordered = TRUE)
  clinical$strict_MIBC <- clinical$t_stage4 %in% c("pT2", "pT3", "pT4")
  clinical$P01_quartile <- make_quartile(clinical$P01_z)

  raw_audit <- bind_rows_fill(list(
    raw_value_count(clinical[[sex_col]], "sex", sex_col),
    raw_value_count(clinical[[t_col]], "pathologic_T", t_col),
    raw_value_count(clinical[[n_col]], "pathologic_N", n_col),
    raw_value_count(clinical[[m_col]], "pathologic_M", m_col),
    raw_value_count(clinical[[os_event_col]], "os_event", os_event_col),
    raw_value_count(clinical[[dss_event_col]], "dss_event", dss_event_col)
  ))
  write_tsv(raw_audit, file.path(DIRS$audit, "A10_raw_clinical_value_counts.tsv"))

  sex_counts <- c(Female = sum(clinical$sex == "Female", na.rm = TRUE), Male = sum(clinical$sex == "Male", na.rm = TRUE))
  t4_counts <- c(pT1 = sum(clinical$t_stage4 == "pT1", na.rm = TRUE), pT2 = sum(clinical$t_stage4 == "pT2", na.rm = TRUE), pT3 = sum(clinical$t_stage4 == "pT3", na.rm = TRUE), pT4 = sum(clinical$t_stage4 == "pT4", na.rm = TRUE))
  n3_counts <- c(N0 = sum(clinical$node_stage3 == "N0", na.rm = TRUE), N1 = sum(clinical$node_stage3 == "N1", na.rm = TRUE), N2 = sum(clinical$node_stage3 == "N2", na.rm = TRUE))
  tn_counts <- c(G1 = sum(clinical$tn_burden == "G1_LowT_N0", na.rm = TRUE), G2 = sum(clinical$tn_burden == "G2_HighT_N0", na.rm = TRUE), G3 = sum(clinical$tn_burden == "G3_Node_Positive", na.rm = TRUE))
  m_counts <- c(M0 = sum(clinical$m_stage == "M0", na.rm = TRUE), M1 = sum(clinical$m_stage == "M1", na.rm = TRUE))
  clinical_qc <- data.frame(
    check = c("all57", "Female", "Male", "pT1", "pT2", "pT3", "pT4", "N0", "N1", "N2", "TN_G1", "TN_G2", "TN_G3", "M0", "M1", "strict_pT2_pT4", "OS_events", "DSS_events"),
    expected = c(57, 8, 49, 2, 12, 35, 8, 22, 9, 26, 8, 14, 35, 56, 1, 55, 26, 17),
    observed = c(nrow(clinical), sex_counts, t4_counts, n3_counts, tn_counts, m_counts, sum(clinical$strict_MIBC), sum(clinical$os_event == 1L), sum(clinical$dss_event == 1L)),
    stringsAsFactors = FALSE
  )
  clinical_qc$passed <- clinical_qc$expected == clinical_qc$observed
  write_tsv(clinical_qc, file.path(DIRS$audit, "A11_locked_MDA57_clinical_QC.tsv"))
  if (!all(clinical_qc$passed)) fail("MDA_FFPE1-57 fixed clinical QC failed. See A11_locked_MDA57_clinical_QC.tsv.")
  if (sum(is.finite(clinical$age)) != 57L || any(!is.finite(clinical$P01_z))) fail("Missing age or P01 values among the 57 patients.")
  log_msg("Locked MDA57 QC passed: Female/Male=8/49; pT=2/12/35/8; N=22/9/26; TN=8/14/35; M0/M1=56/1; OS/DSS events=26/17.")

  policy_correlation_57 <- stats::cor(clinical[, c("P00_z", "P01_z", "P11_z")], method = "spearman", use = "pairwise.complete.obs")
  write_tsv(data.frame(policy = rownames(policy_correlation_57), as.data.frame(policy_correlation_57, check.names = FALSE), stringsAsFactors = FALSE), file.path(DIRS$projection, "P00_P01_P11_Spearman_correlation_MDA57.tsv"))
  iqr_rho <- suppressWarnings(stats::cor(clinical$P01_z, clinical$P01_z_IQR, method = "spearman", use = "complete.obs"))
  write_tsv(data.frame(comparison = "P01 median-probe vs highest-IQR probe within MDA57", spearman_rho = iqr_rho, stringsAsFactors = FALSE), file.path(SENS_DIRS$iqr_probe, "tables", "P01_aggregation_concordance_MDA57.tsv"))
  write_tsv(clinical, file.path(DIRS$tables, "GSE48276_FINAL_MDA57_MATRIX_WITH_P00_P01_P11.tsv"))
  saveRDS(clinical, file.path(DIRS$tables, "GSE48276_FINAL_MDA57_MATRIX_WITH_P00_P01_P11.rds"), compress = TRUE)

  log_msg("[STEP 6/10] Run prespecified P01 pathology associations")
  pathology_specs <- list(
    T = list(group = "t_pathology3", levels = c("pT1/pT2", "pT3", "pT4"), ordered = TRUE, title = "P01 across pathological T categories"),
    N = list(group = "node_stage3", levels = c("N0", "N1", "N2"), ordered = TRUE, title = "P01 across pathological N categories"),
    TN = list(group = "tn_burden", levels = c("G1_LowT_N0", "G2_HighT_N0", "G3_Node_Positive"), ordered = TRUE, title = "P01 across combined TN burden")
  )
  pathology_results <- list(); adjusted_results <- list(); extreme_rows <- list(); counter <- 0L
  for (name in names(pathology_specs)) {
    spec <- pathology_specs[[name]]
    result <- pathology_analysis(clinical, "P01_z", spec$group, spec$levels, spec$ordered, paste0("P01_", name))
    adjusted <- robust_lm_contrasts(clinical, "P01_z", spec$group, spec$levels, paste0("P01_", name))
    pathology_results[[paste0(name, "_summary")]] <- result$summary
    pathology_results[[paste0(name, "_overall")]] <- result$overall
    pathology_results[[paste0(name, "_pairs")]] <- result$pairs
    adjusted_results[[paste0(name, "_contrasts")]] <- adjusted$contrasts
    adjusted_results[[paste0(name, "_overall")]] <- adjusted$overall
    subtitle <- paste0("Kruskal-Wallis P = ", format_p(result$overall$kruskal_p[1L]), "; trend rho = ", sprintf("%.2f", result$overall$spearman_rho[1L]), ", P = ", format_p(result$overall$spearman_p[1L]))
    plot_pathology(clinical, "P01_z", spec$group, spec$levels, spec$title, subtitle, file.path(DIRS$figures_main, paste0("Figure_pathology_", name, "_P01")))
  }
  pathology_summary <- bind_rows_fill(pathology_results[grepl("_summary$", names(pathology_results))])
  pathology_overall <- bind_rows_fill(pathology_results[grepl("_overall$", names(pathology_results))])
  pathology_pairs <- bind_rows_fill(pathology_results[grepl("_pairs$", names(pathology_results))])
  pathology_adjusted <- bind_rows_fill(adjusted_results[grepl("_contrasts$", names(adjusted_results))])
  pathology_adjusted_overall <- bind_rows_fill(adjusted_results[grepl("_overall$", names(adjusted_results))])
  write_tsv(pathology_summary, file.path(DIRS$pathology, "P01_pathology_group_summaries.tsv"))
  write_tsv(pathology_overall, file.path(DIRS$pathology, "P01_pathology_overall_tests.tsv"))
  write_tsv(pathology_pairs, file.path(DIRS$pathology, "P01_pathology_all_pairwise_tests.tsv"))
  write_tsv(pathology_adjusted, file.path(DIRS$pathology, "P01_pathology_age_sex_adjusted_contrasts_HC3.tsv"))
  write_tsv(pathology_adjusted_overall, file.path(DIRS$pathology, "P01_pathology_age_sex_adjusted_overall_HC3.tsv"))

  policy_pathology_extremes <- list()
  for (policy in c("P00", "P01", "P11")) {
    score_col <- paste0(policy, "_z")
    folder <- if (policy == "P00") SENS_DIRS$P00 else if (policy == "P11") SENS_DIRS$P11 else DIRS$pathology
    for (name in names(pathology_specs)) {
      spec <- pathology_specs[[name]]
      result <- pathology_analysis(clinical, score_col, spec$group, spec$levels, spec$ordered, paste0(policy, "_", name))
      adjusted <- robust_lm_contrasts(clinical, score_col, spec$group, spec$levels, paste0(policy, "_", name))
      if (policy != "P01") {
        write_tsv(result$overall, file.path(folder, "tables", paste0(policy, "_pathology_", name, "_overall.tsv")))
        write_tsv(adjusted$contrasts, file.path(folder, "tables", paste0(policy, "_pathology_", name, "_adjusted_contrasts.tsv")))
      }
      target_comparison <- c(T = "pT4 vs pT1/pT2", N = "N2 vs N0", TN = "G3_Node_Positive vs G1_LowT_N0")[[name]]
      one <- adjusted$contrasts[adjusted$contrasts$comparison == target_comparison, , drop = FALSE]
      if (nrow(one) > 0L) {
        one$policy <- policy
        one$family <- name
        policy_pathology_extremes[[length(policy_pathology_extremes) + 1L]] <- one
      }
    }
  }
  policy_pathology_extremes <- bind_rows_fill(policy_pathology_extremes)
  write_tsv(policy_pathology_extremes, file.path(DIRS$tables, "P00_P01_P11_pathology_extreme_adjusted_contrasts.tsv"))
  plot_difference_forest(policy_pathology_extremes, "P00/P01/P11 pathology sensitivity", file.path(DIRS$figures_main, "Figure_pathology_policy_sensitivity_forest"))

  quartile_counts <- as.data.frame(table(clinical$P01_quartile, clinical$t_pathology3, useNA = "no"), stringsAsFactors = FALSE)
  names(quartile_counts) <- c("P01_quartile", "T_stage", "n")
  quartile_counts$within_quartile_fraction <- ave(quartile_counts$n, quartile_counts$P01_quartile, FUN = function(x) x / sum(x))
  p_quartile <- ggplot2::ggplot(quartile_counts, ggplot2::aes(x = P01_quartile, y = within_quartile_fraction, fill = T_stage)) +
    ggplot2::geom_col(color = "white", linewidth = 0.4) +
    ggplot2::scale_fill_manual(values = c("pT1/pT2" = "#4477AA", "pT3" = "#EEAA55", "pT4" = "#AA3377")) +
    ggplot2::labs(title = "Pathological T composition across P01 quartiles", x = "P01 quartile", y = "Within-quartile fraction", fill = "Pathological T") + base_theme()
  safe_ggsave(p_quartile, file.path(DIRS$figures_main, "Figure_P01_quartile_T_composition"), 7.2, 5.3)

  log_msg("[STEP 7/10] Run demographic-first DSS/OS Cox models M1-M5")
  main_covariates <- list(
    M1 = character(),
    M2 = c("age_z"),
    M3 = c("sex"),
    M4 = c("age_z", "sex"),
    M5 = c("age_z", "sex", "t_group", "node_group")
  )
  main_results <- list(); main_terms <- list(); main_ph <- list(); index <- 0L
  for (endpoint in c("DSS", "OS")) {
    for (model_id in names(main_covariates)) {
      index <- index + 1L
      fit <- fit_cox_model(clinical, endpoint, "P01_z", main_covariates[[model_id]], model_id, "MDA57", "P01", "median_probe")
      main_results[[index]] <- fit$summary; main_terms[[index]] <- fit$terms; main_ph[[index]] <- fit$ph
    }
  }
  main_cox <- bind_rows_fill(main_results)
  main_cox$model_role <- ifelse(main_cox$model_id == "M4", "PRIMARY_DEMOGRAPHIC", ifelse(main_cox$model_id == "M5", "PATHOLOGY_SENSITIVITY", "SEQUENTIAL_SUPPORT"))
  main_cox$label <- paste(main_cox$endpoint, main_cox$model_id, sep = " — ")
  main_terms <- bind_rows_fill(main_terms)
  ph_table <- bind_rows_fill(main_ph)
  write_tsv(main_cox, file.path(DIRS$survival, "P01_DSS_OS_demographic_first_M1_M5.tsv"))
  write_tsv(main_terms, file.path(DIRS$survival, "P01_DSS_OS_demographic_first_all_terms.tsv"))
  write_tsv(ph_table, file.path(DIRS$survival, "P01_DSS_OS_PH_tests.tsv"))
  plot_forest(main_cox, "label", "GSE48276 continuous P01 survival models", file.path(DIRS$figures_main, "Figure_survival_P01_M1_M5"), color = "#CC6677")

  bootstrap_table <- bind_rows_fill(unlist(lapply(c("DSS", "OS"), function(endpoint) {
    list(
      bootstrap_cox(clinical, endpoint, "P01_z", main_covariates$M4, "MDA57", "P01", "median_probe_M4", B = 1000L, seed = 20260804L + ifelse(endpoint == "DSS", 1L, 2L)),
      bootstrap_cox(clinical, endpoint, "P01_z", main_covariates$M5, "MDA57", "P01", "median_probe_M5", B = 1000L, seed = 20260804L + ifelse(endpoint == "DSS", 11L, 12L))
    )
  }), recursive = FALSE))
  write_tsv(bootstrap_table, file.path(DIRS$survival, "P01_M4_M5_bootstrap_1000.tsv"))

  lrt_table <- bind_rows_fill(unlist(lapply(c("DSS", "OS"), function(endpoint) {
    list(
      likelihood_ratio_add_score(clinical, endpoint, "P01_z", main_covariates$M4, "M4_add_P01"),
      likelihood_ratio_add_score(clinical, endpoint, "P01_z", main_covariates$M5, "M5_add_P01")
    )
  }), recursive = FALSE))
  write_tsv(lrt_table, file.path(DIRS$survival, "P01_incremental_likelihood_ratio_tests.tsv"))

  firth_table <- bind_rows_fill(unlist(lapply(c("DSS", "OS"), function(endpoint) {
    list(
      fit_firth_cox(clinical, endpoint, "P01_z", main_covariates$M4, "M4_Firth", "MDA57", "P01", "median_probe"),
      fit_firth_cox(clinical, endpoint, "P01_z", main_covariates$M5, "M5_Firth", "MDA57", "P01", "median_probe")
    )
  }), recursive = FALSE))
  write_tsv(firth_table, file.path(SENS_DIRS$firth, "tables", "P01_M4_M5_Firth_optional.tsv"))

  km_rows <- list()
  for (endpoint in c("DSS", "OS")) {
    km_rows[[endpoint]] <- km_plot(clinical, endpoint, "P01_quartile", paste0(endpoint, " by P01 quartile (descriptive)"), file.path(SENS_DIRS$quartile, "figures", paste0("KM_", endpoint, "_P01_quartiles")))
  }
  write_tsv(bind_rows_fill(km_rows), file.path(SENS_DIRS$quartile, "tables", "P01_quartile_KM_logrank.tsv"))

  log_msg("[STEP 8/10] Run policy, aggregation and fixed-cohort sensitivities")
  sensitivity_results <- list(); sensitivity_index <- 0L
  add_sensitivity <- function(fit) {
    sensitivity_index <<- sensitivity_index + 1L
    sensitivity_results[[sensitivity_index]] <<- fit$summary
  }
  for (policy in c("P00", "P11")) {
    score_col <- paste0(policy, "_z")
    folder <- SENS_DIRS[[policy]]
    rows <- list()
    for (endpoint in c("DSS", "OS")) {
      for (model_id in c("M4", "M5")) {
        fit <- fit_cox_model(clinical, endpoint, score_col, main_covariates[[model_id]], model_id, "MDA57", policy, "median_probe")
        rows[[paste(endpoint, model_id, sep = "_")]] <- fit$summary; add_sensitivity(fit)
      }
    }
    write_tsv(bind_rows_fill(rows), file.path(folder, "tables", paste0(policy, "_M4_M5_DSS_OS.tsv")))
  }
  iqr_rows <- list()
  for (endpoint in c("DSS", "OS")) for (model_id in c("M4", "M5")) {
    fit <- fit_cox_model(clinical, endpoint, "P01_z_IQR", main_covariates[[model_id]], model_id, "MDA57", "P01", "highest_IQR_probe")
    iqr_rows[[paste(endpoint, model_id, sep = "_")]] <- fit$summary; add_sensitivity(fit)
  }
  write_tsv(bind_rows_fill(iqr_rows), file.path(SENS_DIRS$iqr_probe, "tables", "P01_IQR_M4_M5_DSS_OS.tsv"))

  cohort_definitions <- list(
    M0_only = clinical$m_stage == "M0",
    strict_pT2_pT4 = clinical$strict_MIBC
  )
  cohort_labels <- c(M0_only = "M0 only", strict_pT2_pT4 = "Strict pT2-pT4")
  cohort_folders <- list(M0_only = SENS_DIRS$m0_only, strict_pT2_pT4 = SENS_DIRS$strict_mibc)
  for (cohort_name in names(cohort_definitions)) {
    use <- cohort_definitions[[cohort_name]]; use[is.na(use)] <- FALSE
    d <- clinical[use, , drop = FALSE]
    rows <- list()
    for (endpoint in c("DSS", "OS")) for (model_id in c("M4", "M5")) {
      fit <- fit_cox_model(d, endpoint, "P01_z", main_covariates[[model_id]], model_id, cohort_labels[[cohort_name]], "P01", "median_probe")
      rows[[paste(endpoint, model_id, sep = "_")]] <- fit$summary; add_sensitivity(fit)
    }
    write_tsv(bind_rows_fill(rows), file.path(cohort_folders[[cohort_name]], "tables", paste0("P01_", cohort_name, "_M4_M5_DSS_OS.tsv")))
  }
  sensitivity_cox <- bind_rows_fill(sensitivity_results)
  sensitivity_cox$label <- paste(sensitivity_cox$endpoint, sensitivity_cox$model_id, sensitivity_cox$cohort, sensitivity_cox$policy, sensitivity_cox$aggregation, sep = " — ")
  write_tsv(sensitivity_cox, file.path(DIRS$tables, "P01_policy_aggregation_cohort_sensitivity_Cox.tsv"))
  plot_forest(sensitivity_cox, "label", "GSE48276 P01 sensitivity analyses", file.path(DIRS$figures_main, "Figure_survival_sensitivity_forest"), color = "#AA3377")

  policy_main <- sensitivity_cox[sensitivity_cox$cohort == "MDA57" & sensitivity_cox$aggregation == "median_probe" & sensitivity_cox$model_id == "M4", , drop = FALSE]
  p01_m4 <- main_cox[main_cox$model_id == "M4", , drop = FALSE]
  common <- intersect(names(p01_m4), names(policy_main))
  policy_forest <- rbind(p01_m4[, common, drop = FALSE], policy_main[, common, drop = FALSE])
  policy_forest$label <- paste(policy_forest$endpoint, policy_forest$policy, sep = " — ")
  plot_forest(policy_forest, "label", "P00/P01/P11 demographic-model sensitivity", file.path(DIRS$figures_main, "Figure_policy_survival_forest_M4"), color = "#AA3377")

  log_msg("[STEP 9/10] Evaluate M4 PH trigger and run time-varying P01 only if prespecified trigger is met")
  trigger_rows <- list(); tv_tables <- list(); tv_curves <- list()
  for (endpoint in c("DSS", "OS")) {
    endpoint_ph <- ph_table[ph_table$endpoint == endpoint & ph_table$model_id == "M4" & ph_table$policy == "P01" & ph_table$term %in% c("P01_z", "GLOBAL"), , drop = FALSE]
    trigger <- nrow(endpoint_ph) > 0L && any(is.finite(endpoint_ph$PH_p) & endpoint_ph$PH_p < 0.05)
    trigger_rows[[endpoint]] <- data.frame(
      endpoint = endpoint, M4_PH_triggered = trigger,
      triggering_terms = ifelse(trigger, paste(endpoint_ph$term[is.finite(endpoint_ph$PH_p) & endpoint_ph$PH_p < 0.05], collapse = "|"), "NONE"),
      minimum_M4_PH_p = if (nrow(endpoint_ph) > 0L && any(is.finite(endpoint_ph$PH_p))) min(endpoint_ph$PH_p, na.rm = TRUE) else NA_real_,
      action = ifelse(trigger, "RUN_TIME_VARYING_P01", "NOT_TRIGGERED_PH_NOT_VIOLATED"), stringsAsFactors = FALSE
    )
    if (trigger) {
      tv <- fit_time_varying_demographic(clinical, endpoint, "P01_z")
      tv_tables[[endpoint]] <- tv$table; tv_curves[[endpoint]] <- tv$curve
      if (nrow(tv$curve) > 0L) plot_time_varying(tv$curve, paste0(endpoint, " time-varying P01 HR"), file.path(SENS_DIRS$nonph, "figures", paste0("Time_varying_HR_", endpoint)))
    } else {
      tv_tables[[endpoint]] <- data.frame(endpoint = endpoint, status = "NOT_TRIGGERED_PH_NOT_VIOLATED", message = "M4 P01/global PH tests did not cross 0.05.", stringsAsFactors = FALSE)
      tv_curves[[endpoint]] <- data.frame()
    }
  }
  write_tsv(bind_rows_fill(trigger_rows), file.path(SENS_DIRS$nonph, "tables", "P01_M4_PH_trigger_decision.tsv"))
  write_tsv(bind_rows_fill(tv_tables), file.path(SENS_DIRS$nonph, "tables", "P01_time_varying_Cox.tsv"))
  write_tsv(bind_rows_fill(tv_curves), file.path(SENS_DIRS$nonph, "tables", "P01_time_specific_HR.tsv"))

  log_msg("[STEP 10/10] Final QC, Chinese report and session information")
  main_m4m5 <- main_cox[main_cox$model_id %in% c("M4", "M5"), , drop = FALSE]
  pdf_count <- length(list.files(DIRS$figures_main, pattern = "\\.pdf$", ignore.case = TRUE))
  png_count <- length(list.files(DIRS$figures_main, pattern = "\\.png$", ignore.case = TRUE))
  qc <- data.frame(
    check = c(
      "series_matrix_29377x116", "B01_57_unique_GSM", "expression_clinical_57_match",
      "marker_policy_120_positive_120_negative", "P01_complete_57", "sex_8_49",
      "T_2_12_35_8", "N_22_9_26", "TN_8_14_35", "M0_M1_56_1",
      "OS_DSS_events_26_17", "pathology_three_families_n57", "M4_M5_DSS_OS_success_n57",
      "bootstrap_1000_majority", "sensitivity_models_nonzero_n", "all_sensitivity_folders_created",
      "main_figures_pdf_png_created"
    ),
    passed = c(
      nrow(expr) == 29377L && ncol(expr) == 116L,
      nrow(clinical) == 57L && !anyDuplicated(clinical$gsm_id),
      sum(clinical$gsm_id %in% colnames(expr)) == 57L,
      all(score_median$coverage$positive_total == 120L & score_median$coverage$negative_total == 120L),
      sum(is.finite(clinical$P01_z)) == 57L,
      identical(as.integer(sex_counts), c(8L, 49L)),
      identical(as.integer(t4_counts), c(2L, 12L, 35L, 8L)),
      identical(as.integer(n3_counts), c(22L, 9L, 26L)),
      identical(as.integer(tn_counts), c(8L, 14L, 35L)),
      identical(as.integer(m_counts), c(56L, 1L)),
      sum(clinical$os_event == 1L) == 26L && sum(clinical$dss_event == 1L) == 17L,
      nrow(pathology_overall) == 3L && all(pathology_overall$n == 57L),
      nrow(main_m4m5) == 4L && all(main_m4m5$fit_status %in% c("SUCCESS", "SUCCESS_WITH_WARNING")) && all(main_m4m5$n == 57L),
      nrow(bootstrap_table) == 4L && all(bootstrap_table$requested == 1000L) && all(bootstrap_table$successful >= 800L),
      nrow(sensitivity_cox) > 0L && all(sensitivity_cox$n > 0L),
      all(dir.exists(unlist(SENS_DIRS, use.names = FALSE))),
      pdf_count >= 9L && png_count >= 9L
    ),
    stringsAsFactors = FALSE
  )
  write_tsv(qc, file.path(OUT, "GSE48276_SECOND_PART_QC.tsv"))
  writeLines(capture.output(sessionInfo()), file.path(OUT, "GSE48276_SECOND_PART_sessionInfo.txt"), useBytes = TRUE)

  get_row <- function(endpoint, model_id) main_cox[main_cox$endpoint == endpoint & main_cox$model_id == model_id, , drop = FALSE]
  dss_m4 <- get_row("DSS", "M4"); dss_m5 <- get_row("DSS", "M5")
  os_m4 <- get_row("OS", "M4"); os_m5 <- get_row("OS", "M5")
  get_path <- function(group) pathology_overall[pathology_overall$group_variable == group, , drop = FALSE]
  p_t <- get_path("t_pathology3"); p_n <- get_path("node_stage3"); p_tn <- get_path("tn_burden")
  report <- c(
    "GSE48276 P01 external validation, part 2, V1.0 completed.",
    "",
    paste0("Run directory: ", OUT),
    paste0("Series Matrix: ", series_file),
    paste0("Family SOFT: ", soft_file),
    paste0("B01 clinical matrix: ", b01_file),
    paste0("Marker source: ", marker_loaded$source),
    "",
    "1. Fixed cohort and projection",
    "Analysis was restricted to the original MDA_FFPE1-57 prognostic validation cohort; PHIL43 and the 16 NAC extension cases were not included.",
    "Fixed QC: female=8, male=49; pT1/pT2/pT3/pT4=2/12/35/8; N0/N1/N2=22/9/26; M0/M1=56/1; OS/DSS events=26/17.",
    paste0("P01 positive marker coverage: ", score_median$coverage$positive_available[score_median$coverage$policy == "P01"], "/120; negative coverage: ", score_median$coverage$negative_available[score_median$coverage$policy == "P01"], "/120."),
    paste0("Spearman rho for P01 median probe aggregation versus highest-IQR probe: ", signif(iqr_rho, 4L), "."),
    "P00/P01/P11 were separately re-standardized as Z scores within MDA57; Z scores across all 116 samples were used only for expression platform auditing.",
    "",
    "2. Pathological associations",
    paste0("Three-level T-stage overall P=", if (nrow(p_t)) format_p(p_t$kruskal_p[1L]) else "NA", "; trend P=", if (nrow(p_t)) format_p(p_t$spearman_p[1L]) else "NA", "."),
    paste0("Three-level N-status overall P=", if (nrow(p_n)) format_p(p_n$kruskal_p[1L]) else "NA", "; trend P=", if (nrow(p_n)) format_p(p_n$spearman_p[1L]) else "NA", "."),
    paste0("TN burden overall P=", if (nrow(p_tn)) format_p(p_tn$kruskal_p[1L]) else "NA", "; trend P=", if (nrow(p_tn)) format_p(p_tn$spearman_p[1L]) else "NA", "."),
    "T stage is fixed as pT1/pT2, pT3 and pT4; N status is fixed as N0, N1 and N2; TN burden is not AJCC stage.",
    "",
    "3. Primary survival results",
    paste0("DSS M4 (P01 + age + sex): HR=", signif(dss_m4$HR[1L], 4L), ", 95%CI ", signif(dss_m4$CI_low[1L], 4L), "-", signif(dss_m4$CI_high[1L], 4L), ", P=", format_p(dss_m4$p_wald[1L]), "."),
    paste0("DSS M5 (add T/N): HR=", signif(dss_m5$HR[1L], 4L), ", 95%CI ", signif(dss_m5$CI_low[1L], 4L), "-", signif(dss_m5$CI_high[1L], 4L), ", P=", format_p(dss_m5$p_wald[1L]), "."),
    paste0("OS M4: HR=", signif(os_m4$HR[1L], 4L), ", 95%CI ", signif(os_m4$CI_low[1L], 4L), "-", signif(os_m4$CI_high[1L], 4L), ", P=", format_p(os_m4$p_wald[1L]), "."),
    paste0("OS M5: HR=", signif(os_m5$HR[1L], 4L), ", 95%CI ", signif(os_m5$CI_low[1L], 4L), "-", signif(os_m5$CI_high[1L], 4L), ", P=", format_p(os_m5$p_wald[1L]), "."),
    "DSS is the prioritized endpoint for this cohort; OS is also reported. M4 is the primary demographic model, and M5 is the pathology-extended sensitivity model. The primary model is not changed based on P values.",
    "",
    "4. Prespecified sensitivity analyses",
    "01_IQR_highest_probe: highest-IQR probe per gene across 116 samples;",
    "02_policy_P00 and 03_policy_P11: LAMC2/ITGB4 marker policy sensitivity;",
    "04_M0_only: exclude the sole M1 case, retaining 56 M0 cases;",
    "05_strict_pT2_pT4: exclude two pT1 cases, retaining 55 strict pathological MIBC cases;",
    "06_Firth_optional: run only if coxphf is installed in scrna_full; no automatic installation;",
    "07_nonPH_timevarying: triggered only if the M4 P01 or global PH test has P<0.05;",
    "08_quartile_descriptive: quartile KM analysis for descriptive purposes only, not primary inference.",
    "",
    "5. Methodological boundaries",
    "P01 is a fixed single-cell reference-guided within-sample rank projection, not a cell proportion.",
    "The primary analysis does not change markers, aggregation rules, cohorts or models based on pathological outcomes, survival outcomes, P values or probe results.",
    "Survival time is described only as overall survival or disease-specific survival; the public materials do not fully establish time zero, so it must not be called postoperative survival.",
    "",
    paste0("All final QC checks passed: ", all(qc$passed)),
    "Upload only the final ZIP generated in this run directory."
  )
  writeLines(report, REPORT_FILE, useBytes = TRUE)
  final_status <- if (all(qc$passed)) "SUCCESS_GSE48276_P01_SECOND_PART" else "COMPLETED_WITH_QC_WARNINGS_GSE48276_P01_SECOND_PART"
  writeLines(final_status, STATUS_FILE, useBytes = TRUE)
  log_msg(final_status)
  log_msg("End: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
}

tryCatch(
  main(),
  error = function(e) {
    error_text <- conditionMessage(e)
    writeLines(paste0("FAILED_GSE48276_P01_SECOND_PART: ", error_text), STATUS_FILE, useBytes = TRUE)
    calls <- capture.output(sys.calls())
    writeLines(c(
      "GSE48276 P01 part 2 run failed.", paste0("Error: ", error_text), "", "Call stack: ", calls, "",
      "No previous results were deleted or modified.", "No R packages were automatically installed or upgraded.", "Neither R nor the terminal was closed.", paste0("Log: ", LOG_FILE)
    ), REPORT_FILE, useBytes = TRUE)
    log_msg("[ERROR] ", error_text)
    stop(error_text, call. = FALSE)
  }
)

