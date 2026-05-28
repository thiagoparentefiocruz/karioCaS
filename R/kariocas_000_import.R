# ==============================================================================
# PRIVATE HELPERS — import_karioCaS()
# ==============================================================================

#' @noRd
.imp_setup <- function(project_dir) {
  if (!dir.exists(project_dir)) stop("Project directory not found: ", project_dir)
  input_dir  <- file.path(project_dir, "000_mpa_original")
  output_dir <- file.path(project_dir, "000_karioCaS_input_matrix")
  log_dir    <- file.path(project_dir, "logs")
  if (!dir.exists(input_dir))  stop("Input folder missing: ", input_dir)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  if (!dir.exists(log_dir))    dir.create(log_dir,    recursive = TRUE)
  log_file <- file.path(log_dir, "log_000_data_import.txt")
  writeLines(c(
    "====================================================",
    "LOG: 000_DATA_IMPORT (Structured)",
    paste0("PROJECT DIR: ", project_dir),
    "===================================================="
  ), con = log_file)
  log_msg <- function(...) {
    msg <- paste0(...)
    message(msg)
    write(
      paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", msg),
      file = log_file, append = TRUE
    )
  }
  list(
    input_dir  = input_dir,
    output_dir = output_dir,
    log_msg    = log_msg
  )
}

#' @noRd
.imp_normalize_cs <- function(cs_raw, cs_str) {
  if (cs_raw == 0) return(0)
  if (cs_raw < 10 && nchar(cs_str) == 2 && substr(cs_str, 1, 1) == "0") {
    return(cs_raw * 10)
  }
  cs_raw
}

#' @noRd
.imp_read_one_file <- function(f, fname_regex, log_msg) {
  fname   <- basename(f)
  matches <- stringr::str_match(fname, fname_regex)
  if (is.na(matches[1, 1])) {
    log_msg("WARNING: Skipping file with unexpected name format: ", fname)
    return(NULL)
  }
  samp_name  <- matches[1, 2]
  cs_str     <- matches[1, 3]
  cs_val_num <- .imp_normalize_cs(as.numeric(cs_str), cs_str)
  log_msg(" -> Reading: ", fname,
          " (interpreted as Sample: ", samp_name, " | CS: ", cs_val_num, ")")
  raw <- tryCatch(
    utils::read.table(f, sep = "\t", stringsAsFactors = FALSE,
                      quote = "", comment.char = "", fill = TRUE),
    error = function(e) {
      log_msg("ERROR: Failed to read file: ", fname, " - ", e$message)
      NULL
    }
  )
  if (is.null(raw) || nrow(raw) == 0 || ncol(raw) < 2) {
    log_msg("WARNING: File is empty or invalid, skipping: ", fname)
    return(NULL)
  }
  raw    <- raw[, seq_len(2)]
  colnames(raw) <- c("Taxonomy", "Counts")
  raw$Counts <- as.numeric(raw$Counts)
  col_id <- paste0(samp_name, "_CS", sprintf("%02d", cs_val_num))
  data.frame(
    Taxonomy = raw$Taxonomy,
    Counts   = raw$Counts,
    Sample   = samp_name,
    CS       = cs_val_num,
    Col_ID   = col_id,
    stringsAsFactors = FALSE
  )
}

#' @noRd
.imp_read_files <- function(input_dir, log_msg) {
  log_msg(">>> [Bioconductor Integration] Scanning .mpa files in: ", input_dir)
  files <- list.files(
    input_dir, full.names = TRUE,
    pattern = "\\.mpa$|\\.txt$|\\.report$"
  )
  if (length(files) == 0) {
    log_msg("CRITICAL ERROR: No files found in ", input_dir)
    stop("No files found.")
  }
  log_msg("FOUND: ", length(files), " potential files.")
  fname_regex   <- "^(.*)_CS([0-9]+)\\.(mpa|txt|report)$"
  long_data_list <- Filter(Negate(is.null), lapply(files, function(f) {
    .imp_read_one_file(f, fname_regex, log_msg)
  }))
  if (length(long_data_list) == 0) {
    log_msg("CRITICAL: No valid data could be assembled.")
    stop("No valid data.")
  }
  log_msg(">>> Assembling Data Matrix...")
  dplyr::bind_rows(long_data_list)
}

#' @noRd
.imp_parse_taxonomy <- function(tax_strings, log_msg) {
  log_msg(">>> Parsing Taxonomy (8 Levels)...")
  extract_rank <- function(x, pattern) {
    part <- stringr::str_extract(x, paste0(pattern, "__[^|]+"))
    stringr::str_remove(part, paste0(pattern, "__"))
  }
  determine_rank <- function(x) {
    dplyr::case_when(
      stringr::str_detect(x, "\\|s__") ~ "Species",
      stringr::str_detect(x, "\\|g__") ~ "Genus",
      stringr::str_detect(x, "\\|f__") ~ "Family",
      stringr::str_detect(x, "\\|o__") ~ "Order",
      stringr::str_detect(x, "\\|c__") ~ "Class",
      stringr::str_detect(x, "\\|p__") ~ "Phylum",
      stringr::str_detect(x, "\\|k__") ~ "Kingdom",
      stringr::str_detect(x, "d__")    ~ "Domain",
      TRUE                              ~ "Unknown"
    )
  }
  data.frame(
    Taxonomy_Full = tax_strings,
    Domain        = extract_rank(tax_strings, "d"),
    Kingdom       = extract_rank(tax_strings, "k"),
    Phylum        = extract_rank(tax_strings, "p"),
    Class         = extract_rank(tax_strings, "c"),
    Order         = extract_rank(tax_strings, "o"),
    Family        = extract_rank(tax_strings, "f"),
    Genus         = extract_rank(tax_strings, "g"),
    Species       = extract_rank(tax_strings, "s"),
    Rank          = determine_rank(tax_strings),
    stringsAsFactors = FALSE
  )
}

#' @noRd
.imp_build_tse <- function(df_big, row_data_df, log_msg) {
  mat_df <- df_big |>
    dplyr::select("Taxonomy", "Col_ID", "Counts") |>
    tidyr::pivot_wider(
      names_from  = "Col_ID",
      values_from = "Counts",
      values_fill = 0
    )
  count_mat          <- as.matrix(mat_df[, -1])
  rownames(count_mat) <- mat_df$Taxonomy
  meta_map     <- dplyr::distinct(dplyr::select(df_big, "Col_ID", "Sample", "CS"))
  meta_ordered <- meta_map[match(colnames(count_mat), meta_map$Col_ID), ]
  col_data     <- S4Vectors::DataFrame(
    Sample_ID        = meta_ordered$Sample,
    Confidence_Score = meta_ordered$CS,
    row.names        = meta_ordered$Col_ID
  )
  TreeSummarizedExperiment::TreeSummarizedExperiment(
    assays  = list(counts = count_mat),
    rowData = S4Vectors::DataFrame(row_data_df),
    colData = col_data
  )
}

#' @noRd
.imp_save_outputs <- function(tse, row_data_df, count_mat,
                              output_dir, log_msg) {
  audit_file <- file.path(output_dir, "karioCaS_matrix_audit.tsv")
  readr::write_tsv(
    dplyr::bind_cols(row_data_df, as.data.frame(count_mat)),
    audit_file
  )
  log_msg("SAVED TSV AUDIT: ", audit_file)
  save_path <- file.path(output_dir, "karioCaS_TSE.rds")
  saveRDS(tse, save_path)
  log_msg("SAVED TSE OBJECT: ", save_path)
  log_msg("FINAL DIMENSIONS: ", nrow(tse), " Taxa x ", ncol(tse), " Observations")
}

# ==============================================================================
# EXPORTED FUNCTION
# ==============================================================================

#' Import Kraken MPA Reports to TreeSummarizedExperiment (Step 000)
#'
#' Reads Kraken2 MPA-style reports from the \code{000_mpa_original} folder.
#' Parses filenames like \code{SAMPLE_CSXX.mpa} (e.g., \code{PILO_CS09.mpa}).
#' Generates a detailed log file for traceability and saves a
#' \code{TreeSummarizedExperiment} object for downstream analysis.
#'
#' @param project_dir Path to the project root. Must contain a
#'   \code{000_mpa_original} subfolder with \code{.mpa} files.
#'
#' @return Invisibly returns a \code{TreeSummarizedExperiment} object.
#' @export
#' @importFrom utils read.table
#' @importFrom dplyr bind_rows mutate select distinct left_join case_when
#' @importFrom tidyr pivot_wider
#' @importFrom stringr str_match str_extract str_remove str_detect
#' @importFrom TreeSummarizedExperiment TreeSummarizedExperiment
#' @importFrom SummarizedExperiment colData rowData
#' @importFrom S4Vectors DataFrame
#' @importFrom readr write_tsv
#' @examples
#' toy_project <- system.file("extdata", "your_project_name", package = "karioCaS")
#' import_karioCaS(project_dir = toy_project)
import_karioCaS <- function(project_dir) {
  setup       <- .imp_setup(project_dir)
  df_big      <- .imp_read_files(setup$input_dir, setup$log_msg)
  mat_df      <- df_big |>
    dplyr::select("Taxonomy", "Col_ID", "Counts") |>
    tidyr::pivot_wider(
      names_from  = "Col_ID",
      values_from = "Counts",
      values_fill = 0
    )
  count_mat           <- as.matrix(mat_df[, -1])
  rownames(count_mat) <- mat_df$Taxonomy
  row_data_df <- .imp_parse_taxonomy(rownames(count_mat), setup$log_msg)
  tse         <- .imp_build_tse(df_big, row_data_df, setup$log_msg)
  .imp_save_outputs(tse, row_data_df, count_mat, setup$output_dir, setup$log_msg)
  setup$log_msg("SUCCESS: Import completed.")
  invisible(tse)
}