#' Import Kraken MPA Reports to TreeSummarizedExperiment (Step 000)
#'
#' Reads Kraken2 MPA-style reports directly from the '000_mpa_original' folder.
#' Parses filenames like 'Sample_CSXX.mpa' (e.g., PILO_CS09.mpa).
#' Generates a detailed log file for traceability.
#'
#' @param project_dir Path to the project root. The script expects a '000_mpa_original' folder inside.
#'
#' @return A TreeSummarizedExperiment object (invisibly).
#' @export
#' @importFrom utils read.table
#' @importFrom dplyr bind_rows mutate select distinct left_join case_when %>%
#' @importFrom tidyr pivot_wider
#' @importFrom stringr str_match str_extract str_remove str_detect
#' @importFrom TreeSummarizedExperiment TreeSummarizedExperiment
#' @importFrom SummarizedExperiment colData rowData
#' @importFrom S4Vectors DataFrame
#' @importFrom readr write_tsv
#' @examples
#' # Pega o caminho interno do pacote onde está o toy dataset
#' toy_project <- system.file("extdata", "your_project_name", package = "karioCaS")
#' 
#' # Roda a função de verdade para o fiscal ver
#' import_karioCaS(project_dir = toy_project)

import_karioCaS <- function(project_dir) {

  # ==============================================================================
  # 1. SETUP & LOGGING INIT
  # ==============================================================================
  if (!dir.exists(project_dir)) stop("Project directory not found: ", project_dir)

  input_dir  <- file.path(project_dir, "000_mpa_original")
  output_dir <- file.path(project_dir, "000_karioCaS_input_matrix")
  log_dir    <- file.path(project_dir, "logs")

  if (!dir.exists(input_dir)) stop("Input folder missing: ", input_dir)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)

  # Initialize Log File
  log_file <- file.path(log_dir, "log_000_data_import.txt")

  # Helper function to write to both console and file
  log_msg <- function(...) {
    msg <- paste0(...)
    message(msg)
    cat(paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", msg, "\n"),
        file = log_file, append = TRUE)
  }

  # Header
  cat("====================================================\n", file = log_file)
  cat("LOG: 000_DATA_IMPORT (Structured)\n", file = log_file, append = TRUE)
  cat("PROJECT DIR: ", project_dir, "\n", file = log_file, append = TRUE)
  cat("====================================================\n", file = log_file, append = TRUE)

  log_msg(">>> [Bioconductor Integration] Scanning .mpa files in: ", input_dir)

  files <- list.files(input_dir, full.names = TRUE, pattern = "\\.mpa$|\\.txt$|\\.report$")
  if (length(files) == 0) {
    log_msg("CRITICAL ERROR: No files found in ", input_dir)
    stop("No files found.")
  }

  log_msg("FOUND: ", length(files), " potential files.")

  long_data_list <- list()
  fname_regex <- "^(.*)_CS([0-9]+)\\.(mpa|txt|report)$"

  # ==============================================================================
  # 2. READ LOOP
  # ==============================================================================
  for (f in files) {
    fname <- basename(f)
    matches <- stringr::str_match(fname, fname_regex)

    if (is.na(matches[1, 1])) {
      log_msg("WARNING: Skipping file with unexpected name format: ", fname)
      next
    }

    samp_name <- matches[1, 2]
    cs_str    <- matches[1, 3]
    cs_raw    <- as.numeric(cs_str)

    # Normalize CS logic
    if (cs_raw < 10 && nchar(cs_str) == 2 && substr(cs_str, 1, 1) == "0") {
      cs_val_num <- cs_raw * 10
    } else {
      cs_val_num <- cs_raw
    }
    if (cs_raw == 0) cs_val_num <- 0

    log_msg(" -> Reading: ", fname, " (interpreted as Sample: ", samp_name, " | CS: ", cs_val_num, ")")

    raw <- tryCatch(
      utils::read.table(f, sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "", fill = TRUE),
      error = function(e) {
        log_msg("ERROR: Failed to read file: ", fname, " - ", e$message)
        return(NULL)
      }
    )

    if (is.null(raw) || nrow(raw) == 0) {
      log_msg("WARNING: File is empty or skipped: ", fname)
      next
    }

    if (ncol(raw) < 2) {
      log_msg("WARNING: File has fewer than 2 columns, skipping: ", fname)
      next
    }

    raw <- raw[, 1:2]
    colnames(raw) <- c("Taxonomy", "Counts")
    raw$Counts <- as.numeric(raw$Counts)

    unique_col_id <- paste0(samp_name, "_CS", sprintf("%02d", cs_val_num))

    long_data_list[[unique_col_id]] <- data.frame(
      Taxonomy = raw$Taxonomy,
      Counts = raw$Counts,
      Sample = samp_name,
      CS = cs_val_num,
      Col_ID = unique_col_id,
      stringsAsFactors = FALSE
    )
  }

  if (length(long_data_list) == 0) {
    log_msg("CRITICAL: No valid data could be assembled.")
    stop("No valid data.")
  }

  log_msg(">>> Assembling Data Matrix...")
  df_big <- dplyr::bind_rows(long_data_list)

  # ==============================================================================
  # 3. CONSTRUCT ASSAY MATRIX
  # ==============================================================================
  mat_df <- df_big %>%
    dplyr::select(Taxonomy, Col_ID, Counts) %>%
    tidyr::pivot_wider(names_from = Col_ID, values_from = Counts, values_fill = 0)

  count_mat <- as.matrix(mat_df[, -1])
  rownames(count_mat) <- mat_df$Taxonomy

  # ==============================================================================
  # 4. PARSE TAXONOMY
  # ==============================================================================
  log_msg(">>> Parsing Taxonomy (8 Levels)...")
  tax_strings <- rownames(count_mat)

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
      TRUE ~ "Unknown"
    )
  }

  row_data_df <- data.frame(
    Taxonomy_Full = tax_strings,
    Domain  = extract_rank(tax_strings, "d"),
    Kingdom = extract_rank(tax_strings, "k"),
    Phylum  = extract_rank(tax_strings, "p"),
    Class   = extract_rank(tax_strings, "c"),
    Order   = extract_rank(tax_strings, "o"),
    Family  = extract_rank(tax_strings, "f"),
    Genus   = extract_rank(tax_strings, "g"),
    Species = extract_rank(tax_strings, "s"),
    Rank    = determine_rank(tax_strings),
    stringsAsFactors = FALSE
  )

  # ==============================================================================
  # 5. SAVE OUTPUTS & LOG
  # ==============================================================================
  audit_file <- file.path(output_dir, "karioCaS_matrix_audit.tsv")
  audit_df <- dplyr::bind_cols(row_data_df, as.data.frame(count_mat))
  readr::write_tsv(audit_df, audit_file)
  log_msg("SAVED TSV AUDIT: ", audit_file)

  # Create TSE
  meta_map <- df_big %>% dplyr::select(Col_ID, Sample, CS) %>% dplyr::distinct()
  meta_ordered <- meta_map[match(colnames(count_mat), meta_map$Col_ID), ]

  col_data <- S4Vectors::DataFrame(
    Sample_ID = meta_ordered$Sample,
    Confidence_Score = meta_ordered$CS,
    row.names = meta_ordered$Col_ID
  )

  tse <- TreeSummarizedExperiment::TreeSummarizedExperiment(
    assays = list(counts = count_mat),
    rowData = S4Vectors::DataFrame(row_data_df),
    colData = col_data
  )

  save_path <- file.path(output_dir, "karioCaS_TSE.rds")
  saveRDS(tse, save_path)

  log_msg("SAVED TSE OBJECT: ", save_path)
  log_msg("FINAL DIMENSIONS: ", nrow(tse), " Taxa x ", ncol(tse), " Observations")
  log_msg("SUCCESS: Import completed.")

  return(invisible(tse))
}
