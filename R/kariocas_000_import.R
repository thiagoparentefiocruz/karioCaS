#' Import Kraken MPA Reports to TreeSummarizedExperiment (Step 000)
#'
#' Reads Kraken2 MPA-style reports directly from the '000_mpa_original' folder.
#' Parses filenames like 'Sample_CSXX.mpa' (e.g., PILO_CS09.mpa) to extract
#' Sample Name and Confidence Score.
#' Parses Taxonomy into 8 ranks: Domain (d__), Kingdom (k__), Phylum (p__)... Species (s__).
#' Determines the specific 'Rank' of each row.
#'
#' @param project_dir Path to the project root. The script expects a '000_mpa_original' folder inside.
#'
#' @return A TreeSummarizedExperiment object (invisibly).
#' @export
#' @importFrom utils read.table
#' @importFrom dplyr bind_rows mutate select distinct left_join case_when
#' @importFrom tidyr pivot_wider
#' @importFrom stringr str_match str_extract str_remove str_detect
#' @importFrom TreeSummarizedExperiment TreeSummarizedExperiment
#' @importFrom SummarizedExperiment colData rowData
#' @importFrom S4Vectors DataFrame
#' @importFrom readr write_tsv

import_karioCaS <- function(project_dir) {

  # ==============================================================================
  # 1. DEFINE PATHS & VALIDATION
  # ==============================================================================
  if (!dir.exists(project_dir)) stop("Project directory not found: ", project_dir)

  input_dir  <- file.path(project_dir, "000_mpa_original")
  output_dir <- file.path(project_dir, "000_karioCaS_input_matrix")

  if (!dir.exists(input_dir)) stop("Input folder missing: ", input_dir)

  if (!dir.exists(output_dir)) {
    message(">>> Creating output directory: ", output_dir)
    dir.create(output_dir, recursive = TRUE)
  }

  message(">>> [Bioconductor Integration] Scanning .mpa files in: ", input_dir)

  files <- list.files(input_dir, full.names = TRUE, pattern = "\\.mpa$|\\.txt$|\\.report$")
  if (length(files) == 0) stop("No files found in ", input_dir)

  long_data_list <- list()
  fname_regex <- "^(.*)_CS([0-9]+)\\.(mpa|txt|report)$"

  # ==============================================================================
  # 2. READ LOOP
  # ==============================================================================
  for (f in files) {
    fname <- basename(f)
    matches <- stringr::str_match(fname, fname_regex)

    if (is.na(matches[1, 1])) {
      warning("Skipping file with unexpected name format: ", fname)
      next
    }

    samp_name <- matches[1, 2]
    cs_str    <- matches[1, 3]
    cs_raw    <- as.numeric(cs_str)

    # Normalize CS logic (09 -> 90, 00 -> 0)
    if (cs_raw < 10 && nchar(cs_str) == 2 && substr(cs_str, 1, 1) == "0") {
      cs_val_num <- cs_raw * 10
    } else {
      cs_val_num <- cs_raw
    }
    if (cs_raw == 0) cs_val_num <- 0

    raw <- tryCatch(
      utils::read.table(f, sep = "\t", stringsAsFactors = FALSE, quote = ""),
      error = function(e) NULL
    )

    if (is.null(raw) || nrow(raw) == 0) next

    colnames(raw) <- c("Taxonomy", "Counts")
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

  if (length(long_data_list) == 0) stop("No valid data could be read.")

  message(">>> Assembling Data Matrix...")
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
  # 4. PARSE TAXONOMY (FULL HIERARCHY)
  # ==============================================================================
  message(">>> Parsing Taxonomy (8 Levels)...")
  tax_strings <- rownames(count_mat)

  extract_rank <- function(x, pattern) {
    part <- stringr::str_extract(x, paste0(pattern, "__[^|]+"))
    stringr::str_remove(part, paste0(pattern, "__"))
  }

  # Determine the lowest Rank for each row
  # Logic: Check the last tag in the string
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
    Kingdom = extract_rank(tax_strings, "k"), # Kingdom is separate now!
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
  # 5. SAVE AUDIT TSV (COMPLETE STRUCTURE)
  # ==============================================================================
  # We want the audit file to look like the matrix: Taxonomy | Ranks... | Samples...

  # Combine row data with counts
  audit_df <- dplyr::bind_cols(row_data_df, as.data.frame(count_mat))

  audit_file <- file.path(output_dir, "karioCaS_matrix_audit.tsv")
  message(">>> Saving Human-Readable Audit File: ", audit_file)
  readr::write_tsv(audit_df, audit_file)

  # ==============================================================================
  # 6. CREATE & SAVE TSE
  # ==============================================================================
  # Metadata
  meta_map <- df_big %>%
    dplyr::select(Col_ID, Sample, CS) %>%
    dplyr::distinct()

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

  message("SUCCESS: TreeSummarizedExperiment saved at: ", save_path)

  return(invisible(tse))
}
