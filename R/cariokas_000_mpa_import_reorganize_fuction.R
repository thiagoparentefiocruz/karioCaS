#' Import and Process MPA Data (Project Structure Based)
#'
#' This function enforces a specific project structure. It takes the project root directory,
#' looks specifically for a folder named "000_mpa_original", processes the .mpa files therein,
#' and saves logs to "log_archives".
#'
#' @param project_dir Path to the root of the project (e.g., "/Users/name/project_x").
#' @param log_filename Name of the log file (Default: "log_000_data_import.txt").
#'
#' @return A data.frame containing the processed count matrix.
#' @export
#' @importFrom readr read_tsv write_rds write_tsv
#' @importFrom dplyr mutate case_when select relocate %>%
#' @importFrom tidyr pivot_wider
#' @importFrom stringr str_remove str_extract str_sub
#' @importFrom purrr map_dfr set_names

import_mpa_data <- function(project_dir, log_filename = "log_000_data_import.txt") {

  # ==============================================================================
  # 1. DEFINE PATHS BASED ON CONVENTION
  # ==============================================================================

  # Enforce specific folder names
  input_dir <- file.path(project_dir, "000_mpa_original")
  log_dir   <- file.path(project_dir, "log_archives")
  path_log  <- file.path(log_dir, log_filename)

  # ==============================================================================
  # 2. VALIDATION (FAIL FAST)
  # ==============================================================================

  # Check if Project Dir exists
  if (!dir.exists(project_dir)) {
    stop("CRITICAL ERROR: The Project Directory does not exist:\n", project_dir)
  }

  # Check if 000_mpa_original exists (Enforcing Structure)
  if (!dir.exists(input_dir)) {
    stop("STRUCTURE ERROR: The folder '000_mpa_original' was not found inside the project directory.\n",
         "Expected location: ", input_dir, "\n",
         "Please create this folder and place your .mpa files inside it.")
  }

  # Create Log directory if it doesn't exist
  if (!dir.exists(log_dir)) {
    dir.create(log_dir, recursive = TRUE)
  }

  # ==============================================================================
  # 3. INITIALIZE LOGGING
  # ==============================================================================

  log_con <- file(path_log, open = "wt")
  sink(log_con, type = "output")
  sink(log_con, type = "message")

  # Safety mechanism to always close logs
  on.exit({
    sink(type = "output")
    sink(type = "message")
    close(log_con)
  }, add = TRUE)

  cat("====================================================\n")
  cat("LOG: MPA DATA IMPORT FUNCTION (STRUCTURED)\n")
  cat("DATE:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("PROJECT DIR:", project_dir, "\n")
  cat("INPUT DIR:", input_dir, "\n")
  cat("====================================================\n\n")

  tryCatch({

    # ==============================================================================
    # 4. LIST FILES (WITH FILTER)
    # ==============================================================================

    # FILTER: Only read files ending in .mpa to avoid reading logs/RDS/TSV
    mpa_files <- list.files(path = input_dir,
                            pattern = "\\.mpa$",
                            full.names = TRUE,
                            ignore.case = TRUE)

    cat("FOUND:", length(mpa_files), ".mpa files.\n")

    if (length(mpa_files) == 0) {
      stop("No .mpa files found inside '000_mpa_original'. Please check your folder.")
    }

    # Print file list for audit
    cat(basename(mpa_files))
    cat("\n")

    # ==============================================================================
    # 5. READ DATA
    # ==============================================================================
    cat("Reading files...\n")

    combined_mpas <- mpa_files %>%
      set_names(basename(.)) %>%
      purrr::map_dfr(~ readr::read_tsv(.,
                                       col_names = c("Taxonomy", "Counts"),
                                       col_types = "cd",
                                       progress = FALSE,
                                       show_col_types = FALSE),
                     .id = "Sample")

    # ==============================================================================
    # 6. TRANSFORM TO MATRIX
    # ==============================================================================
    cat("Creating Matrix and Formatting Taxonomy...\n")

    final_matrix <- combined_mpas %>%
      dplyr::mutate(Sample = stringr::str_remove(Sample, "\\.mpa$")) %>%
      tidyr::pivot_wider(names_from = Sample, values_from = Counts, values_fill = 0) %>%
      dplyr::mutate(
        ultimo_segmento = stringr::str_extract(Taxonomy, "[^|]+$"),
        prefixo_tax = stringr::str_sub(ultimo_segmento, 1, 1),
        Tax_Level = dplyr::case_when(
          prefixo_tax == "d" ~ "Domain",
          prefixo_tax == "k" ~ "Kingdom",
          prefixo_tax == "p" ~ "Phylum",
          prefixo_tax == "c" ~ "Class",
          prefixo_tax == "o" ~ "Order",
          prefixo_tax == "f" ~ "Family",
          prefixo_tax == "g" ~ "Genus",
          prefixo_tax == "s" ~ "Species",
          TRUE ~ "Other"
        )
      ) %>%
      dplyr::select(-ultimo_segmento, -prefixo_tax) %>%
      dplyr::relocate(Tax_Level, .after = Taxonomy)

    # ==============================================================================
    # 7. SAVE CHECKPOINTS (RDS/TSV)
    # ==============================================================================
    # Saving directly into 000_mpa_original (input_dir)
    path_rds <- file.path(input_dir, "combined_mpas_matrix_taxon_level.rds")
    path_tsv <- file.path(input_dir, "combined_mpas_matrix_taxon_level.tsv")

    cat("\nSAVING CHECKPOINTS:\n")
    cat("  -> RDS:", path_rds, "\n")
    cat("  -> TSV:", path_tsv, "\n")

    readr::write_rds(final_matrix, path_rds)
    readr::write_tsv(final_matrix, path_tsv)

    # ==============================================================================
    # 8. FINAL VALIDATION
    # ==============================================================================
    n_expected <- length(mpa_files)
    n_observed <- ncol(final_matrix)

    # Matrix columns = Taxonomy + Tax_Level + Samples (N) -> Total = N + 2
    if ((n_observed - 2) == n_expected) {
      cat("\nSUCCESS: Import completed successfully.\n")
      cat("Final Dimensions:", nrow(final_matrix), "x", ncol(final_matrix), "\n")
      return(final_matrix)
    } else {
      stop("Integrity Check Failed: Matrix columns do not match input file count.")
    }

  }, error = function(e) {
    cat("\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
    cat("CRITICAL ERROR DETECTED:\n")
    cat(e$message, "\n")
    cat("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
    stop(e$message)
  })
}
