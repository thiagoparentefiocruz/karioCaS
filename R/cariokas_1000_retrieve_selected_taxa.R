#' Retrieve Selected Taxa with Domain-Specific Thresholds (Final Mosaic Step)
#'
#' This function creates a "biological mosaic" for each sample. It allows the user to
#' define specific Confidence Scores (CS) and minimum read counts for each Domain individually,
#' ensuring the optimal parameters are used for Bacteria, Archaea, Eukaryota, and Viruses.
#' It also supports filtering by a specific Taxonomic Level.
#'
#' @param project_dir Path to the project root.
#' @param tax_level Specific taxonomic level to filter (e.g., "Species", "Genus"). If NULL, keeps all levels.
#' @param CS_A Character. CS suffix for Archaea (e.g., "06").
#' @param reads_min_A Integer. Min reads for Archaea.
#' @param CS_B Character. CS suffix for Bacteria (e.g., "08").
#' @param reads_min_B Integer. Min reads for Bacteria.
#' @param CS_E Character. CS suffix for Eukaryota (e.g., "06").
#' @param reads_min_E Integer. Min reads for Eukaryota.
#' @param CS_V Character. CS suffix for Viruses (e.g., "04").
#' @param reads_min_V Integer. Min reads for Viruses.
#'
#' @return Returns TRUE invisibly. Generates mosaic files in "1000_final_selection".
#' @export
#' @importFrom readr read_rds write_tsv write_delim
#' @importFrom dplyr select filter rename matches bind_rows mutate case_when
#' @importFrom stringr str_detect str_remove str_extract
#' @importFrom tidyselect all_of

retrieve_selected_taxa <- function(project_dir,
                                   tax_level = NULL,
                                   CS_A, reads_min_A,
                                   CS_B, reads_min_B,
                                   CS_E, reads_min_E,
                                   CS_V, reads_min_V) {

  # ==============================================================================
  # 1. SETUP & DIRECTORIES
  # ==============================================================================

  input_dir   <- file.path(project_dir, "000_mpa_original")
  output_dir  <- file.path(project_dir, "1000_final_selection")
  log_dir     <- file.path(project_dir, "log_archives")
  log_file    <- file.path(log_dir, "log_1000_mosaic_retrieval.txt")

  if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  if(!dir.exists(log_dir))    dir.create(log_dir, recursive = TRUE)

  # ==============================================================================
  # 2. LOGGING
  # ==============================================================================

  log_con <- file(log_file, open = "wt")
  sink(log_con, type = "output")
  sink(log_con, type = "message")

  on.exit({
    sink(type = "output")
    sink(type = "message")
    close(log_con)
  }, add = TRUE)

  cat("====================================================\n")
  cat("LOG: 1000_RETRIEVE_SELECTED_TAXA (MOSAIC MODE)\n")
  cat("DATE:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("TAX LEVEL:", ifelse(is.null(tax_level), "ALL", tax_level), "\n")
  cat("PARAMS:\n")
  cat("  Archaea   -> CS:", CS_A, "| Min:", reads_min_A, "\n")
  cat("  Bacteria  -> CS:", CS_B, "| Min:", reads_min_B, "\n")
  cat("  Eukaryota -> CS:", CS_E, "| Min:", reads_min_E, "\n")
  cat("  Viruses   -> CS:", CS_V, "| Min:", reads_min_V, "\n")
  cat("====================================================\n\n")

  tryCatch({

    # ==============================================================================
    # 3. INPUT DATA
    # ==============================================================================
    cat("STEP 1: Loading Matrix...\n")

    file_rds <- file.path(input_dir, "combined_mpas_matrix_taxon_level.rds")

    if (!file.exists(file_rds)) {
      cat("  -> RDS not found. Calling internal Import Function...\n")
      import_mpa_data(project_dir = project_dir)
    }

    INPUT_MATRIX <- readr::read_rds(file_rds)
    cat("  -> Matrix loaded. Dimensions:", dim(INPUT_MATRIX), "\n")

    # Identificar colunas de taxonomia
    tax_col_name <- if("Taxonomia" %in% names(INPUT_MATRIX)) "Taxonomia" else "Taxonomy"
    tax_lvl_col  <- if("Nivel_Final" %in% names(INPUT_MATRIX)) "Nivel_Final" else "Tax_Level"

    # 3.1. Filter by Tax_Level (if requested)
    if (!is.null(tax_level)) {
      cat("  -> Filtering for Tax_Level:", tax_level, "\n")
      INPUT_MATRIX <- INPUT_MATRIX %>% dplyr::filter(!!as.symbol(tax_lvl_col) == tax_level)
      cat("     Rows remaining:", nrow(INPUT_MATRIX), "\n")
    }

    # 3.2. Identify Domains in the data
    # Create a helper column for simplified Domain ID
    INPUT_MATRIX <- INPUT_MATRIX %>%
      dplyr::mutate(Domain_Code = dplyr::case_when(
        stringr::str_detect(!!as.symbol(tax_col_name), "d__Archaea")   ~ "A",
        stringr::str_detect(!!as.symbol(tax_col_name), "d__Bacteria")  ~ "B",
        stringr::str_detect(!!as.symbol(tax_col_name), "d__Eukaryota") ~ "E",
        stringr::str_detect(!!as.symbol(tax_col_name), "d__Viruses")   ~ "V",
        TRUE ~ "Other"
      ))

    # ==============================================================================
    # 4. IDENTIFY SAMPLES
    # ==============================================================================
    # Extract unique sample names (basenames) from columns
    # We look for columns ending in _CS... and remove the suffix
    all_cols <- colnames(INPUT_MATRIX)
    sample_cols_raw <- all_cols[stringr::str_detect(all_cols, "_CS[0-9.]+$")]
    sample_basenames <- unique(stringr::str_remove(sample_cols_raw, "_CS[0-9.]+$"))

    cat("  -> Unique Samples found:", length(sample_basenames), "\n\n")

    # ==============================================================================
    # 5. MOSAIC CONSTRUCTION LOOP
    # ==============================================================================
    cat("STEP 2: Constructing Mosaic Profiles...\n")

    # Config Map for iteration
    domain_config <- list(
      "A" = list(CS = CS_A, Min = reads_min_A, Name = "Archaea"),
      "B" = list(CS = CS_B, Min = reads_min_B, Name = "Bacteria"),
      "E" = list(CS = CS_E, Min = reads_min_E, Name = "Eukaryota"),
      "V" = list(CS = CS_V, Min = reads_min_V, Name = "Viruses")
    )

    domain_order <- c("A", "B", "E", "V") # Order requested: Archaea, Bacteria, Eukaryota, Viruses

    for (samp in sample_basenames) {
      cat("  Processing Sample:", samp, "...\n")

      mosaic_parts <- list()

      for (dom_code in domain_order) {
        cfg <- domain_config[[dom_code]]

        # 1. Define Target Column for this Domain/Sample
        target_col <- paste0(samp, "_CS", cfg$CS)

        if (!target_col %in% colnames(INPUT_MATRIX)) {
          cat("     [WARNING] Column", target_col, "not found for", cfg$Name, "- Skipping Domain.\n")
          next
        }

        # 2. Filter Rows for this Domain AND Thresholds
        part_df <- INPUT_MATRIX %>%
          dplyr::filter(Domain_Code == dom_code) %>%
          dplyr::select(dplyr::all_of(c(tax_col_name, target_col))) %>%
          dplyr::rename(Taxonomy = !!tax_col_name, Counts = !!target_col) %>%
          dplyr::filter(Counts >= cfg$Min)

        if (nrow(part_df) > 0) {
          # cat("     ->", cfg$Name, "(CS", cfg$CS, ", >", cfg$Min, "):", nrow(part_df), "taxa.\n")
          mosaic_parts[[dom_code]] <- part_df
        }
      }

      # Combine all parts
      if (length(mosaic_parts) > 0) {
        final_df <- dplyr::bind_rows(mosaic_parts)

        # Rename 'Counts' to the Sample Name (as requested)
        final_df <- final_df %>%
          dplyr::rename(!!samp := Counts)

        cat("     -> Final Mosaic Size:", nrow(final_df), "taxa.\n")

        # Define Filenames
        lvl_tag <- if(is.null(tax_level)) "AllLevels" else tax_level
        base_name <- paste0(samp, "_Mosaic_", lvl_tag)
        path_mpa  <- file.path(output_dir, paste0(base_name, ".mpa"))
        path_tsv  <- file.path(output_dir, paste0(base_name, ".tsv"))

        readr::write_delim(final_df, path_mpa, delim = "\t")
        readr::write_tsv(final_df, path_tsv)

      } else {
        cat("     -> [WARNING] No taxa retained for this sample across any domain.\n")
      }
    }

    cat("\n====================================================\n")
    cat("SUCCESS: Mosaic extraction completed.\n")
    cat("====================================================\n")

    return(invisible(TRUE))

  }, error = function(e) {
    cat("\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
    cat("CRITICAL ERROR:\n")
    cat(e$message, "\n")
    cat("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
    stop(e$message)
  })
}
