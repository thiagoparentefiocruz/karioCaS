#' Import and Process MPA Data (Hybrid Approach)
#'
#' Scans the '000_mpa_original' folder, parses filenames to extract Sample and CS,
#' reads the raw .mpa data (Taxonomy, Counts), parses the taxonomic string into
#' tidy columns (Domain...Species), infers the Rank, and saves two outputs:
#' 1. A Tidy Long RDS (for karioCaS internal use).
#' 2. A Wide TSV (for human inspection).
#'
#' @param project_dir Path to the root of the project.
#' @param log_filename Name of the log file (Default: "log_000_data_import.txt").
#'
#' @return A data.frame containing the processed tidy matrix (invisibly).
#' @export
#' @importFrom readr read_tsv write_rds write_tsv
#' @importFrom dplyr mutate select relocate filter bind_rows rename case_when group_by summarize ungroup
#' @importFrom tidyr pivot_wider
#' @importFrom stringr str_remove str_extract str_detect str_sub str_match
#' @importFrom purrr map_dfr
#' @importFrom fs dir_exists file_exists dir_create

import_mpa_data <- function(project_dir, log_filename = "log_000_data_import.txt") {

  # ==============================================================================
  # 1. DEFINE PATHS & SETUP
  # ==============================================================================
  input_dir  <- file.path(project_dir, "000_mpa_original")
  output_dir <- file.path(project_dir, "000_karioCaS_input_matrix")
  log_dir    <- file.path(project_dir, "log_archives")
  path_log   <- file.path(log_dir, log_filename)

  # Check Input Directory
  if (!dir.exists(input_dir)) {
    stop("CRITICAL ERROR: Input directory not found: ", input_dir,
         "\nPlease ensure the folder structure follows: YOUR_BASE_DIR/000_mpa_original/")
  }

  # Create Output and Log Directories
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)

  # Initialize Log
  log_con <- file(path_log, open = "wt")
  sink(log_con, type = "message")
  sink(log_con, type = "output")

  cat("==============================================================================\n")
  cat("LOG STARTED: Import MPA Data (Hybrid Tidy/Wide)\n")
  cat("Timestamp:", as.character(Sys.time()), "\n")
  cat("User:", Sys.info()[["user"]], "\n")
  cat("Project Directory:", project_dir, "\n")
  cat("Input Directory:", input_dir, "\n")
  cat("Output Directory:", output_dir, "\n")
  cat("==============================================================================\n\n")

  tryCatch({
    # ==============================================================================
    # 2. SCAN FILES
    # ==============================================================================
    # Expects pattern: SAMPLE_CSXX.mpa
    mpa_files <- list.files(input_dir, pattern = "\\.mpa$", full.names = TRUE, recursive = TRUE)

    if (length(mpa_files) == 0) {
      stop("No .mpa files found in ", input_dir)
    }

    cat("Found", length(mpa_files), "files to process.\n\n")

    # ==============================================================================
    # 3. READ, PARSE & STACK (LONG FORMAT CREATION)
    # ==============================================================================
    cat(">>> Processing files (Reading & Parsing metadata)...\n")

    # Iterate and bind rows immediately (More efficient than full joins)
    df_long_raw <- purrr::map_dfr(mpa_files, function(f_path) {

      # A. Extract Metadata from Filename
      f_name <- basename(f_path)

      # Regex Breakdown:
      # ^(.*)       -> Captures Sample Name (Group 1) - anything before _CS
      # _CS         -> Literal separator
      # ([0-9]+)    -> Captures the CS code XX (Group 2)
      # \.mpa$      -> File extension
      match_meta <- stringr::str_match(f_name, "^(.*)_CS([0-9]+)\\.mpa$")

      if (any(is.na(match_meta))) {
        cat("  [WARNING] Skipping file (invalid naming format):", f_name, "\n")
        return(NULL)
      }

      sample_id <- match_meta[1, 2]
      cs_code   <- match_meta[1, 3] # This is "00", "05", "10"

      # Logic: XX represents a decimal without the dot.
      # 00 = 0.0 -> 0%
      # 05 = 0.5 -> 50%
      # 10 = 1.0 -> 100%
      # Formula: Numeric Value * 10
      cs_numeric <- as.numeric(cs_code) * 10

      # B. Read Data (Assuming 2 cols: Taxonomy, Counts - No Header or Specific Header)
      # We force col_names to ensure consistency even if headers vary
      raw_data <- readr::read_tsv(f_path, col_names = c("Taxonomy", "Counts"),
                                  col_types = "ci", show_col_types = FALSE)

      # Remove header row if it exists (checks if 'Counts' column is numeric, if parsed as char header)
      # Since we forced 'ci' types, readr usually handles this, but robust cleaning:
      if(nrow(raw_data) > 0 && raw_data$Taxonomy[1] == "Taxonomy") {
        raw_data <- raw_data[-1, ]
      }

      # C. Add Metadata Columns
      raw_data %>%
        dplyr::mutate(
          sample = sample_id,
          CS = cs_numeric,
          Original_CS_Code = cs_code # Keep for reconstructing wide header later
        )
    })

    cat(">>> Raw data stacked. Total rows:", nrow(df_long_raw), "\n")

    # ==============================================================================
    # 4. TAXONOMY PARSING & ENRICHMENT (UTILS LOGIC)
    # ==============================================================================
    cat(">>> Parsing taxonomic strings and inferring Ranks...\n")

    df_tidy_long <- df_long_raw %>%
      dplyr::mutate(
        # --- TAXONOMY EXTRACTION ---
        Domain  = stringr::str_extract(Taxonomy, "d__[^|]+") %>% stringr::str_remove("d__"),
        Kingdom = stringr::str_extract(Taxonomy, "k__[^|]+") %>% stringr::str_remove("k__"),
        Phylum  = stringr::str_extract(Taxonomy, "p__[^|]+") %>% stringr::str_remove("p__"),
        Class   = stringr::str_extract(Taxonomy, "c__[^|]+") %>% stringr::str_remove("c__"),
        Order   = stringr::str_extract(Taxonomy, "o__[^|]+") %>% stringr::str_remove("o__"),
        Family  = stringr::str_extract(Taxonomy, "f__[^|]+") %>% stringr::str_remove("f__"),
        Genus   = stringr::str_extract(Taxonomy, "g__[^|]+") %>% stringr::str_remove("g__"),
        Species = stringr::str_extract(Taxonomy, "s__[^|]+") %>% stringr::str_remove("s__"),

        # --- RANK INFERENCE ---
        # Logic: Determine the deepest non-NA level found in the string
        Rank = dplyr::case_when(
          !is.na(Species) ~ "Species",
          !is.na(Genus)   ~ "Genus",
          !is.na(Family)  ~ "Family",
          !is.na(Order)   ~ "Order",
          !is.na(Class)   ~ "Class",
          !is.na(Phylum)  ~ "Phylum",
          !is.na(Kingdom) ~ "Kingdom",
          !is.na(Domain)  ~ "Domain",
          TRUE            ~ "Unknown"
        )
      ) %>%
      # Reorder Columns for final Long Output
      dplyr::select(
        Taxonomy, Domain, Kingdom, Phylum, Class, Order, Family, Genus, Species,
        Rank, sample, CS, Counts, Original_CS_Code
      )

    # ==============================================================================
    # 5. GENERATE OUTPUTS
    # ==============================================================================

    # OUTPUT 1: RDS (Tidy Long Format) - For Machine
    # Remove the helper column 'Original_CS_Code' for the RDS
    df_rds_output <- df_tidy_long %>% dplyr::select(-Original_CS_Code)

    path_rds <- file.path(output_dir, "karioCaS_input_matrix.rds")
    cat(">>> Saving RDS (Tidy Long):", path_rds, "\n")
    readr::write_rds(df_rds_output, path_rds)

    # OUTPUT 2: TSV (Wide Format) - For Human
    cat(">>> Generating Wide Format for TSV...\n")

    # Create the column header expected by humans: SAMPLE_CSXX
    df_wide <- df_tidy_long %>%
      dplyr::mutate(
        Wide_Header = paste0(sample, "_CS", Original_CS_Code)
      ) %>%
      dplyr::select(-sample, -CS, -Original_CS_Code) %>%
      # Pivot to Wide
      tidyr::pivot_wider(
        names_from = Wide_Header,
        values_from = Counts,
        values_fill = 0
      ) %>%
      # Arrange rows for readability (e.g., by Domain then Rank)
      dplyr::arrange(Domain, Rank)

    path_tsv <- file.path(output_dir, "karioCaS_input_matrix.tsv")
    cat(">>> Saving TSV (Wide):", path_tsv, "\n")
    readr::write_tsv(df_wide, path_tsv)

    # ==============================================================================
    # 6. FINAL STATS & CLOSE
    # ==============================================================================
    cat("\nSUCCESS: Import and Processing completed.\n")
    cat("Total Samples Processed:", length(unique(df_tidy_long$sample)), "\n")
    cat("Unique Taxa Found:", length(unique(df_tidy_long$Taxonomy)), "\n")
    cat("RDS Dimensions:", nrow(df_rds_output), "x", ncol(df_rds_output), "\n")
    cat("TSV Dimensions:", nrow(df_wide), "x", ncol(df_wide), "\n")

    return(invisible(df_rds_output))

  }, error = function(e) {
    cat("\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
    cat("CRITICAL ERROR DETECTED:\n")
    cat(conditionMessage(e), "\n")
    cat("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
    stop(e)
  }, finally = {
    sink(type = "message")
    sink(type = "output")
    close(log_con)
    message("Log saved to: ", path_log)
  })
}
