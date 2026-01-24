#' Generate UpSet Plots per Sample/Domain/Level (Step 002)
#'
#' "karioCaS never are upset!"
#' This function processes the unified matrix and generates UpSet plots to visualize
#' intersections between Confidence Scores (CS).
#'
#' @param project_dir Path to the project root (e.g., "/Users/name/project_x").
#' @param import_script_path Deprecated. The import function is now internal to the package.
#'
#' @return Returns TRUE invisibly if successful. Generates PDFs in "002_UpSetComparison_Plots".
#' @export
#' @importFrom readr read_rds
#' @importFrom dplyr rename_with filter mutate case_when distinct select pull %>%
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom stringr str_detect str_remove str_extract
#' @importFrom tidyselect any_of
#' @importFrom UpSetR upset
#' @importFrom grid grid.text gpar
#' @importFrom grDevices pdf dev.off

upset_karioCaS_never_are <- function(project_dir,
                                     import_script_path = NULL) {

  # ==============================================================================
  # 1. SETUP & DIRECTORIES
  # ==============================================================================

  input_dir   <- file.path(project_dir, "000_mpa_original")
  output_dir  <- file.path(project_dir, "002_UpSetComparison_Plots")
  log_dir     <- file.path(project_dir, "log_archives")
  log_file    <- file.path(log_dir, "log_002_upsetplots.txt")

  # Create directories
  if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  if(!dir.exists(log_dir))    dir.create(log_dir, recursive = TRUE)

  # ==============================================================================
  # 2. START LOGGING
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
  cat("LOG: 002_UPSET_PLOTS_GENERATION (karioCaS never are upset)\n")
  cat("DATE/TIME:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("PROJECT DIR:", project_dir, "\n")
  cat("====================================================\n\n")

  tryCatch({

    # ==============================================================================
    # 3. INPUT LOGIC (SMART CHECK)
    # ==============================================================================
    cat("STEP 1: Checking Input Data...\n")

    file_rds <- file.path(input_dir, "combined_mpas_matrix_taxon_level.rds")
    INPUT_MATRIX <- NULL

    # 3.1. Check if RDS exists
    if (file.exists(file_rds)) {
      cat("  -> Input file found:", file_rds, "\n")
      cat("  -> Loading RDS object...\n")
      INPUT_MATRIX <- readr::read_rds(file_rds)

    } else {
      cat("  -> [WARNING] Input file NOT found at:", file_rds, "\n")
      cat("  -> Attempting to run internal Import Function (Step 000)...\n")

      # 3.2. Fallback to Internal Function
      # PACKAGE ADAPTATION: Calling the internal function directly
      tryCatch({
        INPUT_MATRIX <- import_mpa_data(project_dir = project_dir)
      }, error = function(e) {
        stop("CRITICAL: Could not run 'import_mpa_data'. Please ensure the previous step logic is correct.\nError: ", e$message)
      })
    }

    # 3.3. Validate Matrix
    cat("\n  -> MATRIX VALIDATION:\n")
    cat("     Observations (Taxa):", nrow(INPUT_MATRIX), "\n")
    cat("     Variáveis (Cols):", ncol(INPUT_MATRIX), "\n")

    # ==============================================================================
    # 4. DATA PREPARATION
    # ==============================================================================
    cat("STEP 2: Processing Data Structure...\n")

    # Standardize Headers
    df_std <- INPUT_MATRIX %>%
      dplyr::rename_with(~ "Taxonomy",  .cols = tidyselect::any_of(c("Taxonomia", "Taxonomy"))) %>%
      dplyr::rename_with(~ "Tax_Level", .cols = tidyselect::any_of(c("tax_level", "Tax_Level", "Nivel_Final")))

    # Pivot to Long
    cat("  -> Pivoting matrix to long format...\n")
    df_long <- df_std %>%
      tidyr::pivot_longer(
        cols = -c(Taxonomy, Tax_Level),
        names_to = "Raw_Sample_Col",
        values_to = "Counts"
      ) %>%
      dplyr::filter(Counts > 0)

    # Extract Metadata
    cat("  -> Extracting Metadata (Domains, Samples, CS)...\n")
    df_master <- df_long %>%
      dplyr::mutate(
        Domain = dplyr::case_when(
          stringr::str_detect(Taxonomy, "d__Viruses")   ~ "Viruses",
          stringr::str_detect(Taxonomy, "d__Bacteria")  ~ "Bacteria",
          stringr::str_detect(Taxonomy, "d__Archaea")   ~ "Archaea",
          stringr::str_detect(Taxonomy, "d__Eukaryota") ~ "Eukaryota",
          TRUE ~ "Unknown"
        ),
        Sample_Name = stringr::str_remove(Raw_Sample_Col, "_CS[0-9.]+$"),
        CS_Set_Label = stringr::str_extract(Raw_Sample_Col, "CS[0-9.]+$")
      ) %>%
      dplyr::mutate(
        Taxon_Clean_Name = stringr::str_extract(Taxonomy, "[^|]+$"),
        Taxon_Clean_Name = stringr::str_remove(Taxon_Clean_Name, "^[kpcofgs]__")
      )

    unique_samples <- unique(df_master$Sample_Name)
    unique_domains <- unique(df_master$Domain)
    target_levels  <- c("Species", "Genus")

    cat("  -> Unique Samples:", length(unique_samples), "\n")
    cat("  -> Unique Domains:", length(unique_domains), "\n\n")

    # ==============================================================================
    # 5. GENERATION LOOP (UPSET)
    # ==============================================================================
    cat("STEP 3: Generating UpSet Plots...\n")

    for (samp in unique_samples) {
      cat("  Processing Sample:", samp, "...\n")

      df_samp <- df_master %>% dplyr::filter(Sample_Name == samp)

      for (dom in unique_domains) {
        if(dom == "Unknown") next

        df_dom <- df_samp %>% dplyr::filter(Domain == dom)
        if(nrow(df_dom) == 0) next

        for (lvl in target_levels) {

          df_lvl <- df_dom %>% dplyr::filter(Tax_Level == lvl)

          if(nrow(df_lvl) == 0) next

          # --- BINARY MATRIX PREP ---
          binary_matrix <- df_lvl %>%
            dplyr::select(Taxon_Clean_Name, CS_Set_Label) %>%
            dplyr::distinct() %>%
            dplyr::mutate(Presence = 1L) %>%
            tidyr::pivot_wider(
              names_from = CS_Set_Label,
              values_from = Presence,
              values_fill = 0L
            ) %>%
            as.data.frame()

          upset_cols <- colnames(binary_matrix)[-1] # Remove Taxon column

          if(length(upset_cols) < 2) {
            cat("    [SKIP]", dom, lvl, "- Not enough CS sets (<2) for comparison.\n")
            next
          }

          # --- PLOT GENERATION ---
          cat("    -> Plotting:", dom, "-", lvl, "\n")

          current_out_dir <- file.path(output_dir, samp, dom)
          if(!dir.exists(current_out_dir)) dir.create(current_out_dir, recursive = TRUE)

          file_name <- paste0(samp, "_", dom, "_", lvl, "_UpSet.pdf")
          full_path <- file.path(current_out_dir, file_name)

          # A4 LANDSCAPE SETUP
          grDevices::pdf(file = full_path, width = 11.69, height = 8.27, onefile = FALSE)

          tryCatch({
            print(
              UpSetR::upset(binary_matrix,
                            sets = upset_cols,
                            order.by = "freq",
                            empty.intersections = NULL,
                            mainbar.y.label = paste(lvl, "Intersections"),
                            sets.x.label = "Total Taxa per CS",
                            text.scale = c(1.3, 1.3, 1, 1, 1.5, 1.2),
                            mb.ratio = c(0.6, 0.4)
              )
            )
            grid::grid.text(paste(samp, "-", dom, lvl), x = 0.65, y = 0.95, gp = grid::gpar(fontsize = 15))

          }, error = function(e) {
            cat("      [ERROR PLOTTING]:", e$message, "\n")
          })

          grDevices::dev.off() # Close PDF
        }
      }
    }

    cat("\n====================================================\n")
    cat("SUCCESS: All UpSet plots generated.\n")
    cat("====================================================\n")

    return(invisible(TRUE))

  }, error = function(e) {
    cat("\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
    cat("CRITICAL ERROR DETECTED:\n")
    cat(e$message, "\n")
    cat("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
    stop(e$message)
  })
}
