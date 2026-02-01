#' Generate UpSet Plots per Sample/Domain/Level (Step 002)
#'
#' "karioCaS never are upset!"
#'
#' This function generates UpSet plots to visualize the intersection of taxa
#' recovered across different Confidence Scores (CS). It helps identify which taxa
#' are robust (present in all CS) versus fragile (lost at higher stringency).
#'
#' @param project_dir Path to the project root.
#'
#' @export
#' @importFrom readr read_rds
#' @importFrom dplyr filter mutate select distinct case_when pull
#' @importFrom tidyr pivot_wider
#' @importFrom UpSetR upset
#' @importFrom grDevices pdf dev.off
#' @importFrom grid grid.text gpar

upset_kariocas <- function(project_dir) {

  # ==============================================================================
  # 1. SETUP & DIRECTORIES
  # ==============================================================================
  input_dir   <- file.path(project_dir, "000_karioCaS_input_matrix")
  output_dir  <- file.path(project_dir, "002_UpSetComparison_Plots")

  if (!dir.exists(input_dir)) {
    stop("CRITICAL ERROR: Input directory not found: ", input_dir,
         "\nPlease ensure step 000 (import) was run successfully.")
  }

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    message("Created output directory: ", output_dir)
  }

  rds_file <- file.path(input_dir, "karioCaS_input_matrix.rds")
  if (!file.exists(rds_file)) stop("Input RDS file not found: ", rds_file)

  # Check for UpSetR dependency explicitly
  if (!requireNamespace("UpSetR", quietly = TRUE)) {
    stop("Package 'UpSetR' is required for this script. Please install it using install.packages('UpSetR').")
  }

  # ==============================================================================
  # 2. DATA LOADING & PREP
  # ==============================================================================
  message(">>> Loading Tidy Data from: ", rds_file)
  df_long <- readr::read_rds(rds_file)

  # Create a clean Taxon_Name for unique identification
  df_proc <- df_long %>%
    dplyr::mutate(
      Taxon_Name = dplyr::case_when(
        Rank == "Species" ~ Species,
        Rank == "Genus"   ~ Genus,
        Rank == "Family"  ~ Family,
        TRUE ~ NA_character_
      )
    )

  SAMPLES <- unique(df_proc$sample)
  DOMAINS <- names(kariocas_colors$domains) # Use standardized domains
  LEVELS  <- c("Species", "Genus", "Family") # Standard levels for UpSet

  message(">>> Starting UpSet Analysis for ", length(SAMPLES), " samples.")

  # ==============================================================================
  # 3. ANALYSIS LOOP
  # ==============================================================================
  for (samp in SAMPLES) {
    message("  Processing Sample: ", samp)

    # Create sample subfolder
    current_out_dir <- file.path(output_dir, samp)
    if (!dir.exists(current_out_dir)) dir.create(current_out_dir)

    for (dom in DOMAINS) {
      for (lvl in LEVELS) {

        # 3.1 DATA FILTERING ---------------------------------------------------
        # Filter for current Sample, Domain, and Rank
        # We also filter out 0 counts because UpSet is about PRESENCE
        df_sub <- df_proc %>%
          dplyr::filter(sample == samp, Domain == dom, Rank == lvl, Counts > 0) %>%
          dplyr::select(CS, Taxon_Name) %>%
          dplyr::distinct() # Ensure unique Taxon per CS

        if (nrow(df_sub) == 0) {
          # Skip if no data
          next
        }

        # 3.2 BINARY MATRIX GENERATION -----------------------------------------
        # UpSetR requires a binary matrix (Rows=Taxa, Cols=Sets/CS)
        # We use a dummy "Present" column to pivot
        binary_matrix <- df_sub %>%
          dplyr::mutate(
            Present = 1,
            # Format CS as string labels for the Sets (e.g., "CS00", "CS90")
            # We assume CS is numeric (0, 10, ...). We format to 2 digits.
            CS_Label = sprintf("CS%02d", as.numeric(CS))
          ) %>%
          dplyr::select(-CS) %>%
          tidyr::pivot_wider(
            names_from = CS_Label,
            values_from = Present,
            values_fill = 0
          ) %>%
          as.data.frame() # UpSetR works best with base data.frames

        # The first column is usually Taxon_Name (from pivot), usually we need to remove it
        # or set it as row names, but UpSetR takes the whole DF and we specify columns.
        # Let's clean it up:

        # Identify the intersection columns (the sets)
        upset_cols <- colnames(binary_matrix)[colnames(binary_matrix) != "Taxon_Name"]

        # Validation: Need at least 2 sets (CS levels) to show an intersection plot
        if (length(upset_cols) < 2) {
          message("    Skipping ", dom, "-", lvl, ": Not enough CS levels for intersection.")
          next
        }

        # 3.3 PLOTTING ---------------------------------------------------------
        file_name <- paste0(samp, "_", dom, "_", lvl, "_UpSet.pdf")
        full_path <- file.path(current_out_dir, file_name)

        # Open PDF Device (Using karioCaS dimensions)
        grDevices::pdf(file = full_path,
                       width = kariocas_dims$width,
                       height = kariocas_dims$height,
                       onefile = FALSE)

        tryCatch({
          # UpSetR Plot
          print(
            UpSetR::upset(
              binary_matrix,
              sets = rev(upset_cols), # Reverse to show CS00 first (usually intuitive)
              keep.order = TRUE,      # Keep the order of sets we defined
              order.by = "freq",      # Order intersections by frequency
              empty.intersections = NULL,

              # Labels
              mainbar.y.label = paste(lvl, "Intersections"),
              sets.x.label = paste("Total", lvl, "per CS"),

              # Scales (Text sizes) - Tuned for A4
              text.scale = c(1.5, 1.5, 1.2, 1.2, 1.5, 1.3),
              mb.ratio = c(0.6, 0.4),

              # karioCaS Visual Identity (Colors from Style)
              main.bar.color = kariocas_colors$upset$main,
              sets.bar.color = kariocas_colors$upset$sets,

              # Grid & Background
              matrix.color = "grey20",
              shade.color = "grey90"
            )
          )

          # Add Title manually (UpSetR doesn't support main title natively well)
          grid::grid.text(
            label = paste(samp, "-", dom, "|", lvl, "Intersection Analysis"),
            x = 0.5, y = 0.98,
            gp = grid::gpar(fontsize = 16, fontface = "bold")
          )

        }, error = function(e) {
          message("    Error plotting ", file_name, ": ", e$message)
        })

        grDevices::dev.off() # Close device

      } # End Level
    } # End Domain
  } # End Sample

  message("\nSUCCESS: UpSet analysis completed. Files saved in: ", output_dir)
  return(invisible(NULL))
}
