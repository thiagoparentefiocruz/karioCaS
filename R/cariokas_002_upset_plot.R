#' Generate UpSet Plots per Sample/Domain/Level (Step 002)
#'
#' "karioCaS never are upset!"
#'
#' @param project_dir Path to the project root.
#' @param import_script_path Deprecated.
#'
#' @export
#' @importFrom readr read_rds
#' @importFrom dplyr rename_with filter mutate case_when distinct select pull
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom UpSetR upset
#' @importFrom grid grid.text gpar
#' @importFrom grDevices pdf dev.off

upset_karioCaS_never_are <- function(project_dir, import_script_path = NULL) {

  # ==============================================================================
  # 1. SETUP
  # ==============================================================================
  input_dir   <- file.path(project_dir, "000_mpa_original")
  output_dir  <- file.path(project_dir, "002_UpSetComparison_Plots")

  if (!dir.exists(input_dir)) stop("Input directory not found.")
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  rds_file <- file.path(input_dir, "combined_mpas_matrix_taxon_level.rds")
  df_long <- load_and_process_mpa(rds_file)

  SAMPLES <- unique(df_long$Sample)
  DOMAINS <- c("Bacteria", "Archaea", "Viruses", "Eukaryota")
  LEVELS  <- c("Species", "Genus", "Family")

  cat("\n=== Starting UpSet Analysis (Step 002) ===\n")

  for (samp in SAMPLES) {
    cat("  Processing Sample:", samp, "\n")
    current_out_dir <- file.path(output_dir, samp)
    if (!dir.exists(current_out_dir)) dir.create(current_out_dir)

    df_samp <- df_long %>% dplyr::filter(Sample == samp)

    for (dom in DOMAINS) {
      for (lvl in LEVELS) {

        # Prepare Binary Matrix
        df_sub <- df_samp %>%
          dplyr::filter(Domain == dom, Level == lvl) %>%
          dplyr::mutate(Present = 1) %>%
          dplyr::select(Taxon, CS, Present) %>%
          dplyr::distinct() %>%
          tidyr::pivot_wider(names_from = CS, values_from = Present, values_fill = 0)

        if (nrow(df_sub) < 2) next

        binary_matrix <- as.data.frame(df_sub)
        upset_cols <- colnames(binary_matrix)[-1] # Exclude Taxon column

        file_name <- paste0(samp, "_", dom, "_", lvl, "_UpSet.pdf")
        full_path <- file.path(current_out_dir, file_name)

        # Plotting (Using Centralized Color)
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
                          mb.ratio = c(0.6, 0.4),
                          # Applying karioCaS visual identity:
                          main.bar.color = karioCaS_cols$main,
                          sets.bar.color = "grey30"
            )
          )
          grid::grid.text(paste(samp, "-", dom, lvl), x = 0.65, y = 0.95, gp = grid::gpar(fontsize = 15))

        }, error = function(e) {
          cat("      [ERROR PLOTTING]:", e$message, "\n")
        })

        grDevices::dev.off()
      }
    }
  }
  cat("\nSUCCESS: All UpSet plots generated.\n")
  return(invisible(TRUE))
}
