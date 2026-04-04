#' Generate UpSet Plots per Sample/Domain/Level (Step 002)
#'
#' "karioCaS never are upset!"
#' Generates UpSet plots with detailed logging.
#'
#' @param project_dir Path to the project root.
#'
#' @return Generates PDF plots saved in the project directory and returns a data frame invisibly.
#' @export
#' @importFrom dplyr filter mutate select distinct case_when pull
#' @importFrom tidyr pivot_wider
#' @importFrom UpSetR upset
#' @importFrom grDevices pdf dev.off
#' @importFrom grid grid.text gpar
#' @examples
#' # Get the path to the included toy dataset
#' toy_project <- system.file("extdata", "your_project_name", package = "karioCaS")
#'
#' # upset_kariocas(project_dir = toy_project)

upset_kariocas <- function(project_dir) {

  # ==============================================================================
  # 1. SETUP & LOGGING
  # ==============================================================================
  output_dir  <- file.path(project_dir, "002_UpSetComparison_Plots")
  log_dir     <- file.path(project_dir, "logs")

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)

  log_file <- file.path(log_dir, "log_002_upsetplots.txt")

  log_msg <- function(...) {
    msg <- paste0(...)
    message(msg)
    cat(paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", msg, "\n"),
        file = log_file, append = TRUE)
  }

  cat("====================================================\n", file = log_file)
  cat("LOG: 002_UPSET_PLOTS (karioCaS never are upset)\n", file = log_file, append = TRUE)
  cat("PROJECT DIR: ", project_dir, "\n", file = log_file, append = TRUE)
  cat("====================================================\n", file = log_file, append = TRUE)

  if (!requireNamespace("UpSetR", quietly = TRUE)) {
    stop("Package 'UpSetR' is required.")
  }

  # ==============================================================================
  # 2. DATA LOADING
  # ==============================================================================
  log_msg(">>> Loading Data...")
  df_long <- .get_tidy_data(project_dir)

  SAMPLES <- unique(df_long$sample)
  DOMAINS <- names(get_kariocas_colors("domains"))
  LEVELS  <- c("Species", "Genus", "Family")

  log_msg(">>> Starting UpSet Analysis for ", length(SAMPLES), " samples.")

  # ==============================================================================
  # 3. ANALYSIS LOOP
  # ==============================================================================
  for (samp in SAMPLES) {
    log_msg("------------------------------------------------")
    log_msg("  Processing Sample: ", samp)

    current_out_dir <- file.path(output_dir, samp)
    if (!dir.exists(current_out_dir)) dir.create(current_out_dir)

    for (dom in DOMAINS) {
      for (lvl in LEVELS) {

        # 3.1 DATA FILTERING
        df_sub <- df_long %>%
          dplyr::filter(sample == samp, Domain == dom, Rank == lvl, Counts > 0) %>%
          dplyr::select(CS, Taxon_Name) %>%
          dplyr::distinct()

        if (nrow(df_sub) == 0) {
          log_msg("    Skipping ", dom, "-", lvl, ": No data found.")
          next
        }

        # 3.2 BINARY MATRIX
        binary_matrix <- df_sub %>%
          dplyr::mutate(
            Present = 1,
            CS_Label = sprintf("CS%02d", as.numeric(CS))
          ) %>%
          dplyr::select(-CS) %>%
          tidyr::pivot_wider(
            names_from = CS_Label,
            values_from = Present,
            values_fill = 0
          ) %>%
          as.data.frame()

        upset_cols <- colnames(binary_matrix)[colnames(binary_matrix) != "Taxon_Name"]

        if (length(upset_cols) < 2) {
          log_msg("    Skipping ", dom, "-", lvl, ": Not enough intersection levels (Found: ", length(upset_cols), ")")
          next
        }

        # 3.3 PLOTTING
        file_name <- paste0(samp, "_", dom, "_", lvl, "_UpSet.pdf")
        full_path <- file.path(current_out_dir, file_name)

        grDevices::pdf(file = full_path,
                       width = get_kariocas_dims()$width,
                       height = get_kariocas_dims()$height,
                       onefile = FALSE)

        tryCatch({
          print(
            UpSetR::upset(
              binary_matrix,
              sets = rev(upset_cols),
              keep.order = TRUE,
              order.by = "freq",
              empty.intersections = NULL,
              mainbar.y.label = paste(lvl, "Intersections"),
              sets.x.label = paste("Total", lvl, "per CS"),
              text.scale = c(1.5, 1.5, 1.2, 1.2, 1.5, 1.3),
              mb.ratio = c(0.6, 0.4),
              main.bar.color = get_kariocas_colors("upset")$main,
              sets.bar.color = get_kariocas_colors("upset")$sets,
              matrix.color = "grey20",
              shade.color = "grey90"
            )
          )

          grid::grid.text(
            label = paste(samp, "-", dom, "|", lvl, "Intersection Analysis"),
            x = 0.5, y = 0.98,
            gp = grid::gpar(fontsize = 16, fontface = "bold")
          )

          log_msg("    -> Generated: ", file_name)

        }, error = function(e) {
          log_msg("    ERROR plotting ", file_name, ": ", e$message)
        })

        grDevices::dev.off()
      }
    }
  }

  log_msg("SUCCESS: UpSet analysis completed.")
  return(invisible(NULL))
}
