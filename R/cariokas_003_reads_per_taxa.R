#' Generate Cutoff Saturation Analysis (Step 003)
#'
#' Generates publication-quality saturation curves.
#'
#' @param project_dir Path to the project root.
#' @param analysis_level Taxonomic level to analyze. Default: "Species".
#' @param cutoff_mode Analysis mode: "Extinction", "Stabilization", or "Both".
#' @param import_script_path Deprecated.
#'
#' @export
#' @importFrom readr read_rds
#' @importFrom dplyr rename_with filter mutate case_when select bind_rows
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_line geom_point scale_color_manual scale_shape_manual scale_linetype_manual scale_y_continuous labs theme scale_x_log10 scale_x_continuous ggsave
#' @importFrom patchwork plot_layout plot_annotation
#' @importFrom scales percent_format label_number

reads_per_taxa <- function(project_dir,
                           analysis_level = "Species",
                           cutoff_mode = "Both",
                           import_script_path = NULL) {

  # ==============================================================================
  # 1. SETUP
  # ==============================================================================
  input_dir  <- file.path(project_dir, "000_mpa_original")
  output_dir <- file.path(project_dir, "003_cutoffs")

  if (!dir.exists(input_dir)) stop("Input directory not found.")
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  rds_file <- file.path(input_dir, "combined_mpas_matrix_taxon_level.rds")
  df_long <- load_and_process_mpa(rds_file)

  # Prepare Cutoffs
  cutoffs <- c(0, 1, 5, 10, 50, 100, 500, 1000, 5000, 10000)
  DOMAINS <- c("Bacteria", "Archaea", "Viruses", "Eukaryota")
  SAMPLES <- unique(df_long$Sample)
  CS_LIST <- unique(df_long$CS)

  modes_to_run <- if(cutoff_mode == "Both") c("Extinction", "Stabilization") else cutoff_mode

  cat("\n=== Starting Saturation Analysis (Step 003) ===\n")

  for (samp in SAMPLES) {
    cat("  Processing:", samp, "\n")
    samp_out_dir <- file.path(output_dir, samp)
    if (!dir.exists(samp_out_dir)) dir.create(samp_out_dir)

    df_samp <- df_long %>%
      dplyr::filter(Sample == samp, Level == analysis_level)

    for (cs in CS_LIST) {
      df_cs <- df_samp %>% dplyr::filter(CS == cs)

      # Calculation Logic (Kept Intact)
      results_list <- list()
      total_reads_initial <- sum(df_cs$Counts)
      total_taxa_initial  <- dplyr::n_distinct(df_cs$Taxon)

      for (co in cutoffs) {
        filtered <- df_cs %>% dplyr::filter(Counts > co)
        retained_reads <- sum(filtered$Counts)
        retained_taxa  <- dplyr::n_distinct(filtered$Taxon)
        results_list[[as.character(co)]] <- data.frame(
          Cutoff = co,
          Retained_Reads_Pct = (retained_reads / total_reads_initial),
          Retained_Taxa_Pct  = (retained_taxa / total_taxa_initial),
          Domain = unique(df_cs$Domain) # Simplify for plotting structure
        )
      }
      df_res <- dplyr::bind_rows(results_list)

      # --- PLOTTING (Refactored) ---
      for (current_mode in modes_to_run) {
        plot_list <- list()

        for (dom in DOMAINS) {
          # Mock data structure for plotting consistency across domains
          # In a real run, you'd filter df_res by Domain, but here assuming structure for brevity
          # Logic adapted to match 001/003 structure:
          dat_dom <- df_res[df_res$Domain == dom, ] # Simplified filter if data allows

          # Since the original script calculated per CS/Sample globally,
          # we assume df_res contains the domain info.
          # If original script logic was "Process whole CS then split by Domain", we follow:
          df_dom_calc <- df_cs %>% dplyr::filter(Domain == dom)

          # Recalculate specific for Domain to be precise
          res_dom_list <- list()
          t_reads <- sum(df_dom_calc$Counts)
          t_taxa  <- dplyr::n_distinct(df_dom_calc$Taxon)

          if(t_taxa == 0) {
            plot_list[[dom]] <- patchwork::wrap_plots(ggplot2::ggplot() + ggplot2::theme_void())
            next
          }

          for (co in cutoffs) {
            f <- df_dom_calc %>% dplyr::filter(Counts > co)
            res_dom_list[[as.character(co)]] <- data.frame(
              Cutoff = co,
              Reads = sum(f$Counts) / t_reads,
              Taxa  = dplyr::n_distinct(f$Taxon) / t_taxa,
              Domain = dom
            )
          }
          plot_data <- dplyr::bind_rows(res_dom_list) %>%
            tidyr::pivot_longer(cols = c("Reads", "Taxa"), names_to = "Metric", values_to = "Value")

          # Plot Construction
          p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Cutoff, y = Value, color = Domain, linetype = Metric)) +
            ggplot2::geom_line(linewidth = 1.2) +
            ggplot2::geom_point(size = 2) +
            # Use Centralized Styles
            ggplot2::scale_color_manual(values = karioCaS_cols$domain) +
            ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1.05)) +
            ggplot2::scale_x_continuous(trans = "log1p", breaks = cutoffs) +
            ggplot2::labs(title = dom, y = "Retention (%)", x = "Read Count Cutoff (>N)") +
            theme_karioCaS() + # Use Central Theme
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

          if(current_mode == "Stabilization") {
            p <- p + ggplot2::coord_cartesian(ylim = c(0.8, 1)) # Zoom in
          }

          plot_list[[dom]] <- p
        }

        final_layout <- (plot_list[["Bacteria"]] | plot_list[["Archaea"]]) /
          (plot_list[["Eukaryota"]] | plot_list[["Viruses"]]) +
          patchwork::plot_annotation(
            title = paste(samp, "-", cs, "| Saturation:", current_mode),
            subtitle = paste("Retention of Reads vs", analysis_level),
            theme = ggplot2::theme(plot.title = ggplot2::element_text(face="bold", size=14, hjust=0.5))
          )

        file_name <- paste0(samp, "_", cs, "_Cutoff_", analysis_level, "_", current_mode, ".pdf")
        out_path <- file.path(samp_out_dir, file_name)
        ggplot2::ggsave(out_path, final_layout, width = 11.69, height = 8.27, units = "in")
      }
    }
  }
  cat("\nSUCCESS: Saturation plots generated.\n")
  return(TRUE)
}
