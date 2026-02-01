#' Run Confidence Score Retention Analysis (Step 001)
#'
#' Executes read retention analysis based on Confidence Score (Kraken/Bracken).
#' Calculates retention relative to the baseline (CS00) for both Reads and Distinct Taxa count.
#'
#' @param project_dir Root path of the project.
#'
#' @export
#' @importFrom dplyr filter select group_by summarise mutate left_join full_join n_distinct arrange pull case_when
#' @importFrom tidyr replace_na pivot_longer
#' @importFrom ggplot2 ggplot aes geom_line geom_point scale_color_manual scale_shape_manual scale_linetype_manual labs ggsave scale_x_continuous coord_cartesian guides guide_legend
#' @importFrom patchwork plot_layout plot_annotation
#' @importFrom scales percent
#' @importFrom ggtext element_markdown

taxa_retention <- function(project_dir) {

  # ==============================================================================
  # 1. SETUP & DIRECTORIES
  # ==============================================================================
  # Input logic is handled by .get_tidy_data helper
  output_dir <- file.path(project_dir, "001_taxa_retention")

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # ==============================================================================
  # 2. DATA LOADING
  # ==============================================================================
  message(">>> Loading Data (Auto-detected format)...")
  df_long <- .get_tidy_data(project_dir)

  # Define Ranks including Kingdom as requested
  target_ranks <- c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

  # Clean and Prepare Data
  # We simply filter for the valid ranks and ensure Rank is a Factor in correct order.
  # Note: Taxon_Name already exists in df_long coming from .get_tidy_data!
  df_proc <- df_long %>%
    dplyr::filter(Rank %in% target_ranks) %>%
    dplyr::mutate(Rank = factor(Rank, levels = target_ranks))

  # Validation
  if (nrow(df_proc) == 0) stop("No data found after filtering Ranks.")

  SAMPLES <- unique(df_proc$sample)
  DOMAINS <- names(kariocas_colors$domains) # Uses package colors

  message(">>> Starting Retention Analysis...")

  # ==============================================================================
  # 3. ANALYSIS LOOP
  # ==============================================================================
  for (samp in SAMPLES) {
    message("  Processing Sample: ", samp)

    plot_list <- list()

    # 3.1 Calculate Statistics per Domain -> Rank -> CS
    # We calculate Total Reads and Distinct Taxa Count
    df_samp <- df_proc %>%
      dplyr::filter(sample == samp)

    # Pre-calculate totals
    df_stats <- df_samp %>%
      dplyr::group_by(Domain, Rank, CS) %>%
      dplyr::summarise(
        total_reads = sum(Counts),
        n_taxa = dplyr::n_distinct(Taxon_Name),
        .groups = "drop"
      )

    # 3.2 Calculate Retention % (Baseline = CS 00)
    # Identify Baseline values (CS == 0)
    df_baseline <- df_stats %>%
      dplyr::filter(CS == 0) %>%
      dplyr::rename(base_reads = total_reads, base_taxa = n_taxa) %>%
      dplyr::select(Domain, Rank, base_reads, base_taxa)

    # Join and Calculate %
    df_plot_data <- df_stats %>%
      dplyr::left_join(df_baseline, by = c("Domain", "Rank")) %>%
      dplyr::mutate(
        pct_reads = (total_reads / base_reads) * 100,
        pct_taxa  = (n_taxa / base_taxa) * 100
      ) %>%
      # Remove NaN/Inf if baseline was 0
      dplyr::mutate(
        pct_reads = ifelse(base_reads == 0, 0, pct_reads),
        pct_taxa  = ifelse(base_taxa == 0, 0, pct_taxa)
      )

    # 3.3 PLOTTING
    # We generate plots for each Domain
    for (dom in DOMAINS) {

      df_dom <- df_plot_data %>% dplyr::filter(Domain == dom)

      # Handle Empty Data
      if (nrow(df_dom) == 0) {
        # Create Empty Plot placeholder
        p <- plot_kariocas_empty(
          title_text = dom,
          subtitle_text = "No data found",
          x_label = kariocas_labels$y_confidence,
          y_label = "Retention (%)"
        )
        plot_list[[dom]] <- p
        next
      }

      # Prepare for plotting: Melt metrics (Reads vs Taxa) to use as Facets or Colors
      # In the original script logic, we plot distinct lines for Ranks
      # Let's pivot to have "Metric" (Reads vs Taxa)

      # We create two plots per domain? Or combined?
      # The user script had `mode_tag`. Let's stick to the visual style requested.
      # The script seems to generate TWO files: one for Reads, one for Taxa.
      # Let's Loop over Modes
    } # End Domain Pre-Check

    # Actually, the logic needs to generate the plot FOR ALL DOMAINS inside a loop of MODES
    # to save separate PDFs like "CS_Retention_Reads.pdf" and "CS_Retention_Taxa.pdf"

    modes <- c("reads", "taxa")

    for (mode in modes) {
      plot_list <- list()

      for (dom in DOMAINS) {
        df_dom <- df_plot_data %>% dplyr::filter(Domain == dom)

        if (nrow(df_dom) == 0) {
          plot_list[[dom]] <- plot_kariocas_empty(dom, "No Data")
          next
        }

        # Select Y variable based on mode
        if (mode == "reads") {
          df_dom$Y_Val <- df_dom$pct_reads
          y_lab <- "Reads Retention (%)"
          subtitle_stats <- "Sum of reads normalized by CS00"
        } else {
          df_dom$Y_Val <- df_dom$pct_taxa
          y_lab <- "Taxa Retention (%)"
          subtitle_stats <- "Count of distinct taxa normalized by CS00"
        }

        # Plot Construction
        # X axis = CS, Y axis = Retention %, Color = Rank
        p <- ggplot2::ggplot(df_dom, ggplot2::aes(x = CS, y = Y_Val, color = Rank, group = Rank)) +
          ggplot2::geom_line(linewidth = 1) +
          ggplot2::geom_point(size = 2) +

          # Colors
          ggplot2::scale_color_manual(values = kariocas_colors$ranks) +

          # Labels
          ggplot2::labs(
            title = dom,
            subtitle = subtitle_stats,
            x = "Kraken Confidence Score (%)",
            y = y_lab,
            color = "Rank"
          ) +

          # Scale Y (Log10 hybrid or continuous?)
          # Retention starts at 100%. Let's use standard log scale if needed,
          # but usually retention plots are linear 0-100.
          # The snippet used `scale_y_kariocas_log10`. Let's respect it but check ranges.
          # Since it's percentage 0-100, standard linear is often better,
          # but if it drops to 0.001%, log is better.
          scale_y_kariocas_log10(limits = c(0.01, 105)) +

          ggplot2::scale_x_continuous(breaks = seq(0, 100, 20), limits = c(0, 100)) +

          theme_kariocas() +
          ggplot2::theme(
            legend.position = "bottom",
            plot.subtitle = ggplot2::element_text(size = 9, color = "grey40")
          )

        plot_list[[dom]] <- p
      }

      # Assemble 4-panel
      final_layout <- (plot_list[["Bacteria"]] | plot_list[["Archaea"]]) /
        (plot_list[["Eukaryota"]] | plot_list[["Viruses"]]) +
        patchwork::plot_annotation(
          title = paste(samp, "- Retention Analysis"),
          subtitle = paste("Metric:", tools::toTitleCase(mode), "| Relative to Baseline (CS00)"),
          theme = ggplot2::theme(plot.title = ggplot2::element_text(face="bold", size=16, hjust=0.5))
        ) +
        patchwork::plot_layout(guides = "collect") &
        ggplot2::theme(legend.position = "bottom")

      # Save
      file_name <- paste0(samp, "_CS_Retention_", mode, ".pdf")
      save_path <- file.path(output_dir, file_name)

      # Standard A4 Landscape
      ggplot2::ggsave(save_path, final_layout, width = kariocas_dims$width, height = kariocas_dims$height)
    }
  }

  message("\nSUCCESS: Retention analysis completed. Files saved in: ", output_dir)
  return(invisible(NULL))
}
