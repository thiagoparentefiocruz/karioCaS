#' Generate Read Cutoff Saturation Analysis (Step 003)
#'
#' Performs a saturation analysis by progressively filtering taxa with low read counts.
#'
#' @param project_dir Path to the project root.
#' @param analysis_level Taxonomic rank to analyze (default: "Species").
#' @param cutoff_mode Visualization mode: "Extinction" (Full 0-100%), "Stabilization", or "Both".
#'
#' @export
#' @importFrom readr read_rds
#' @importFrom dplyr filter mutate select group_by summarise arrange pull n distinct case_when left_join bind_rows
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_line geom_point scale_color_manual scale_shape_manual scale_linetype_manual scale_y_continuous labs coord_cartesian guides guide_legend element_text
#' @importFrom patchwork plot_layout plot_annotation
#' @importFrom scales percent label_number
#' @importFrom ggtext element_markdown

reads_per_taxa <- function(project_dir,
                           analysis_level = "Species",
                           cutoff_mode = "Both") {

  # ==============================================================================
  # 1. SETUP & DIRECTORIES
  # ==============================================================================
  input_dir  <- file.path(project_dir, "000_karioCaS_input_matrix")
  output_dir <- file.path(project_dir, "003_cutoffs")

  if (!dir.exists(input_dir)) {
    stop("CRITICAL ERROR: Input directory not found: ", input_dir)
  }

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  rds_file <- file.path(input_dir, "karioCaS_input_matrix.rds")
  if (!file.exists(rds_file)) stop("Input RDS file not found: ", rds_file)

  # ==============================================================================
  # 2. DATA LOADING & PREP
  # ==============================================================================
  message(">>> Loading Tidy Data from: ", rds_file)
  df_long <- readr::read_rds(rds_file)

  df_proc <- df_long %>%
    dplyr::filter(Rank == analysis_level) %>%
    dplyr::mutate(Taxon_Name = get(analysis_level))

  SAMPLES <- unique(df_proc$sample)
  DOMAINS <- names(kariocas_colors$domains)
  CS_LIST <- unique(df_proc$CS)

  # Template for calculation (Dense to make smooth curves)
  base_steps <- c(1, 3, 5)
  multipliers <- c(1, 10, 100, 1000, 10000, 100000, 1000000)
  cutoff_template <- sort(unique(c(0, as.vector(outer(base_steps, multipliers, "*")))))

  message(">>> Starting Cutoff Analysis for ", length(SAMPLES), " samples.")
  message(">>> Level: ", analysis_level, " | Mode: ", cutoff_mode)

  # ==============================================================================
  # 3. ANALYSIS LOOP
  # ==============================================================================
  for (samp in SAMPLES) {
    message("  Processing Sample: ", samp)

    for (cs in CS_LIST) {

      plot_list <- list()
      modes_to_run <- if (cutoff_mode == "Both") c("Extinction", "Stabilization") else cutoff_mode

      df_curr <- df_proc %>% dplyr::filter(sample == samp, CS == cs)

      if (nrow(df_curr) == 0) next

      for (mode in modes_to_run) {

        for (dom in DOMAINS) {

          # 3.1 CALCULATION ENGINE ---------------------------------------------
          df_dom <- df_curr %>% dplyr::filter(Domain == dom, Counts > 0)

          total_reads_dom <- sum(df_dom$Counts)
          total_taxa_dom  <- nrow(df_dom)

          if (total_reads_dom == 0) {
            plot_list[[dom]] <- plot_kariocas_empty(
              title_text = dom,
              subtitle_text = "No reads detected",
              x_label = kariocas_labels$x_log10_reads,
              y_label = "**% Retained**"
            )
            next
          }

          # Calculate Curve Points
          max_val <- max(df_dom$Counts)
          relevant_cutoffs <- cutoff_template[cutoff_template <= max_val]
          # Add one point beyond max to ensure drop to zero is visible
          if(length(relevant_cutoffs) > 0) {
            next_cut <- cutoff_template[which(cutoff_template > max_val)[1]]
            if(!is.na(next_cut)) relevant_cutoffs <- c(relevant_cutoffs, next_cut)
          }

          stats_list <- lapply(relevant_cutoffs, function(k) {
            survivors <- df_dom %>% dplyr::filter(Counts >= k)
            data.frame(
              Cutoff = k,
              Ret_Reads_Pct = sum(survivors$Counts) / total_reads_dom,
              Ret_Taxa_Pct  = nrow(survivors) / total_taxa_dom
            )
          })

          df_stats <- do.call(rbind, stats_list)

          # 3.2 PLOT PREP ------------------------------------------------------
          df_plot <- df_stats %>%
            tidyr::pivot_longer(
              cols = c(Ret_Reads_Pct, Ret_Taxa_Pct),
              names_to = "Metric",
              values_to = "Pct"
            ) %>%
            dplyr::mutate(Metric_Key = dplyr::case_when(
              Metric == "Ret_Reads_Pct" ~ "Level Reads",
              Metric == "Ret_Taxa_Pct"  ~ "Level Taxa"
            ))

          df_plot$Metric_Key <- factor(df_plot$Metric_Key, levels = c("Level Taxa", "Level Reads"))
          legend_labels <- c(analysis_level, "Reads")

          # --- DYNAMIC LIMIT LOGIC (Fixing Issue 1 & 2) ---
          # 1. Determine the "Visual Ceiling": Max Data OR 10k, whichever is smaller
          # This ensures small datasets don't have empty space up to 10k.
          # And large datasets get cut at 10k.
          if (mode == "Extinction") {
            display_limit <- min(max_val, 10000)

            # 2. Filter DATA to stop the curve at the limit (Prevents bleeding)
            df_plot <- df_plot %>% dplyr::filter(Cutoff <= display_limit)
          } else {
            # For Stabilization, we usually show full range or let ggplot decide
            display_limit <- max_val
          }

          # --- AXIS BREAKS LOGIC (Fixing Issue 3) ---
          # Generate clean Powers of 10 breaks: 1, 10, 100, 1000, 10000
          # We intersect this with the range [1, display_limit]
          all_powers <- c(1, 10, 100, 1000, 10000, 100000, 1000000)
          clean_breaks <- all_powers[all_powers <= display_limit]

          # Ensure we have at least min and max if breaks are too sparse?
          # Actually, powers of 10 are usually enough and cleanest.

          subtitle_stats <- paste0(label_kariocas_auto(total_reads_dom), " Reads | ",
                                   label_kariocas_auto(total_taxa_dom), " ", analysis_level)

          # 3.3 PLOT CONSTRUCTION ----------------------------------------------
          p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = Cutoff, y = Pct, group = Metric_Key)) +
            ggplot2::geom_line(ggplot2::aes(color = Metric_Key, linetype = Metric_Key), linewidth = 1) +
            ggplot2::geom_point(ggplot2::aes(color = Metric_Key, shape = Metric_Key), size = 3) +

            # X Axis: Powers of 10 breaks + K notation
            scale_x_kariocas_log10(
              breaks = clean_breaks,
              labels = label_k_number # Use K notation (1K, 10K)
            ) +

            # Y Axis: 0-100 Integers
            ggplot2::scale_y_continuous(
              breaks = seq(0, 1, 0.25),
              labels = function(x) x * 100
            ) +

            # Manual Styles
            ggplot2::scale_color_manual(values = kariocas_colors$special, labels = legend_labels) +
            ggplot2::scale_shape_manual(values = kariocas_shapes, labels = legend_labels) +
            ggplot2::scale_linetype_manual(values = kariocas_linetypes, labels = legend_labels) +

            ggplot2::labs(
              title = dom,
              subtitle = subtitle_stats,
              x = kariocas_labels$x_log10_reads,
              y = "**% Retained**"
            ) +

            theme_kariocas() +

            ggplot2::guides(
              color = ggplot2::guide_legend(nrow = 1, title = NULL),
              shape = ggplot2::guide_legend(nrow = 1, title = NULL),
              linetype = ggplot2::guide_legend(nrow = 1, title = NULL)
            ) +

            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5))

          # 3.4 APPLY LIMITS ---------------------------------------------------
          if (mode == "Extinction") {
            # Force X axis to stop exactly at display_limit (dynamic)
            p <- p + ggplot2::coord_cartesian(xlim = c(NA, display_limit), ylim = c(0, 1.05), clip = "off")
          } else {
            # Stabilization: Full X, Y full 0-100 (per request)
            p <- p + ggplot2::coord_cartesian(ylim = c(0, 1.05), clip = "off")
          }

          plot_list[[dom]] <- p

        } # End Domain Loop

        # 3.5 ASSEMBLE PANEL ---------------------------------------------------
        final_layout <- (plot_list[["Bacteria"]] | plot_list[["Archaea"]]) /
          (plot_list[["Eukaryota"]] | plot_list[["Viruses"]]) +
          patchwork::plot_annotation(
            title = paste(samp, "- CS", sprintf("%02d", cs), "| Saturation:", mode),
            subtitle = paste("Retention of Reads vs", analysis_level),
            theme = ggplot2::theme(
              plot.title = ggplot2::element_text(face = "bold", size = 16, hjust = 0.5),
              plot.subtitle = ggplot2::element_text(size = 12, hjust = 0.5, color = "grey30")
            )
          ) +
          patchwork::plot_layout(guides = "collect") &
          ggplot2::theme(legend.position = "bottom")

        file_name <- paste0(samp, "_CS", sprintf("%02d", cs), "_Cutoff_", analysis_level, "_", mode, ".pdf")
        save_path <- file.path(output_dir, file_name)

        ggplot2::ggsave(save_path, final_layout, width = kariocas_dims$width, height = kariocas_dims$height)

      } # End Mode Loop
    } # End CS Loop
  } # End Sample Loop

  message("\nSUCCESS: Cutoff analysis completed. Files saved in: ", output_dir)
  return(invisible(NULL))
}
