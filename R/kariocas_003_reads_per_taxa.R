#' Generate Read Cutoff Saturation Analysis (Step 003)
#'
#' Performs a saturation analysis by progressively filtering taxa with low read counts.
#' Generates two views: "Saturation" (Log scale overview) and "Rare_Taxa" (Linear scale focus on low counts).
#'
#' @param project_dir Path to the project root.
#' @param analysis_level Taxonomic rank to analyze (default: "Species").
#' @param x_max_bac Integer. Max X-axis value for Bacteria in Rare Taxa plot (default: 10).
#' @param x_max_arc Integer. Max X-axis value for Archaea in Rare Taxa plot (default: 10).
#' @param x_max_euk Integer. Max X-axis value for Eukaryota in Rare Taxa plot (default: 10).
#' @param x_max_vir Integer. Max X-axis value for Viruses in Rare Taxa plot (default: 10).
#'
#' @return Generates PDF plots saved in the project directory and returns a data frame invisibly.
#' @export
#' @importFrom dplyr filter mutate select group_by summarise arrange pull n distinct case_when left_join bind_rows
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_line geom_point scale_color_manual scale_shape_manual scale_linetype_manual scale_y_continuous scale_x_continuous labs coord_cartesian guides guide_legend element_text
#' @importFrom patchwork plot_layout plot_annotation
#' @importFrom scales percent label_number breaks_pretty
#' @importFrom ggtext element_markdown
#' @examples
#' # Get the path to the included toy dataset
#' toy_project <- system.file("extdata", "your_project_name", package = "karioCaS")
#'
#' # Example 1: Basic usage (defaults to Species and x_max = 10 for all domains)
#' # reads_per_taxa(
#' #   project_dir = toy_project
#' # )
#' 
#' # Example 2: Analyzing Genus level and adjusting the X-axis zoom for Rare Taxa
#' # reads_per_taxa(
#' #   project_dir = toy_project,
#' #   analysis_level = "Genus",
#' #   x_max_bac = 50,
#' #   x_max_vir = 5
#' # )

reads_per_taxa <- function(project_dir,
                           analysis_level = "Species",
                           x_max_bac = 10,
                           x_max_arc = 10,
                           x_max_euk = 10,
                           x_max_vir = 10) {

  # ==============================================================================
  # 1. SETUP & LOGGING
  # ==============================================================================
  output_dir <- file.path(project_dir, "003_cutoffs")
  log_dir    <- file.path(project_dir, "logs")

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)

  log_file <- file.path(log_dir, "log_003_cutoff_analysis.txt")

  log_msg <- function(...) {
    msg <- paste0(...)
    message(msg)
    cat(paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", msg, "\n"),
        file = log_file, append = TRUE)
  }

  cat("====================================================\n", file = log_file)
  cat("LOG: 003_CUTOFF_ANALYSIS (Saturation + Rare Taxa)\n", file = log_file, append = TRUE)
  cat("PROJECT DIR: ", project_dir, "\n", file = log_file, append = TRUE)
  cat("ANALYSIS LEVEL: ", analysis_level, "\n", file = log_file, append = TRUE)
  cat("RARE TAXA MAX: B=", x_max_bac, "|A=", x_max_arc, "|E=", x_max_euk, "|V=", x_max_vir, "\n", file = log_file, append = TRUE)
  cat("====================================================\n", file = log_file, append = TRUE)

  # ==============================================================================
  # 2. DATA LOADING & PREP
  # ==============================================================================
  log_msg(">>> Loading Data (Auto-detected format)...")
  df_long <- .get_tidy_data(project_dir)

  # Filter for the specific level requested
  df_proc <- df_long %>%
    dplyr::filter(Rank == analysis_level)

  if (nrow(df_proc) == 0) {
    log_msg("CRITICAL ERROR: No data found for Rank: ", analysis_level)
    stop("No data for specified rank.")
  }

  SAMPLES <- unique(df_proc$sample)
  DOMAINS <- names(kariocas_colors$domains)
  CS_LIST <- unique(df_proc$CS)

  # Template for Saturation calculation (Log Scale)
  base_steps <- c(1, 3, 5)
  multipliers <- c(1, 10, 100, 1000, 10000, 100000, 1000000)
  cutoff_template_log <- sort(unique(c(0, as.vector(outer(base_steps, multipliers, "*")))))

  # Map for Rare Taxa Limits
  rare_limits <- list(
    "Bacteria"  = x_max_bac,
    "Archaea"   = x_max_arc,
    "Eukaryota" = x_max_euk,
    "Viruses"   = x_max_vir
  )

  log_msg(">>> Starting Cutoff Analysis for ", length(SAMPLES), " samples.")

  # ==============================================================================
  # 3. ANALYSIS LOOP
  # ==============================================================================
  for (samp in SAMPLES) {
    log_msg("------------------------------------------------")
    log_msg("  Processing Sample: ", samp)

    for (cs in CS_LIST) {

      plot_list <- list()
      df_curr <- df_proc %>% dplyr::filter(sample == samp, CS == cs)

      if (nrow(df_curr) == 0) {
        log_msg("    Skipping CS", cs, ": No data.")
        next
      }

      # Run both modes: Saturation (Standard) and Rare_Taxa (New)
      modes_to_run <- c("Saturation", "Rare_Taxa")

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
              x_label = if(mode=="Saturation") kariocas_labels$x_log10_reads else "**Reads**",
              y_label = "**% Retained**"
            )
            next
          }

          # Define Cutoff Points based on Mode
          max_val <- max(df_dom$Counts)

          if (mode == "Saturation") {
            # Standard Log Logic
            relevant_cutoffs <- cutoff_template_log[cutoff_template_log <= max_val]
            # Add one point beyond max to ensure drop to zero is visible
            if(length(relevant_cutoffs) > 0) {
              next_cut <- cutoff_template_log[which(cutoff_template_log > max_val)[1]]
              if(!is.na(next_cut)) relevant_cutoffs <- c(relevant_cutoffs, next_cut)
            }
          } else {
            # Rare Taxa (Linear 1 to X_MAX)
            dom_limit <- rare_limits[[dom]]
            relevant_cutoffs <- seq(1, dom_limit, by = 1)
          }

          # Calculate Stats
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

          subtitle_stats <- paste0(label_kariocas_auto(total_reads_dom), " Reads | ",
                                   label_kariocas_auto(total_taxa_dom), " ", analysis_level)

          # Style Setup
          shapes_vec     <- get_kariocas_shapes("ranks")
          spec_colors    <- get_kariocas_colors("special")
          spec_linetypes <- get_kariocas_linetypes()
          labels_vec     <- get_kariocas_labels()

          # 3.3 PLOT CONSTRUCTION ----------------------------------------------
          p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = Cutoff, y = Pct, group = Metric_Key)) +
            ggplot2::geom_line(ggplot2::aes(color = Metric_Key, linetype = Metric_Key), linewidth = 1) +
            ggplot2::geom_point(ggplot2::aes(color = Metric_Key, shape = Metric_Key), size = 3) +

            # Y Axis: 0-100 Integers
            ggplot2::scale_y_continuous(
              breaks = seq(0, 1, 0.25),
              labels = function(x) x * 100
            ) +

            # Manual Styles
            ggplot2::scale_color_manual(values = spec_colors, labels = legend_labels) +
            ggplot2::scale_shape_manual(values = shapes_vec, labels = legend_labels) +
            ggplot2::scale_linetype_manual(values = kariocas_linetypes, labels = legend_labels) +

            theme_kariocas() +

            ggplot2::guides(
              color = ggplot2::guide_legend(nrow = 1, title = NULL),
              shape = ggplot2::guide_legend(nrow = 1, title = NULL),
              linetype = ggplot2::guide_legend(nrow = 1, title = NULL)
            ) +

            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5))

          # 3.4 MODE SPECIFIC SCALES & LIMITS ----------------------------------

          if (mode == "Saturation") {
            # Log Scale logic
            display_limit <- min(max_val, 10000)

            # Breaks logic
            all_powers <- c(1, 10, 100, 1000, 10000, 100000, 1000000)
            clean_breaks <- all_powers[all_powers <= display_limit]

            # CLIPPING ON (Standard behavior to prevent leaks)
            p <- p +
              scale_x_kariocas_log10(breaks = clean_breaks, labels = label_k_number) +
              ggplot2::coord_cartesian(xlim = c(NA, display_limit), ylim = c(0, 1.05), clip = "on") +
              ggplot2::labs(title = dom, subtitle = subtitle_stats, x = kariocas_labels$x_log10_reads, y = "**% Retained**")

          } else {
            # Rare Taxa Logic (Linear Scale)
            dom_limit <- rare_limits[[dom]]

            # Automatic Integer Ticks
            if(dom_limit <= 20) {
              x_breaks <- seq(1, dom_limit)
            } else {
              x_breaks <- scales::breaks_pretty(n = 10)(c(1, dom_limit))
              x_breaks <- unique(round(x_breaks)) # Ensure integers
              x_breaks <- x_breaks[x_breaks >= 1 & x_breaks <= dom_limit]
            }

            # CLIPPING ON
            p <- p +
              ggplot2::scale_x_continuous(breaks = x_breaks) +
              ggplot2::coord_cartesian(xlim = c(1, dom_limit), ylim = c(0, 1.05), clip = "on") +
              ggplot2::labs(title = dom, subtitle = subtitle_stats, x = "**Reads**", y = "**% Retained**")
          }

          plot_list[[dom]] <- p

        } # End Domain Loop

        # 3.5 ASSEMBLE PANEL ---------------------------------------------------
        final_layout <- (plot_list[["Bacteria"]] | plot_list[["Archaea"]]) /
          (plot_list[["Eukaryota"]] | plot_list[["Viruses"]]) +
          patchwork::plot_annotation(
            title = paste(samp, "- CS", sprintf("%02d", cs), "| Analysis:", mode),
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
        log_msg("    -> Saved [", mode, "]: ", file_name)

      } # End Mode Loop
    } # End CS Loop
  } # End Sample Loop

  log_msg("SUCCESS: Cutoff analysis completed.")
  return(invisible(NULL))
}
