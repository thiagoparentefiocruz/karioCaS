#' Run Confidence Score Retention Analysis (Step 001)
#'
#' Executes read retention analysis based on Confidence Score (Kraken/Bracken).
#'
#' @param project_dir Root path of the project.
#'
#' @export
#' @importFrom readr read_rds
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
  input_dir  <- file.path(project_dir, "000_karioCaS_input_matrix")
  output_dir <- file.path(project_dir, "001_cs_retention")

  if (!dir.exists(input_dir)) {
    stop("CRITICAL ERROR: Input directory not found: ", input_dir,
         "\nPlease ensure step 000 (import) was run successfully.")
  }

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  rds_file <- file.path(input_dir, "karioCaS_input_matrix.rds")
  if (!file.exists(rds_file)) stop("Input RDS file not found: ", rds_file)

  # ==============================================================================
  # 2. DATA LOADING
  # ==============================================================================
  message(">>> Loading Tidy Data from: ", rds_file)
  df_long <- readr::read_rds(rds_file)

  df_proc <- df_long %>%
    dplyr::mutate(
      Taxon_Name = dplyr::case_when(
        Rank == "Domain"  ~ Domain,
        Rank == "Kingdom" ~ Kingdom,
        Rank == "Phylum"  ~ Phylum,
        Rank == "Class"   ~ Class,
        Rank == "Order"   ~ Order,
        Rank == "Family"  ~ Family,
        Rank == "Genus"   ~ Genus,
        Rank == "Species" ~ Species,
        TRUE ~ NA_character_
      )
    )

  SAMPLES <- unique(df_proc$sample)
  DOMAINS <- names(kariocas_colors$domains)

  rank_plural_map <- c("Phylum"="Phyla", "Class"="Classes", "Order"="Orders",
                       "Family"="Families", "Genus"="Genera", "Species"="Species")

  TARGET_RANKS <- c(NA, "Phylum", "Class", "Order", "Family", "Genus", "Species")

  message(">>> Starting Retention Analysis for ", length(SAMPLES), " samples.")

  # ==============================================================================
  # 3. CORE LOGIC LOOP
  # ==============================================================================
  for (samp in SAMPLES) {
    message("  Processing Sample: ", samp)
    samp_out_dir <- file.path(output_dir, samp)
    if (!dir.exists(samp_out_dir)) dir.create(samp_out_dir)

    df_samp <- df_proc %>% dplyr::filter(sample == samp)

    # 3.1 CALCULATE BASELINES (CS == 0) ----------------------------------------
    df_baseline_domain <- df_samp %>%
      dplyr::filter(CS == 0, Rank == "Domain") %>%
      dplyr::group_by(Domain) %>%
      dplyr::summarise(Base_Domain_Reads = sum(Counts), .groups = "drop")

    df_baseline_ranks <- df_samp %>%
      dplyr::filter(CS == 0) %>%
      dplyr::group_by(Domain, Rank) %>%
      dplyr::summarise(
        Base_Rank_Reads = sum(Counts),
        Base_Rank_Taxa  = dplyr::n_distinct(Taxon_Name[Counts > 0]),
        .groups = "drop"
      )

    # 3.2 CALCULATE METRICS ----------------------------------------------------
    df_stats <- df_samp %>%
      dplyr::group_by(Domain, CS, Rank) %>%
      dplyr::summarise(
        Curr_Reads = sum(Counts),
        Curr_Taxa  = dplyr::n_distinct(Taxon_Name[Counts > 0]),
        .groups = "drop"
      ) %>%
      dplyr::left_join(df_baseline_domain, by = "Domain") %>%
      dplyr::left_join(df_baseline_ranks, by = c("Domain", "Rank")) %>%
      dplyr::mutate(
        Ret_Reads_Pct = (Curr_Reads / Base_Rank_Reads) * 100,
        Ret_Taxa_Pct  = (Curr_Taxa / Base_Rank_Taxa) * 100
      ) %>%
      dplyr::mutate(
        Ret_Reads_Pct = ifelse(is.nan(Ret_Reads_Pct) | Ret_Reads_Pct == 0, NA, Ret_Reads_Pct),
        Ret_Taxa_Pct  = ifelse(is.nan(Ret_Taxa_Pct) | Ret_Taxa_Pct == 0, NA, Ret_Taxa_Pct)
      )

    # Global Domain Trend (using Rank==Domain logic)
    df_domain_trend <- df_samp %>%
      dplyr::filter(Rank == "Domain") %>%
      dplyr::group_by(Domain, CS) %>%
      dplyr::summarise(Curr_Domain_Reads = sum(Counts), .groups="drop") %>%
      dplyr::left_join(df_baseline_domain, by = "Domain") %>%
      dplyr::mutate(
        Ret_Total_Reads_Trend = (Curr_Domain_Reads / Base_Domain_Reads) * 100,
        Ret_Total_Reads_Trend = ifelse(Ret_Total_Reads_Trend == 0, NA, Ret_Total_Reads_Trend)
      ) %>%
      dplyr::select(Domain, CS, Ret_Total_Reads_Trend)

    # ==========================================================================
    # 4. VISUALIZATION LOOP
    # ==========================================================================
    for (lvl in TARGET_RANKS) {

      plot_list <- list()
      is_all_levels <- is.na(lvl)
      mode_tag      <- if(is_all_levels) "All_Levels" else rank_plural_map[[lvl]]

      for (dom in DOMAINS) {

        # --- SUBTITLE STATS ---
        base_dom_reads <- df_baseline_domain$Base_Domain_Reads[df_baseline_domain$Domain == dom]
        if(length(base_dom_reads) == 0 || is.na(base_dom_reads)) base_dom_reads <- 0
        reads_str <- label_kariocas_auto(base_dom_reads)

        get_cnt <- function(r) {
          val <- df_baseline_ranks$Base_Rank_Taxa[df_baseline_ranks$Domain == dom & df_baseline_ranks$Rank == r]
          if(length(val) == 0) return(0) else return(val)
        }

        # Check for Data (using Domain trend as minimum requirement)
        has_domain_data <- any(df_domain_trend$Domain == dom & !is.na(df_domain_trend$Ret_Total_Reads_Trend))

        if (!has_domain_data || base_dom_reads == 0) {
          plot_list[[dom]] <- plot_kariocas_empty(
            title_text = dom,
            subtitle_text = "No reads detected at CS00",
            x_label = kariocas_labels$y_confidence,
            y_label = kariocas_labels$y_log10_retained
          )
          next
        }

        # --- PLOT CONSTRUCTION ---
        if (is_all_levels) {
          # ... (All Levels code remains same, it's safer) ...
          subtitle_stats <- paste0(
            "Reads: ", reads_str, " | P: ", get_cnt("Phylum"), " | C: ", get_cnt("Class"),
            " | O: ", get_cnt("Order"), " | F: ", get_cnt("Family"),
            " | G: ", get_cnt("Genus"), " | S: ", get_cnt("Species")
          )

          dat_plot <- df_stats %>%
            dplyr::filter(Domain == dom, Rank %in% names(kariocas_colors$ranks))
          dat_plot$Rank <- factor(dat_plot$Rank, levels = names(kariocas_colors$ranks))

          p <- ggplot2::ggplot(dat_plot, ggplot2::aes(x = CS, y = Ret_Taxa_Pct, color = Rank, shape = Rank)) +
            ggplot2::geom_line(linewidth = 0.8) +
            ggplot2::geom_point(size = 2.5) +
            ggplot2::scale_color_manual(values = kariocas_colors$ranks) +
            ggplot2::scale_shape_manual(values = kariocas_shapes) +
            ggplot2::guides(color = ggplot2::guide_legend(nrow = 1), shape = ggplot2::guide_legend(nrow = 1))

        } else {
          # >>> SPECIFIC RANK MODE WITH FIX <<<
          curr_reads_base <- df_baseline_ranks$Base_Rank_Reads[df_baseline_ranks$Domain == dom & df_baseline_ranks$Rank == lvl]
          curr_reads_str  <- label_kariocas_auto(curr_reads_base)

          subtitle_stats <- paste0(
            reads_str, " Total Reads; ", curr_reads_str, " assigned to ", lvl,
            " | ", get_cnt(lvl), " ", mode_tag
          )

          # 1. Get Domain Trend (Backbone)
          backbone <- df_domain_trend %>% dplyr::filter(Domain == dom)

          # 2. Get Specific Level Stats
          level_stats <- df_stats %>%
            dplyr::filter(Domain == dom, Rank == lvl) %>%
            dplyr::select(Domain, CS, Ret_Lvl_Reads = Ret_Reads_Pct, Ret_Lvl_Taxa = Ret_Taxa_Pct)

          # 3. FULL JOIN (Logic Fix: Keep backbone even if level dies)
          dat_lvl <- dplyr::full_join(backbone, level_stats, by = c("Domain", "CS")) %>%
            dplyr::select(CS, Ret_Total_Reads_Trend, Ret_Lvl_Reads, Ret_Lvl_Taxa) %>%
            tidyr::pivot_longer(cols = -CS, names_to = "Metric", values_to = "Pct")

          # 4. Styling Mappings
          dat_lvl <- dat_lvl %>%
            dplyr::mutate(Metric_Key = dplyr::case_when(
              Metric == "Ret_Lvl_Taxa" ~ "Level Taxa",
              Metric == "Ret_Total_Reads_Trend" ~ "Total Reads",
              Metric == "Ret_Lvl_Reads" ~ "Level Reads"
            ))

          dat_lvl$Metric_Key <- factor(dat_lvl$Metric_Key, levels = c("Level Taxa", "Total Reads", "Level Reads"))
          legend_labels <- c(lvl, "Total Reads", paste(lvl, "Reads"))

          p <- ggplot2::ggplot(dat_lvl, ggplot2::aes(x = CS, y = Pct, group = Metric_Key)) +
            ggplot2::geom_line(ggplot2::aes(color = Metric_Key, linetype = Metric_Key), linewidth = 1) +
            ggplot2::geom_point(ggplot2::aes(color = Metric_Key, shape = Metric_Key), size = 3) +

            ggplot2::scale_color_manual(values = kariocas_colors$special, labels = legend_labels) +
            ggplot2::scale_shape_manual(values = kariocas_shapes, labels = legend_labels) +
            ggplot2::scale_linetype_manual(values = kariocas_linetypes, labels = legend_labels) +

            ggplot2::guides(
              color = ggplot2::guide_legend(nrow = 1, title = NULL),
              shape = ggplot2::guide_legend(nrow = 1, title = NULL),
              linetype = ggplot2::guide_legend(nrow = 1, title = NULL)
            )
        }

        # --- COMMON STYLES ---
        p <- p +
          scale_y_kariocas_log10(limits = c(0.01, 105)) +
          ggplot2::scale_x_continuous(breaks = seq(0, 100, 20), limits = c(0, 100)) +
          ggplot2::labs(
            title = dom,
            subtitle = subtitle_stats,
            x = kariocas_labels$y_confidence,
            y = kariocas_labels$y_log10_retained
          ) +
          theme_kariocas() +
          ggplot2::coord_cartesian(clip = "off")

        plot_list[[dom]] <- p
      }

      # --- ASSEMBLE ---
      final_layout <- (plot_list[["Bacteria"]] | plot_list[["Archaea"]]) /
        (plot_list[["Eukaryota"]] | plot_list[["Viruses"]]) +
        patchwork::plot_annotation(
          title = paste(samp, "- CS Retention Analysis"),
          subtitle = paste("Retention relative to Baseline (CS00) | Mode:", mode_tag),
          theme = ggplot2::theme(plot.title = ggplot2::element_text(face="bold", size=16, hjust=0.5))
        ) +
        patchwork::plot_layout(guides = "collect") &
        ggplot2::theme(legend.position = "bottom")

      file_name <- paste0(samp, "_CS_Retention_", mode_tag, ".pdf")
      save_path <- file.path(samp_out_dir, file_name)

      ggplot2::ggsave(save_path, final_layout, width = kariocas_dims$width, height = kariocas_dims$height)

    }
  }

  message("\nSUCCESS: Retention analysis completed with ggtext styling.")
  return(invisible(NULL))
}
