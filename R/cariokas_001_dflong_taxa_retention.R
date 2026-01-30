#' Run Confidence Score Retention Analysis (Step 001)
#'
#' Executes read retention analysis based on Confidence Score (Kraken/Bracken).
#' Adapted for the new Tidy Data format.
#'
#' @param project_dir Root path of the project.
#' @param import_script_path Deprecated.
#'
#' @export
#' @importFrom readr read_rds
#' @importFrom dplyr rename mutate filter select group_by summarise n distinct left_join bind_rows pull case_when any_of
#' @importFrom tidyr pivot_longer replace_na
#' @importFrom stringr str_detect str_remove str_extract
#' @importFrom ggplot2 ggplot aes geom_line geom_point scale_color_manual scale_shape_manual scale_linetype_manual scale_y_log10 scale_x_continuous labs guides guide_legend element_line coord_cartesian annotate ggsave theme_classic theme element_text element_blank margin
#' @importFrom patchwork plot_layout plot_annotation wrap_plots
#' @import patchwork

taxa_retention <- function(project_dir, import_script_path = NULL) {

  # 1. SETUP & DIRECTORIES
  input_dir   <- file.path(project_dir, "000_mpa_original")
  output_dir  <- file.path(project_dir, "001_cs_retention")

  if (!dir.exists(input_dir)) stop("Input directory not found: ", input_dir)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # 2. LOAD DATA (New Tidy Format)
  rds_file <- file.path(input_dir, "combined_mpas_matrix_taxon_level.rds")

  if (!file.exists(rds_file)) {
    rds_old <- file.path(input_dir, "000_unified_MPA_matrix.rds")
    if(file.exists(rds_old)) { rds_file <- rds_old }
    else { stop("Unified RDS file not found.") }
  }

  df_long <- load_and_process_mpa(rds_file)

  # 3. PRE-PROCESSING (Create Dynamic Taxon Column)
  df_proc <- df_long %>%
    dplyr::mutate(
      Taxon = dplyr::case_when(
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

  # 4. ANALYSIS LOOP
  SAMPLES <- unique(df_proc$sample)
  DOMAINS <- c("Bacteria", "Archaea", "Viruses", "Eukaryota")

  file_suffix_map <- c("Phylum"="Phyla", "Class"="Classes", "Order"="Orders",
                       "Family"="Families", "Genus"="Genera", "Species"="Species")
  TAX_LEVELS <- c(NA, "Phylum", "Class", "Order", "Family", "Genus", "Species")

  # --- STYLES ---
  orig_colors <- c("Phylum"="#000000", "Class"="#E69F00", "Order"="#56B4E9",
                   "Family"="#009E73", "Genus"="#F0E442", "Species"="#D55E00")
  orig_shapes <- c("Phylum"=15, "Class"=16, "Order"=17, "Family"=18, "Genus"=25, "Species"=8)

  spec_colors <- c("Total Reads"="#000000", "Level Reads"="#3C5488", "Level Taxa"="#E64B35")
  spec_linetypes <- c("Total Reads"="longdash", "Level Reads"="longdash", "Level Taxa"="solid")
  spec_shapes <- c("Total Reads"=15, "Level Reads"=17, "Level Taxa"=16)

  cat("\n=== Starting Retention Analysis (Step 001) ===\n")

  for (samp in SAMPLES) {
    cat("  Processing Sample:", samp, "\n")
    samp_out_dir <- file.path(output_dir, samp)
    if (!dir.exists(samp_out_dir)) dir.create(samp_out_dir)

    df_samp <- df_proc %>% dplyr::filter(sample == samp)

    # Baselines (CS == 0)
    df_baseline_domain <- df_samp %>%
      dplyr::filter(CS == 0, Rank == "Domain") %>%
      dplyr::select(Domain, Base_Domain_Reads = Counts)

    df_baseline_levels <- df_samp %>%
      dplyr::filter(CS == 0) %>%
      dplyr::group_by(Domain, Rank) %>%
      dplyr::summarise(
        Base_Lvl_Reads = sum(Counts),
        Base_Lvl_Taxa = dplyr::n_distinct(Taxon[Counts > 0]),
        .groups = "drop"
      )

    # Calculate Metrics
    df_stats <- df_samp %>%
      dplyr::group_by(Domain, CS, Rank) %>%
      dplyr::summarise(
        Curr_Reads = sum(Counts),
        Curr_Taxa = dplyr::n_distinct(Taxon[Counts > 0]),
        .groups = "drop"
      ) %>%
      dplyr::left_join(df_baseline_levels, by = c("Domain", "Rank")) %>%
      dplyr::left_join(df_baseline_domain, by = "Domain") %>%
      dplyr::mutate(
        Ret_Lvl_Reads = (Curr_Reads / Base_Lvl_Reads) * 100,
        Ret_Lvl_Taxa  = (Curr_Taxa / Base_Lvl_Taxa) * 100
      ) %>%
      dplyr::mutate(
        Ret_Lvl_Reads = ifelse(Ret_Lvl_Reads == 0, NA, Ret_Lvl_Reads),
        Ret_Lvl_Taxa  = ifelse(Ret_Lvl_Taxa == 0, NA, Ret_Lvl_Taxa)
      )

    df_domain_trend <- df_stats %>%
      dplyr::filter(Rank == "Domain") %>%
      dplyr::select(Domain, CS, Ret_Total_Reads_Trend = Ret_Lvl_Reads)

    # PLOTTING
    for (lvl in TAX_LEVELS) {
      plot_list <- list()
      if (is.na(lvl)) { ANALYSIS_LEVEL <- NULL; mode_tag <- "All_Levels" }
      else { ANALYSIS_LEVEL <- lvl; mode_tag <- file_suffix_map[[lvl]] }

      for (dom in DOMAINS) {
        base_dom_val <- df_baseline_domain$Base_Domain_Reads[df_baseline_domain$Domain == dom]
        if(length(base_dom_val)==0) base_dom_val <- 0
        reads_str <- format(base_dom_val, big.mark=".", decimal.mark=",")

        get_cnt <- function(l) {
          val <- df_baseline_levels$Base_Lvl_Taxa[df_baseline_levels$Domain == dom & df_baseline_levels$Rank == l]
          if(length(val)==0) 0 else val
        }

        # --- MODE: ALL LEVELS ---
        if (is.null(ANALYSIS_LEVEL)) {
          subtitle_stats <- paste0("Reads: ", reads_str, " | P: ", get_cnt("Phylum"), " | C: ", get_cnt("Class"),
                                   " | O: ", get_cnt("Order"), " | F: ", get_cnt("Family"),
                                   " | G: ", get_cnt("Genus"), " | S: ", get_cnt("Species"))

          dat_plot <- df_stats %>% dplyr::filter(Domain == dom, Rank %in% names(orig_colors))
          dat_plot$Rank <- factor(dat_plot$Rank, levels = c("Phylum", "Class", "Order", "Family", "Genus", "Species"))

          p <- ggplot2::ggplot(dat_plot, ggplot2::aes(x = CS, y = Ret_Lvl_Taxa, color = Rank, shape = Rank)) +
            ggplot2::geom_line(linewidth = 0.8) + ggplot2::geom_point(size = 2.5) +
            ggplot2::scale_color_manual(values = orig_colors) +
            ggplot2::scale_shape_manual(values = orig_shapes) +
            ggplot2::guides(color = ggplot2::guide_legend(nrow = 1), shape = ggplot2::guide_legend(nrow = 1)) +
            ggplot2::labs(subtitle = subtitle_stats)

        } else {
          # --- MODE: SPECIFIC ---
          curr_reads_base <- df_baseline_levels$Base_Lvl_Reads[df_baseline_levels$Domain == dom & df_baseline_levels$Rank == lvl]
          curr_reads_str <- format(curr_reads_base, big.mark=".", decimal.mark=",")
          subtitle_stats <- paste0(reads_str, " Total Reads; ", curr_reads_str, " assigned to ", lvl, " | ", get_cnt(lvl), " ", mode_tag)

          dat_lvl <- df_stats %>% dplyr::filter(Domain == dom, Rank == lvl) %>%
            dplyr::left_join(df_domain_trend, by = c("Domain", "CS")) %>%
            dplyr::select(CS, Ret_Total_Reads_Trend, Ret_Lvl_Reads, Ret_Lvl_Taxa) %>%
            tidyr::pivot_longer(cols = -CS, names_to = "Metric", values_to = "Pct")

          new_levels <- c("Ret_Lvl_Taxa", "Ret_Total_Reads_Trend", "Ret_Lvl_Reads")
          new_labels <- c(lvl, "Total Reads", paste(lvl, "Reads"))
          dat_lvl$Metric <- factor(dat_lvl$Metric, levels = new_levels, labels = new_labels)

          p <- ggplot2::ggplot(dat_lvl, ggplot2::aes(x = CS, y = Pct, group = Metric)) +
            ggplot2::geom_line(ggplot2::aes(color = Metric, linetype = Metric), linewidth = 1) +
            ggplot2::geom_point(ggplot2::aes(color = Metric, shape = Metric), size = 3) +
            ggplot2::scale_color_manual(values = spec_colors) +
            ggplot2::scale_linetype_manual(values = spec_linetypes) +
            ggplot2::scale_shape_manual(values = spec_shapes) +
            ggplot2::guides(color = ggplot2::guide_legend(nrow = 1), shape = ggplot2::guide_legend(nrow = 1), linetype = ggplot2::guide_legend(nrow = 1)) +
            ggplot2::labs(subtitle = subtitle_stats)
        }

        # --- COMMON STYLE ---
        p <- p + ggplot2::labs(title = dom, x = "Confidence Score (%)", y = "% Retained\n(axis scaled to log10)") +
          ggplot2::scale_x_continuous(breaks = seq(0, 100, 20), limits = c(0, 100)) +
          ggplot2::scale_y_log10(limits = c(0.01, 105), breaks = c(0.01, 0.1, 1, 10, 100), labels = c("0.01", "0.1", "1", "10", "100")) +
          ggplot2::theme_classic(base_size = 12) +
          ggplot2::theme(plot.title = ggplot2::element_text(face="bold", size=14, hjust=0.5),
                         plot.subtitle = ggplot2::element_text(hjust=0.5, size=10, color="grey30"),
                         panel.grid.major.y = ggplot2::element_line(color="grey90", linetype="dotted"),
                         legend.position="bottom", legend.title=ggplot2::element_blank()) +
          ggplot2::coord_cartesian(clip="off")
        plot_list[[dom]] <- p
      }
      final_layout <- (plot_list[["Bacteria"]] | plot_list[["Archaea"]]) / (plot_list[["Eukaryota"]] | plot_list[["Viruses"]]) +
        patchwork::plot_annotation(title = paste(samp, "- CS Retention Analysis"), subtitle = paste("Retention relative to Baseline (CS00) | Mode:", mode_tag),
                                   theme = ggplot2::theme(plot.title = ggplot2::element_text(face="bold", size=16, hjust=0.5))) +
        patchwork::plot_layout(guides = "collect") & ggplot2::theme(legend.position = "bottom")

      file_name <- paste0(samp, "_CS_Retention_", mode_tag, ".pdf")
      ggplot2::ggsave(file.path(samp_out_dir, file_name), final_layout, width = 11.69, height = 8.27, units = "in")
    }
  }
  cat("\nSUCCESS: Retention analysis completed.\n")
  return(invisible(NULL))
}
