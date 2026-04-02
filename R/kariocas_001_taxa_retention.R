#' Run Confidence Score Retention Analysis (Step 001)
#'
#' Executes read retention analysis based on Confidence Score (Kraken/Bracken).
#' Generates detailed logs including baseline stats for QC.
#'
#' @param project_dir Root path of the project.
#'
#' @export
#' @importFrom dplyr filter select group_by summarise mutate left_join arrange rename bind_rows pull distinct
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_line geom_point scale_color_manual scale_linetype_manual scale_shape_manual labs ggsave scale_y_continuous scale_x_continuous theme guides guide_legend
#' @importFrom patchwork plot_layout plot_annotation
#' @importFrom scales label_number
#' @importFrom ggtext element_markdown

taxa_retention <- function(project_dir) {

  # ==============================================================================
  # 1. SETUP & LOGGING
  # ==============================================================================
  output_dir <- file.path(project_dir, "001_taxa_retention")
  log_dir    <- file.path(project_dir, "logs")

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)

  log_file <- file.path(log_dir, "log_001_cs_retention.txt")

  log_msg <- function(...) {
    msg <- paste0(...)
    message(msg)
    cat(paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", msg, "\n"),
        file = log_file, append = TRUE)
  }

  cat("====================================================\n", file = log_file)
  cat("LOG: 001_CS_RETENTION_ANALYSIS\n", file = log_file, append = TRUE)
  cat("PROJECT DIR: ", project_dir, "\n", file = log_file, append = TRUE)
  cat("====================================================\n", file = log_file, append = TRUE)

  # ==============================================================================
  # 2. DATA LOADING
  # ==============================================================================
  log_msg(">>> Loading Data (Auto-detected format)...")
  df_long <- .get_tidy_data(project_dir)

  rank_levels_all <- c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

  df_proc <- df_long %>%
    dplyr::filter(Rank %in% rank_levels_all) %>%
    dplyr::mutate(Rank = factor(Rank, levels = rank_levels_all))

  if (nrow(df_proc) == 0) {
    log_msg("CRITICAL ERROR: No data found for specified ranks.")
    stop("No data found.")
  }

  SAMPLES <- unique(df_proc$sample)
  DOMAINS <- names(kariocas_colors$domains)

  fmt_num <- function(x) format(x, big.mark = ",", scientific = FALSE)

  log_msg(">>> Starting Retention Analysis for ", length(SAMPLES), " samples.")

  # ==============================================================================
  # 3. ANALYSIS LOOP
  # ==============================================================================
  for (samp in SAMPLES) {
    log_msg("------------------------------------------------")
    log_msg("  Processing Sample: ", samp)

    df_samp <- df_proc %>% dplyr::filter(sample == samp)

    # 3.1 PRE-CALCULATION
    domain_totals <- df_samp %>%
      dplyr::filter(Rank == "Domain", Lowest_Rank == "Domain") %>%
      dplyr::group_by(Domain, CS) %>%
      dplyr::summarise(Global_Reads = max(Counts), .groups = "drop")

    if (nrow(domain_totals) == 0) {
      domain_totals <- df_samp %>%
        dplyr::filter(Rank == "Domain") %>%
        dplyr::group_by(Domain, CS) %>%
        dplyr::summarise(Global_Reads = max(Counts), .groups = "drop")
    }

    # 3.2 DETAILED STATS
    stats_df <- df_samp %>%
      dplyr::group_by(Domain, Rank, CS) %>%
      dplyr::summarise(
        Rank_Reads = sum(Counts),
        Rank_Taxa  = dplyr::n_distinct(Taxon_Name),
        .groups = "drop"
      ) %>%
      dplyr::left_join(domain_totals, by = c("Domain", "CS"))

    # 3.3 BASELINE LOGGING (QC STATS)
    baseline_df <- stats_df %>%
      dplyr::filter(CS == 0) %>%
      dplyr::rename(
        Base_Global = Global_Reads,
        Base_Reads  = Rank_Reads,
        Base_Taxa   = Rank_Taxa
      ) %>%
      dplyr::select(Domain, Rank, Base_Global, Base_Reads, Base_Taxa)

    # --- ENHANCED LOGGING HERE ---
    # Log the baseline counts for Species (most important rank)
    qc_stats <- baseline_df %>% dplyr::filter(Rank == "Species")
    if (nrow(qc_stats) > 0) {
      log_msg("    > QC STATS (CS00 Baseline):")
      for(i in 1:nrow(qc_stats)) {
        log_msg("      - ", qc_stats$Domain[i], ": ",
                fmt_num(qc_stats$Base_Global[i]), " Total Reads | ",
                fmt_num(qc_stats$Base_Taxa[i]), " Distinct Species")
      }
    }

    df_calc <- stats_df %>%
      dplyr::left_join(baseline_df, by = c("Domain", "Rank")) %>%
      dplyr::mutate(
        Pct_Domain_Retained = ifelse(Base_Global > 0, (Global_Reads / Base_Global) * 100, 0),
        Pct_Rank_Reads      = ifelse(Base_Reads > 0, (Rank_Reads / Base_Reads) * 100, 0),
        Pct_Taxa            = ifelse(Base_Taxa > 0, (Rank_Taxa / Base_Taxa) * 100, 0)
      )

    # ==========================================================================
    # 4. PLOT TYPE A
    # ==========================================================================
    plot_ranks_a <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
    plot_list_all <- list()

    for (dom in DOMAINS) {
      df_dom <- df_calc %>% dplyr::filter(Domain == dom, Rank %in% plot_ranks_a)

      if (nrow(df_dom) == 0) {
        plot_list_all[[dom]] <- plot_kariocas_empty(dom, "No Data")
        next
      }

      base_dom_row <- baseline_df %>% dplyr::filter(Domain == dom, Rank == "Domain")
      base_dom_total <- if(nrow(base_dom_row) > 0) base_dom_row$Base_Global[1] else 0
      base_counts <- baseline_df %>% dplyr::filter(Domain == dom)
      get_cnt <- function(r) { val <- base_counts$Base_Taxa[base_counts$Rank == r]; if(length(val) == 0) 0 else val }

      sub_str <- paste0(
        "Reads: ", fmt_num(base_dom_total), " | ",
        "P: ", get_cnt("Phylum"), " | C: ", get_cnt("Class"), " | O: ", get_cnt("Order"), " | ",
        "F: ", get_cnt("Family"), " | G: ", get_cnt("Genus"), " | S: ", get_cnt("Species")
      )

      # NOTE: Using restored structure
      shapes_vec <- get_kariocas_shapes("ranks")
      colors_vec <- get_kariocas_colors("ranks")
      labels_vec <- get_kariocas_labels()

      p <- ggplot2::ggplot(df_dom, ggplot2::aes(x = CS, y = Pct_Taxa, color = Rank, group = Rank, shape = Rank)) +
        ggplot2::geom_line(linewidth = 1) +
        ggplot2::geom_point(size = 2.5) +
        ggplot2::scale_color_manual(values = colors_vec) +
        ggplot2::scale_shape_manual(values = shapes_vec) +
        scale_y_kariocas_log10(limits = c(0.01, 105)) +
        ggplot2::scale_x_continuous(breaks = seq(0, 100, 20), limits = c(0, 100)) +
        ggplot2::labs(
          title = dom,
          subtitle = sub_str,
          x = labels_vec$y_confidence,
          y = labels_vec$y_log10_retained
        ) +
        theme_kariocas() +
        ggplot2::guides(color = ggplot2::guide_legend(nrow = 1))

      plot_list_all[[dom]] <- p
    }

    layout_all <- (plot_list_all[["Bacteria"]] | plot_list_all[["Archaea"]]) /
      (plot_list_all[["Eukaryota"]] | plot_list_all[["Viruses"]]) +
      patchwork::plot_annotation(
        title = paste(samp, "- Retention (All Levels)"),
        subtitle = "Comparison of taxa loss across ranks",
        theme = ggplot2::theme(plot.title = ggplot2::element_text(face="bold", size=16, hjust=0.5))
      ) +
      patchwork::plot_layout(guides = "collect") &
      ggplot2::theme(legend.position = "bottom")

    fname_all <- paste0(samp, "_CS_Retention_All_Levels.pdf")
    ggplot2::ggsave(file.path(output_dir, fname_all), layout_all, width = kariocas_dims$width, height = kariocas_dims$height)
    log_msg("    -> Generated: ", fname_all)

    # ==========================================================================
    # 5. PLOT TYPE B
    # ==========================================================================
    rank_names_map <- c(
      "Phylum" = "Phyla", "Class" = "Classes", "Order" = "Orders",
      "Family" = "Families", "Genus" = "Genera", "Species" = "Species"
    )

    for (r in names(rank_names_map)) {
      fname_suffix <- rank_names_map[[r]]
      plot_list_rank <- list()

      for (dom in DOMAINS) {
        df_viz <- df_calc %>% dplyr::filter(Domain == dom, Rank == r)

        if (nrow(df_viz) == 0) {
          plot_list_rank[[dom]] <- plot_kariocas_empty(dom, "No Data")
          next
        }

        leg_taxa  <- r
        leg_reads <- paste0(r, "-exclusive Reads")
        leg_total <- "Total Reads"

        df_long_plot <- df_viz %>%
          dplyr::select(CS, Pct_Taxa, Pct_Rank_Reads, Pct_Domain_Retained) %>%
          tidyr::pivot_longer(
            cols = c("Pct_Taxa", "Pct_Rank_Reads", "Pct_Domain_Retained"),
            names_to = "Metric_Type",
            values_to = "Pct_Value"
          ) %>%
          dplyr::mutate(
            Metric_Label = dplyr::case_when(
              Metric_Type == "Pct_Taxa" ~ leg_taxa,
              Metric_Type == "Pct_Rank_Reads" ~ leg_reads,
              Metric_Type == "Pct_Domain_Retained" ~ leg_total
            ),
            Metric_Label = factor(Metric_Label, levels = c(leg_taxa, leg_total, leg_reads))
          )

        base_g <- df_viz$Base_Global[1]
        base_r <- df_viz$Base_Reads[1]
        base_t <- df_viz$Base_Taxa[1]

        sub_str_rank <- paste0(
          fmt_num(base_g), " Total Reads; ",
          fmt_num(base_r), " assigned to ", r, " | ",
          fmt_num(base_t), " ", r
        )

        spec_colors <- get_kariocas_colors("special")
        col_taxa  <- spec_colors[["Level Taxa"]]
        col_total <- spec_colors[["Total Reads"]]
        col_reads <- spec_colors[["Level Reads"]]
        
        shp_vec   <- get_kariocas_shapes("ranks")
        shp_taxa  <- shp_vec[["Level Taxa"]]
        shp_total <- shp_vec[["Total Reads"]]
        shp_reads <- shp_vec[["Level Reads"]]
        
        lt_vec     <- get_kariocas_linetypes()
        labels_vec <- get_kariocas_labels()
        lt_taxa  <- lt_vec[["Level Taxa"]]
        lt_total <- lt_vec[["Total Reads"]]
        lt_reads <- lt_vec[["Level Reads"]]

        p <- ggplot2::ggplot(df_long_plot, ggplot2::aes(x = CS, y = Pct_Value, color = Metric_Label, linetype = Metric_Label, shape = Metric_Label)) +
          ggplot2::geom_line(linewidth = 1) +
          ggplot2::geom_point(size = 3) +
          ggplot2::scale_color_manual(values = setNames(c(col_taxa, col_total, col_reads), c(leg_taxa, leg_total, leg_reads))) +
          ggplot2::scale_linetype_manual(values = setNames(c(lt_taxa, lt_total, lt_reads), c(leg_taxa, leg_total, leg_reads))) +
          ggplot2::scale_shape_manual(values = setNames(c(shp_taxa, shp_total, shp_reads), c(leg_taxa, leg_total, leg_reads))) +
          scale_y_kariocas_log10(limits = c(0.01, 105)) +
          ggplot2::scale_x_continuous(breaks = seq(0, 100, 20), limits = c(0, 100)) +
          ggplot2::labs(
            title = dom,
            subtitle = sub_str_rank,
            x = labels_vec$y_confidence,
            y = labels_vec$y_log10_retained,
            color = NULL, linetype = NULL, shape = NULL
          ) +
          theme_kariocas() +
          ggplot2::theme(legend.position = "bottom")

        plot_list_rank[[dom]] <- p
      }

      layout_rank <- (plot_list_rank[["Bacteria"]] | plot_list_rank[["Archaea"]]) /
        (plot_list_rank[["Eukaryota"]] | plot_list_rank[["Viruses"]]) +
        patchwork::plot_annotation(
          title = paste(samp, "- Retention:", fname_suffix),
          subtitle = "Comparison of Taxa vs Reads Retention",
          theme = ggplot2::theme(plot.title = ggplot2::element_text(face="bold", size=16, hjust=0.5))
        ) +
        patchwork::plot_layout(guides = "collect") &
        ggplot2::theme(legend.position = "bottom")

      fname_rank <- paste0(samp, "_CS_Retention_", fname_suffix, ".pdf")
      ggplot2::ggsave(file.path(output_dir, fname_rank), layout_rank, width = kariocas_dims$width, height = kariocas_dims$height)
      log_msg("    -> Generated: ", fname_rank)
    }
  }

  log_msg("SUCCESS: Retention analysis completed.")
  return(invisible(NULL))
}
