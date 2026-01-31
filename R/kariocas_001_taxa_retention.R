#' Run Confidence Score Retention Analysis (Step 001)
#'
#' Executes read retention analysis based on Confidence Score (Kraken/Bracken).
#' Calculates the percentage of reads and taxa retained at each CS level relative to the
#' baseline (CS00). Generates publication-quality panels using karioCaS styles.
#'
#' @param project_dir Root path of the project.
#'
#' @export
#' @importFrom readr read_rds
#' @importFrom dplyr filter select group_by summarise mutate left_join n_distinct arrange pull case_when
#' @importFrom tidyr replace_na
#' @importFrom ggplot2 ggplot aes geom_line geom_point scale_color_manual scale_shape_manual scale_linetype_manual labs ggsave scale_x_continuous coord_cartesian guides guide_legend
#' @importFrom patchwork plot_layout plot_annotation
#' @importFrom scales percent

taxa_retention <- function(project_dir) {

  # ==============================================================================
  # 1. SETUP & DIRECTORIES
  # ==============================================================================
  # Define standard paths based on karioCaS structure
  input_dir  <- file.path(project_dir, "000_karioCaS_input_matrix")
  output_dir <- file.path(project_dir, "001_cs_retention")

  # Validation
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

  # ==============================================================================
  # 2. DATA LOADING
  # ==============================================================================
  message(">>> Loading Tidy Data from: ", rds_file)
  df_long <- readr::read_rds(rds_file)

  # Prepare Taxon Name Column for Counting Unique Taxa
  # (Since 'Taxonomy' is the full string, we create a cleaner Taxon identifier)
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

  # Analysis Parameters
  SAMPLES <- unique(df_proc$sample)
  DOMAINS <- names(kariocas_colors$domains) # c("Bacteria", "Archaea", "Viruses", "Eukaryota")

  # Map Rank to Plural for filenames/titles
  rank_plural_map <- c("Phylum"="Phyla", "Class"="Classes", "Order"="Orders",
                       "Family"="Families", "Genus"="Genera", "Species"="Species")

  # Loop Targets: NULL means "All Levels Mode", others are "Specific Rank Mode"
  TARGET_RANKS <- c(NA, "Phylum", "Class", "Order", "Family", "Genus", "Species")

  message(">>> Starting Retention Analysis for ", length(SAMPLES), " samples.")

  # ==============================================================================
  # 3. CORE LOGIC LOOP (Sample -> Baseline -> Plotting)
  # ==============================================================================
  for (samp in SAMPLES) {
    message("  Processing Sample: ", samp)

    # Create sample subfolder
    samp_out_dir <- file.path(output_dir, samp)
    if (!dir.exists(samp_out_dir)) dir.create(samp_out_dir)

    df_samp <- df_proc %>% dplyr::filter(sample == samp)

    # 3.1 CALCULATE BASELINES (CS == 0) ----------------------------------------
    # Baseline A: Total Reads per Domain (CS00)
    df_baseline_domain <- df_samp %>%
      dplyr::filter(CS == 0) %>%
      dplyr::group_by(Domain) %>%
      dplyr::summarise(Base_Domain_Reads = sum(Counts), .groups = "drop")

    # Baseline B: Reads and Taxa Count per Rank per Domain (CS00)
    df_baseline_ranks <- df_samp %>%
      dplyr::filter(CS == 0) %>%
      dplyr::group_by(Domain, Rank) %>%
      dplyr::summarise(
        Base_Rank_Reads = sum(Counts),
        Base_Rank_Taxa  = dplyr::n_distinct(Taxon_Name[Counts > 0]),
        .groups = "drop"
      )

    # 3.2 CALCULATE METRICS FOR ALL CS -----------------------------------------
    # Group by Domain/CS/Rank to get current values
    df_stats <- df_samp %>%
      dplyr::group_by(Domain, CS, Rank) %>%
      dplyr::summarise(
        Curr_Reads = sum(Counts),
        Curr_Taxa  = dplyr::n_distinct(Taxon_Name[Counts > 0]),
        .groups = "drop"
      ) %>%
      # Join Baselines
      dplyr::left_join(df_baseline_domain, by = "Domain") %>%
      dplyr::left_join(df_baseline_ranks, by = c("Domain", "Rank")) %>%
      # Calculate Retention Percentages
      dplyr::mutate(
        Ret_Reads_Pct = (Curr_Reads / Base_Rank_Reads) * 100,
        Ret_Taxa_Pct  = (Curr_Taxa / Base_Rank_Taxa) * 100,
        # Domain Global Trend (Total Reads of Domain retained)
        Ret_Global_Reads_Pct = (Curr_Reads / Base_Domain_Reads) * 100
      ) %>%
      # Replace Zeros/NaN with NA to break lines in plot (Business Rule)
      dplyr::mutate(
        Ret_Reads_Pct = ifelse(is.nan(Ret_Reads_Pct) | Ret_Reads_Pct == 0, NA, Ret_Reads_Pct),
        Ret_Taxa_Pct  = ifelse(is.nan(Ret_Taxa_Pct) | Ret_Taxa_Pct == 0, NA, Ret_Taxa_Pct)
      )

    # Separate table for Domain Totals (for Specific Mode background line)
    df_domain_trend <- df_stats %>%
      dplyr::group_by(Domain, CS) %>%
      dplyr::summarise(Ret_Total_Reads_Trend = sum(Curr_Reads) / unique(Base_Domain_Reads) * 100, .groups="drop") %>%
      dplyr::mutate(Ret_Total_Reads_Trend = ifelse(Ret_Total_Reads_Trend == 0, NA, Ret_Total_Reads_Trend))

    # ==========================================================================
    # 4. VISUALIZATION LOOP (Modes: All Levels vs Specific Rank)
    # ==========================================================================
    for (lvl in TARGET_RANKS) {

      plot_list <- list()
      is_all_levels <- is.na(lvl)
      mode_tag      <- if(is_all_levels) "All_Levels" else rank_plural_map[[lvl]]

      for (dom in DOMAINS) {

        # 4.1 SUBTITLE STATS GENERATION ----------------------------------------
        # Extract Baseline stats for the subtitle
        base_dom_reads <- df_baseline_domain$Base_Domain_Reads[df_baseline_domain$Domain == dom]
        if(length(base_dom_reads) == 0) base_dom_reads <- 0
        reads_str <- label_kariocas_auto(base_dom_reads) # Using Style Formatter

        # Helper to get taxa count for subtitle
        get_cnt <- function(r) {
          val <- df_baseline_ranks$Base_Rank_Taxa[df_baseline_ranks$Domain == dom & df_baseline_ranks$Rank == r]
          if(length(val) == 0) return(0) else return(val)
        }

        # Check if we have data to plot
        has_data <- any(df_stats$Domain == dom & !is.na(df_stats$Ret_Taxa_Pct))

        if (!has_data || base_dom_reads == 0) {
          # --> NO DATA SCENARIO
          plot_list[[dom]] <- plot_kariocas_empty(
            title_text = dom,
            subtitle_text = "No reads detected at CS00",
            x_label = kariocas_labels$y_confidence,
            y_label = kariocas_labels$y_log10_retained
          )
          next
        }

        # 4.2 PLOT CONSTRUCTION ------------------------------------------------
        if (is_all_levels) {
          # --- MODE: ALL LEVELS ---
          # Subtitle: Reads: X | P: Y | C: Z ...
          subtitle_stats <- paste0(
            "Reads: ", reads_str, " | P: ", get_cnt("Phylum"), " | C: ", get_cnt("Class"),
            " | O: ", get_cnt("Order"), " | F: ", get_cnt("Family"),
            " | G: ", get_cnt("Genus"), " | S: ", get_cnt("Species")
          )

          # Filter for standard ranks only
          dat_plot <- df_stats %>%
            dplyr::filter(Domain == dom, Rank %in% names(kariocas_colors$ranks))
          # Enforce factor order for Legend
          dat_plot$Rank <- factor(dat_plot$Rank, levels = names(kariocas_colors$ranks))

          p <- ggplot2::ggplot(dat_plot, ggplot2::aes(x = CS, y = Ret_Taxa_Pct, color = Rank, shape = Rank)) +
            ggplot2::geom_line(linewidth = 0.8) +
            ggplot2::geom_point(size = 2.5) +
            ggplot2::scale_color_manual(values = kariocas_colors$ranks) +
            ggplot2::scale_shape_manual(values = kariocas_shapes)

        } else {
          # --- MODE: SPECIFIC RANK ---
          # Subtitle: X Total Reads; Y Reads assigned to Rank | Z Taxa
          curr_reads_base <- df_baseline_ranks$Base_Rank_Reads[df_baseline_ranks$Domain == dom & df_baseline_ranks$Rank == lvl]
          curr_reads_str  <- label_kariocas_auto(curr_reads_base)

          subtitle_stats <- paste0(
            reads_str, " Total Reads; ", curr_reads_str, " assigned to ", lvl,
            " | ", get_cnt(lvl), " ", mode_tag
          )

          # Prepare Data: Join specific rank retention + Global domain trend
          dat_lvl <- df_stats %>%
            dplyr::filter(Domain == dom, Rank == lvl) %>%
            dplyr::left_join(df_domain_trend, by = c("Domain", "CS")) %>%
            dplyr::select(CS, Ret_Total_Reads_Trend, Ret_Lvl_Reads = Ret_Reads_Pct, Ret_Lvl_Taxa = Ret_Taxa_Pct) %>%
            tidyr::pivot_longer(cols = -CS, names_to = "Metric", values_to = "Pct")

          # Map Metrics to Readable Names for Legend
          metric_levels <- c("Ret_Lvl_Taxa", "Ret_Total_Reads_Trend", "Ret_Lvl_Reads")
          metric_labels <- c(lvl, "Total Reads", paste(lvl, "Reads"))

          dat_lvl$Metric <- factor(dat_lvl$Metric, levels = metric_levels, labels = metric_labels)

          p <- ggplot2::ggplot(dat_lvl, ggplot2::aes(x = CS, y = Pct, group = Metric)) +
            ggplot2::geom_line(ggplot2::aes(color = Metric, linetype = Metric), linewidth = 1) +
            ggplot2::geom_point(ggplot2::aes(color = Metric, shape = Metric), size = 3) +
            ggplot2::scale_color_manual(values = kariocas_colors$special) +
            ggplot2::scale_linetype_manual(values = kariocas_linetypes) +
            # Map shapes: Total Reads = Square (Phylum shape/Base), Level Reads = Triangle (Order), Taxa = Circle (Class)
            # (Note: This maps internal style logic to specific plot needs)
            ggplot2::scale_shape_manual(values = c(16, 15, 17))
        }

        # 4.3 APPLY COMMON STYLES ----------------------------------------------
        p <- p +
          scale_y_kariocas_log10(limits = c(0.01, 105)) + # Custom Log Scale from Style
          ggplot2::scale_x_continuous(breaks = seq(0, 100, 20), limits = c(0, 100)) +
          ggplot2::labs(
            title = dom,
            subtitle = subtitle_stats,
            x = kariocas_labels$y_confidence, # "Confidence Score (%)"
            y = kariocas_labels$y_log10_retained # Math Expression
          ) +
          theme_kariocas() +
          ggplot2::coord_cartesian(clip = "off")

        plot_list[[dom]] <- p
      } # End Domain Loop

      # 4.4 ASSEMBLE PANEL -----------------------------------------------------
      # Layout: (Bacteria | Archaea) / (Eukaryota | Viruses)
      final_layout <- (plot_list[["Bacteria"]] | plot_list[["Archaea"]]) /
        (plot_list[["Eukaryota"]] | plot_list[["Viruses"]]) +
        patchwork::plot_annotation(
          title = paste(samp, "- CS Retention Analysis"),
          subtitle = paste("Retention relative to Baseline (CS00) | Mode:", mode_tag),
          theme = ggplot2::theme(plot.title = ggplot2::element_text(face="bold", size=16, hjust=0.5))
        ) +
        patchwork::plot_layout(guides = "collect") &
        ggplot2::theme(legend.position = "bottom")

      # 4.5 SAVE OUTPUT --------------------------------------------------------
      file_name <- paste0(samp, "_CS_Retention_", mode_tag, ".pdf")
      save_path <- file.path(samp_out_dir, file_name)

      ggplot2::ggsave(
        filename = save_path,
        plot = final_layout,
        width = kariocas_dims$width,
        height = kariocas_dims$height,
        units = "in"
      )

    } # End Rank Loop
  } # End Sample Loop

  message("\nSUCCESS: Retention analysis completed. Files saved in: ", output_dir)
  return(invisible(NULL))
}
