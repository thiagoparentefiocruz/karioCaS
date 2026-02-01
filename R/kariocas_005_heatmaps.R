#' Generate Heatmaps of Taxa Abundance with Extinction Patterns (Step 005)
#'
#' Creates Relative Abundance (%) heatmaps displaying the dominance of specific
#' Top Taxa and the aggregated abundance of taxa lost at each filtering step.
#' Uses a vertical faceted layout and a White-to-Terracotta palette.
#'
#' @param project_dir Path to the project root.
#' @param analysis_rank Taxonomic rank to analyze (default: "Class").
#' @param top_n Number of specific taxa to display individually (default: 15).
#'
#' @export
#' @importFrom readr read_rds
#' @importFrom dplyr filter mutate select group_by summarise arrange slice_head pull case_when ungroup left_join bind_rows distinct n_distinct
#' @importFrom tidyr complete pivot_longer pivot_wider replace_na
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_gradientn labs theme element_text element_blank scale_x_discrete scale_y_discrete guides guide_colorbar facet_grid
#' @importFrom forcats fct_reorder fct_inorder

heatmaps_karioCaS <- function(project_dir,
                              analysis_rank = "Class",
                              top_n = 15) {

  # ==============================================================================
  # 1. SETUP & DIRECTORIES
  # ==============================================================================
  input_dir  <- file.path(project_dir, "000_karioCaS_input_matrix")
  output_dir <- file.path(project_dir, "005_heatmaps")

  if (!dir.exists(input_dir)) stop("CRITICAL ERROR: Input directory not found: ", input_dir)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  rds_file <- file.path(input_dir, "karioCaS_input_matrix.rds")
  if (!file.exists(rds_file)) stop("Input RDS file not found: ", rds_file)

  # ==============================================================================
  # 2. DATA LOADING
  # ==============================================================================
  message(">>> Loading Tidy Data from: ", rds_file)
  df_long <- readr::read_rds(rds_file)

  df_proc <- df_long %>%
    dplyr::filter(Rank == analysis_rank) %>%
    dplyr::mutate(Taxon_Name = .data[[analysis_rank]])

  SAMPLES <- unique(df_proc$sample)
  DOMAINS <- names(kariocas_colors$domains)

  # Custom "Barro" Palette (White -> Yellow -> Orange -> Red -> Deep Terracotta)
  barro_palette <- c("#FFFFFF", "#FFEDA0", "#FEB24C", "#F03B20", "#800026")

  message(">>> Starting Heatmap Analysis (Rank: ", analysis_rank, " | Top: ", top_n, ")")

  # ==============================================================================
  # 3. ANALYSIS LOOP
  # ==============================================================================
  for (samp in SAMPLES) {
    message("  Processing Sample: ", samp)

    df_samp <- df_proc %>% dplyr::filter(sample == samp)
    all_cs  <- sort(unique(df_samp$CS))

    # We will collect data frames for all domains and combine them
    # to plot in a single FacetGrid (Vertical Layout)
    domain_data_list <- list()

    for (dom in DOMAINS) {

      df_dom <- df_samp %>% dplyr::filter(Domain == dom)

      if (nrow(df_dom) == 0) next # Skip domain if empty

      # 3.1 IDENTIFY TOP N SURVIVORS (Mean Abundance Strategy)
      # We stick to simple Mean Abundance to match the "Old Script" feel.
      top_taxa_names <- df_dom %>%
        dplyr::group_by(Taxon_Name) %>%
        dplyr::summarise(Mean_Abund = mean(Counts)) %>%
        dplyr::arrange(dplyr::desc(Mean_Abund)) %>%
        dplyr::slice_head(n = top_n) %>%
        dplyr::pull(Taxon_Name)

      # 3.2 SPLIT DATA
      df_elite <- df_dom %>%
        dplyr::filter(Taxon_Name %in% top_taxa_names) %>%
        dplyr::select(CS, Taxon_Name, Counts)

      df_rest <- df_dom %>%
        dplyr::filter(!Taxon_Name %in% top_taxa_names)

      # 3.3 AGGREGATION LOGIC (Diagonal Loss Pattern)
      agg_list <- list()

      for (i in seq_along(all_cs)) {
        curr_cs_val <- all_cs[i]

        if (i < length(all_cs)) {
          next_cs_val <- all_cs[i+1]
          curr_taxa <- unique(df_rest$Taxon_Name[df_rest$CS == curr_cs_val])
          next_taxa <- unique(df_rest$Taxon_Name[df_rest$CS == next_cs_val])

          lost_taxa <- setdiff(curr_taxa, next_taxa)
          label_suffix <- paste0("Recovered only in CS", sprintf("%02d", curr_cs_val))
        } else {
          lost_taxa <- unique(df_rest$Taxon_Name[df_rest$CS == curr_cs_val])
          label_suffix <- "(Low Abundance Survivors)"
        }

        if (length(lost_taxa) > 0) {
          row_label <- paste(length(lost_taxa), analysis_rank, label_suffix)

          agg_data <- df_rest %>%
            dplyr::filter(Taxon_Name %in% lost_taxa) %>%
            dplyr::group_by(CS) %>%
            dplyr::summarise(Counts = sum(Counts), .groups = "drop") %>%
            dplyr::mutate(Taxon_Name = row_label)

          agg_list[[length(agg_list) + 1]] <- agg_data
        }
      }

      df_agg <- dplyr::bind_rows(agg_list)
      df_combined <- dplyr::bind_rows(df_elite, df_agg)

      # 3.4 CALCULATE RELATIVE ABUNDANCE (%)
      # Denominator: Total Reads of the CURRENT CS
      total_reads_per_cs <- df_dom %>%
        dplyr::group_by(CS) %>%
        dplyr::summarise(Total_CS_Reads = sum(Counts), .groups = "drop")

      df_metrics <- df_combined %>%
        dplyr::left_join(total_reads_per_cs, by = "CS") %>%
        dplyr::mutate(Rel_Abund = (Counts / Total_CS_Reads) * 100) %>%
        # Fill missing with 0
        tidyr::complete(Taxon_Name, CS = all_cs, fill = list(Rel_Abund = 0, Counts = 0)) %>%
        dplyr::mutate(
          Domain = dom,
          CS_Label = sprintf("%02d", CS)
        )

      # 3.5 FACTOR ORDERING (Per Domain)
      # Elite: Descending Abundance
      elite_order <- df_elite %>%
        dplyr::group_by(Taxon_Name) %>%
        dplyr::summarise(Tot = sum(Counts)) %>%
        dplyr::arrange(Tot) %>%
        dplyr::pull(Taxon_Name)

      # Aggregated: Diagonal (CS00 at bottom)
      agg_names_ordered <- unique(df_agg$Taxon_Name)
      agg_names_ordered <- rev(agg_names_ordered)

      # Combine levels
      final_levels <- c(agg_names_ordered, elite_order)

      # Apply factor ONLY for this chunk logic, but we need to concatenate later.
      # We will store the Order index to enforce it in the combined plot.
      df_metrics <- df_metrics %>%
        dplyr::mutate(Order_Index = match(Taxon_Name, final_levels))

      domain_data_list[[dom]] <- df_metrics

    } # End Domain Loop

    # Combine all domains
    full_plot_data <- dplyr::bind_rows(domain_data_list)

    if (nrow(full_plot_data) == 0) next

    # Enforce Ordering across the facets
    # We use reorder within the plot, or factor it now.
    # Since Taxon Names are unique per Domain (mostly), we can just factor globally
    # but we need to ensure the group structure holds.
    # Trick: Paste Domain to Name to ensure uniqueness if needed, but usually Names are unique enough.
    # Let's rely on reordering by the Order_Index we created.

    full_plot_data <- full_plot_data %>%
      dplyr::arrange(Domain, Order_Index) %>%
      dplyr::mutate(Taxon_Name = factor(Taxon_Name, levels = unique(Taxon_Name)))

    # 3.6 PLOT CONSTRUCTION (Facet Grid Style)
    p <- ggplot2::ggplot(full_plot_data, ggplot2::aes(x = CS_Label, y = Taxon_Name, fill = Rel_Abund)) +
      ggplot2::geom_tile(color = "white", linewidth = 0.2) +

      # Facet Grid: Vertical Stack
      ggplot2::facet_grid(Domain ~ ., scales = "free_y", space = "free_y") +

      # Palette: White to Terracotta (Raw Percentage)
      ggplot2::scale_fill_gradientn(
        colors = barro_palette,
        name = "Rel. Abund (%)",
        limits = c(0, 100) # Fixed 0-100 scale implies we expect dominance
        # If max is small (e.g. 10%), this might look faint.
        # If so, remove limits to autoscale, but user asked for "absolute" feel.
        # Let's leave limits OFF to allow autoscale per sample max,
        # or Set limits=c(0, NA)
      ) +

      ggplot2::labs(
        title = paste(samp, "- Taxonomic Profile & Loss Patterns"),
        subtitle = paste("Rank:", analysis_rank, "| Top", top_n, "Survivors + Aggregated Loss"),
        x = "Kraken Confidence Score (%)",
        y = NULL
      ) +

      theme_kariocas() +

      ggplot2::theme(
        # X Axis: Upright, Size 9
        axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5, size = 9),
        # Y Axis: Italic
        axis.text.y = ggplot2::element_text(size = 9, face = "italic"),
        # Strips (Facet Headers)
        strip.text = ggplot2::element_text(face = "bold", size = 10),
        strip.background = ggplot2::element_rect(fill = "grey95", color = NA),
        # Clean up
        panel.grid = ggplot2::element_blank(),
        axis.line = ggplot2::element_blank(),
        panel.spacing = ggplot2::unit(0.5, "cm"), # Space between domains
        legend.position = "right"
      )

    # 3.7 SAVE (Vertical A4)
    # Since it's a long vertical list, we might need Portrait orientation or a taller PDF
    file_name <- paste0(samp, "_Heatmap_", analysis_rank, ".pdf")
    save_path <- file.path(output_dir, file_name)

    # We use standard width, but height might need to be dynamic or A4 Portrait
    # Let's stick to kariocas_dims but maybe swapped if list is huge?
    # Keeping A4 Landscape for now as requested "4-in-1" style but vertical stack fits there too.

    ggplot2::ggsave(save_path, p, width = kariocas_dims$width, height = kariocas_dims$height)

  }

  message("\nSUCCESS: Heatmap analysis completed. Files saved in: ", output_dir)
  return(invisible(NULL))
}
