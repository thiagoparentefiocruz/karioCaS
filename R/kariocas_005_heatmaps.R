#' Generate Heatmaps of Taxa Abundance with Extinction Patterns (Step 005)
#'
#' Creates Relative Abundance (%) heatmaps focused on survivors at a specific
#' Confidence Score. "Elite" survivors are clustered by similarity, while lost taxa
#' are aggregated into "Loss Groups" at the bottom.
#' Includes robust data enrichment and count consolidation to prevent abundance dilution.
#' Now explicitly handles domains with zero survivors at target CS (showing only loss groups).
#'
#' @param project_dir Path to the project root.
#' @param analysis_rank Taxonomic rank to analyze. If NULL, defaults to "Genus".
#' @param confidence_score Target CS to define "Survivors" (e.g., 90). If NULL, uses the highest available.
#' @param top_n Number of top survivors to display individually (default: 20).
#'
#' @return Generates PDF plots saved in the project directory and returns a data frame invisibly.
#' @export
#' @importFrom dplyr filter mutate select group_by summarise arrange slice_head pull case_when ungroup left_join bind_rows distinct n_distinct rename
#' @importFrom tidyr complete pivot_longer pivot_wider replace_na
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_gradientn labs theme element_text element_blank scale_x_discrete scale_y_discrete guides guide_colorbar facet_grid unit element_rect
#' @importFrom stats hclust dist as.dendrogram reorder
#' @importFrom forcats fct_reorder fct_inorder
#' @importFrom tibble column_to_rownames
#' @importFrom SummarizedExperiment rowData
#' @examples
#' # Get the path to the included toy dataset
#' toy_project <- system.file("extdata", "your_project_name", package = "karioCaS")
#'
#' # Example 1: Basic usage (defaults to Genus, highest CS, and top 20 taxa)
#' # heatmaps_karioCaS(
#' #   project_dir = toy_project
#' # )
#' 
#' # Example 2: Advanced usage (Targeting Species level at CS 50, showing top 30 taxa)
#' # heatmaps_karioCaS(
#' #   project_dir = toy_project,
#' #   analysis_rank = "Species",
#' #   confidence_score = 50,
#' #   top_n = 30
#' # )

heatmaps_karioCaS <- function(project_dir,
                              analysis_rank = NULL,
                              confidence_score = NULL,
                              top_n = 20) {

  # ==============================================================================
  # 1. SETUP & LOGGING
  # ==============================================================================
  if (is.null(analysis_rank)) {
    analysis_rank <- "Genus"
  }

  output_dir <- file.path(project_dir, "005_heatmaps")
  log_dir    <- file.path(project_dir, "logs")

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)

  log_file <- file.path(log_dir, "log_005_heatmaps.txt")

  log_msg <- function(...) {
    msg <- paste0(...)
    message(msg)
    cat(paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", msg, "\n"),
        file = log_file, append = TRUE)
  }

  cat("====================================================\n", file = log_file)
  cat("LOG: 005_HEATMAP_GENERATION (Loss Groups Fix)\n", file = log_file, append = TRUE)
  cat("PROJECT DIR: ", project_dir, "\n", file = log_file, append = TRUE)
  cat("RANK: ", analysis_rank, " | TARGET CS: ", ifelse(is.null(confidence_score), "MAX", confidence_score), "\n", file = log_file, append = TRUE)
  cat("====================================================\n", file = log_file, append = TRUE)

  # ==============================================================================
  # 2. DATA LOADING & ENRICHMENT
  # ==============================================================================
  log_msg(">>> Loading Data...")
  df_long <- .get_tidy_data(project_dir)

  # 2.1 Enrich Data with Taxonomy Columns
  if (!analysis_rank %in% colnames(df_long)) {
    log_msg(">>> Enriching data with full taxonomy from TSE...")
    tse_path <- file.path(project_dir, "000_karioCaS_input_matrix", "karioCaS_TSE.rds")
    if (!file.exists(tse_path)) stop("TSE file missing for enrichment.")
    tse <- readRDS(tse_path)
    tax_df <- SummarizedExperiment::rowData(tse) %>% as.data.frame()

    key_col <- "Taxonomy_Full"
    if (!key_col %in% colnames(df_long)) {
      if ("Taxonomy" %in% colnames(df_long)) df_long <- df_long %>% dplyr::rename(Taxonomy_Full = Taxonomy)
      else stop("Could not find Taxonomy key column.")
    }
    cols_to_add <- setdiff(colnames(tax_df), colnames(df_long))
    tax_subset <- tax_df %>% dplyr::select(Taxonomy_Full, dplyr::all_of(cols_to_add))
    df_long <- df_long %>% dplyr::left_join(tax_subset, by = "Taxonomy_Full")
  }

  if (!analysis_rank %in% colnames(df_long)) {
    stop("Rank '", analysis_rank, "' could not be found even after enrichment.")
  }

  # 2.2 Filter and Set Taxon_Name Explicitly
  df_proc <- df_long %>%
    dplyr::filter(Rank == analysis_rank) %>%
    dplyr::mutate(Taxon_Name = .data[[analysis_rank]]) %>%
    dplyr::filter(!is.na(Taxon_Name))

  SAMPLES <- unique(df_proc$sample)
  DOMAINS <- names(kariocas_colors$domains)

  # Custom "Barro" Palette (White -> Yellow -> Orange -> Red -> Deep Terracotta)
  barro_palette <- c("#FFFFFF", "#FFEDA0", "#FEB24C", "#F03B20", "#800026")

  log_msg(">>> Starting Heatmap Analysis for ", length(SAMPLES), " samples.")

  # ==============================================================================
  # 3. ANALYSIS LOOP
  # ==============================================================================
  for (samp in SAMPLES) {
    log_msg("------------------------------------------------")
    log_msg("  Processing Sample: ", samp)

    df_samp <- df_proc %>% dplyr::filter(sample == samp)

    # 3.0 DETERMINE TARGET CS
    max_cs_in_data <- max(df_samp$CS, na.rm = TRUE)
    if (!is.null(confidence_score)) {
      target_cs <- confidence_score
      if (target_cs > max_cs_in_data) target_cs <- max_cs_in_data
    } else {
      target_cs <- max_cs_in_data
    }

    df_samp <- df_samp %>% dplyr::filter(CS <= target_cs)
    all_cs  <- sort(unique(df_samp$CS))

    domain_data_list <- list()

    for (dom in DOMAINS) {

      # 3.0.1 CONSOLIDATE DUPLICATES (Normalization Fix)
      df_dom <- df_samp %>%
        dplyr::filter(Domain == dom) %>%
        dplyr::group_by(CS, Taxon_Name) %>%
        dplyr::summarise(Counts = sum(Counts, na.rm = TRUE), .groups = "drop")

      if (nrow(df_dom) == 0) next

      # 3.1 IDENTIFY TOP N *SURVIVORS* AT TARGET CS
      survivors_at_target <- df_dom %>%
        dplyr::filter(CS == target_cs) %>%
        dplyr::arrange(dplyr::desc(Counts))

      elite_taxa_names <- survivors_at_target %>%
        dplyr::slice_head(n = top_n) %>%
        dplyr::pull(Taxon_Name)

      # FIX: DO NOT SKIP if elite_taxa_names is empty.
      # Even with 0 survivors, we must process the loss groups.
      if(length(elite_taxa_names) == 0) {
        log_msg("    Info: ", dom, " has NO survivors at Target CS. Plotting only loss groups.")
      }

      # 3.2 SPLIT DATA
      df_elite <- df_dom %>%
        dplyr::filter(Taxon_Name %in% elite_taxa_names)

      df_rest <- df_dom %>%
        dplyr::filter(!Taxon_Name %in% elite_taxa_names)

      # 3.3 AGGREGATION LOGIC (Cascading Loss)
      agg_list <- list()

      for (i in seq_along(all_cs)) {
        curr_cs_val <- all_cs[i]
        curr_taxa <- unique(df_rest$Taxon_Name[df_rest$CS == curr_cs_val])

        if (curr_cs_val < target_cs) {
          next_cs_val <- all_cs[i+1]
          if(!is.na(next_cs_val)) {
            next_taxa   <- unique(df_rest$Taxon_Name[df_rest$CS == next_cs_val])
            lost_taxa   <- setdiff(curr_taxa, next_taxa)
            label_suffix <- paste0("Recovered only in CS", sprintf("%02d", curr_cs_val))
          } else {
            lost_taxa <- curr_taxa
            label_suffix <- paste0("Lowest abundance ", analysis_rank, " in CS", sprintf("%02d", curr_cs_val))
          }
        } else {
          lost_taxa    <- curr_taxa
          label_suffix <- paste0("Lowest abundance ", analysis_rank, " in CS", sprintf("%02d", curr_cs_val))
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

      # Safety check: if everything is empty (no survivors AND no loss groups?), skip.
      if (nrow(df_combined) == 0) next

      # 3.4 CALCULATE RELATIVE ABUNDANCE (%)
      total_reads_per_cs <- df_dom %>%
        dplyr::group_by(CS) %>%
        dplyr::summarise(Total_CS_Reads = sum(Counts), .groups = "drop")

      df_metrics <- df_combined %>%
        dplyr::left_join(total_reads_per_cs, by = "CS") %>%
        dplyr::mutate(Rel_Abund = (Counts / Total_CS_Reads) * 100) %>%
        tidyr::complete(Taxon_Name, CS = all_cs, fill = list(Rel_Abund = 0, Counts = 0)) %>%
        dplyr::mutate(
          Domain = dom,
          CS_Label = sprintf("%02d", CS)
        )

      # 3.5 FACTOR ORDERING (CLUSTERING + REORDER)

      # A) Cluster Elite Taxa based on Abundance Profile
      if (length(elite_taxa_names) > 2) {
        mat_prep <- df_metrics %>%
          dplyr::filter(Taxon_Name %in% elite_taxa_names) %>%
          dplyr::select(Taxon_Name, CS, Rel_Abund) %>%
          dplyr::distinct(Taxon_Name, CS, .keep_all = TRUE)

        mat_elite <- mat_prep %>%
          tidyr::pivot_wider(names_from = CS, values_from = Rel_Abund, values_fill = 0) %>%
          tibble::column_to_rownames("Taxon_Name") %>%
          as.matrix()

        mat_elite[is.na(mat_elite)] <- 0

        row_means <- rowMeans(mat_elite)
        dist_mat <- stats::dist(mat_elite, method = "euclidean")
        hc <- stats::hclust(dist_mat, method = "complete")

        dend <- stats::as.dendrogram(hc)
        dend_reordered <- stats::reorder(dend, row_means, agglo.FUN = mean)

        elite_order <- labels(dend_reordered)
      } else {
        elite_order <- elite_taxa_names
      }

      # B) Order Aggregated Rows (Bottom Up)
      if(!is.null(df_agg) && nrow(df_agg) > 0) {
        agg_names_ordered <- unique(df_agg$Taxon_Name)
      } else {
        agg_names_ordered <- character(0)
      }

      final_levels <- c(agg_names_ordered, elite_order)

      df_metrics <- df_metrics %>%
        dplyr::mutate(Order_Index = match(Taxon_Name, final_levels))

      domain_data_list[[dom]] <- df_metrics

    } # End Domain Loop

    full_plot_data <- dplyr::bind_rows(domain_data_list)

    if (nrow(full_plot_data) == 0) {
      log_msg("    No data to plot for this sample.")
      next
    }

    full_plot_data <- full_plot_data %>%
      dplyr::arrange(Domain, Order_Index) %>%
      dplyr::mutate(Taxon_Name = factor(Taxon_Name, levels = unique(Taxon_Name)))

    # 3.6 PLOT CONSTRUCTION
    p <- ggplot2::ggplot(full_plot_data, ggplot2::aes(x = CS_Label, y = Taxon_Name, fill = Rel_Abund)) +
      ggplot2::geom_tile(color = "white", linewidth = 0.2) +

      ggplot2::facet_grid(Domain ~ ., scales = "free_y", space = "free_y") +

      ggplot2::scale_fill_gradientn(
        colors = barro_palette,
        name = "Rel. Abund. (%)",
        values = c(0, 0.1, 0.25, 0.6, 1),
        limits = c(0, 100)
      ) +

      ggplot2::guides(fill = ggplot2::guide_colorbar(
        title.position = "right",
        title.hjust = 1,
        title.vjust = 1,
        barwidth = ggplot2::unit(7, "cm"),
        barheight = ggplot2::unit(0.3, "cm")
      )) +

      ggplot2::labs(
        title = paste(analysis_rank, "Relative Abundance -", samp),
        subtitle = paste("Detailed aggregation (Threshold: CS", sprintf("%02d", target_cs), ")"),
        x = "Kraken Confidence Score (%)",
        y = NULL
      ) +

      theme_kariocas() +

      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5, size = 10),
        axis.text.y = ggplot2::element_text(size = 9, face = "italic"),
        axis.ticks.y = ggplot2::element_blank(),
        strip.text.y = ggplot2::element_text(angle = 0, face = "bold", size = 10),
        strip.background = ggplot2::element_rect(fill = "grey95", color = NA),
        panel.grid = ggplot2::element_blank(),
        axis.line = ggplot2::element_blank(),
        panel.spacing = ggplot2::unit(0.5, "cm"),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = ggplot2::element_text(size = 8, face = "plain"),
        legend.box.spacing = ggplot2::unit(0.8, "cm")
      )

    file_name <- paste0(samp, "_Heatmap_", analysis_rank, "_CS", sprintf("%02d", target_cs), ".pdf")
    save_path <- file.path(output_dir, file_name)

    ggplot2::ggsave(save_path, p, width = kariocas_dims$height, height = kariocas_dims$width)
    log_msg("    -> Generated: ", file_name)
  }

  log_msg("SUCCESS: Heatmap analysis completed.")
  return(invisible(NULL))
}
