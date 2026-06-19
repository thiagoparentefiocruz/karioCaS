# ==============================================================================
# PRIVATE HELPERS — heatmaps_karioCaS()
# Not exported. No man/ page needed.
# ==============================================================================

#' @noRd
.hm_setup <- function(project_dir, analysis_rank, confidence_score) {
    output_dir <- file.path(project_dir, "005_heatmaps")
    log_dir <- file.path(project_dir, "logs")
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)
    log_file <- file.path(log_dir, "log_005_heatmaps.txt")
    header <- c(
        "====================================================",
        "LOG: 005_HEATMAP_GENERATION",
        paste0("PROJECT DIR: ", project_dir),
        paste0(
            "RANK: ", analysis_rank, " | TARGET CS: ",
            ifelse(is.null(confidence_score), "MAX", confidence_score)
        ),
        "===================================================="
    )
    writeLines(header, con = log_file)
    log_msg <- function(...) {
        msg <- paste0(...)
        message(msg)
        write(
            paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", msg),
            file = log_file, append = TRUE
        )
    }
    list(output_dir = output_dir, log_file = log_file, log_msg = log_msg)
}

#' @noRd
.hm_load_and_enrich <- function(project_dir, analysis_rank) {
    df_long <- .get_tidy_data(project_dir)
    if (!analysis_rank %in% colnames(df_long)) {
        tse_path <- file.path(
            project_dir, "000_karioCaS_input_matrix", "karioCaS_TSE.rds"
        )
        if (!file.exists(tse_path)) stop("TSE file missing for enrichment.")
        tse <- readRDS(tse_path)
        tax_df <- SummarizedExperiment::rowData(tse) |> as.data.frame()
        key_col <- "Taxonomy_Full"
        if (!key_col %in% colnames(df_long)) {
            if ("Taxonomy" %in% colnames(df_long)) {
                df_long <- dplyr::rename(df_long, Taxonomy_Full = Taxonomy)
            } else {
                stop("Could not find Taxonomy key column.")
            }
        }
        cols_to_add <- setdiff(colnames(tax_df), colnames(df_long))
        tax_sub <- dplyr::select(tax_df, "Taxonomy_Full", dplyr::all_of(cols_to_add))
        df_long <- dplyr::left_join(df_long, tax_sub, by = "Taxonomy_Full")
    }
    if (!analysis_rank %in% colnames(df_long)) {
        stop("Rank '", analysis_rank, "' not found after enrichment.")
    }
    df_long |>
        dplyr::filter(.data$Rank == analysis_rank) |>
        dplyr::mutate(Taxon_Name = .data[[analysis_rank]]) |>
        dplyr::filter(!is.na(.data$Taxon_Name))
}

#' @noRd
.hm_aggregate_loss_groups <- function(df_rest, all_cs, target_cs, analysis_rank) {
    agg_list <- list()
    for (i in seq_along(all_cs)) {
        curr_cs_val <- all_cs[i]
        curr_taxa <- unique(df_rest$Taxon_Name[df_rest$CS == curr_cs_val])
        if (curr_cs_val < target_cs) {
            next_cs_val <- all_cs[i + 1]
            if (!is.na(next_cs_val)) {
                next_taxa <- unique(df_rest$Taxon_Name[df_rest$CS == next_cs_val])
                lost_taxa <- setdiff(curr_taxa, next_taxa)
                label_suffix <- paste0("Recovered only in CS", sprintf("%02d", curr_cs_val))
            } else {
                lost_taxa <- curr_taxa
                label_suffix <- paste0(
                    "Lowest abundance ", analysis_rank,
                    " in CS", sprintf("%02d", curr_cs_val)
                )
            }
        } else {
            lost_taxa <- curr_taxa
            label_suffix <- paste0(
                "Lowest abundance ", analysis_rank,
                " in CS", sprintf("%02d", curr_cs_val)
            )
        }
        if (length(lost_taxa) > 0) {
            row_label <- paste(length(lost_taxa), analysis_rank, label_suffix)
            agg_data <- df_rest |>
                dplyr::filter(.data$Taxon_Name %in% lost_taxa) |>
                dplyr::group_by(.data$CS) |>
                dplyr::summarise(Counts = sum(.data$Counts), .groups = "drop") |>
                dplyr::mutate(Taxon_Name = row_label)
            agg_list[[length(agg_list) + 1]] <- agg_data
        }
    }
    dplyr::bind_rows(agg_list)
}

#' @noRd
.hm_cluster_elite <- function(df_metrics, elite_taxa_names) {
    if (length(elite_taxa_names) <= 2) {
        return(elite_taxa_names)
    }
    mat_prep <- df_metrics |>
        dplyr::filter(.data$Taxon_Name %in% elite_taxa_names) |>
        dplyr::select("Taxon_Name", "CS", "Rel_Abund") |>
        dplyr::distinct(.data$Taxon_Name, .data$CS, .keep_all = TRUE)
    mat_elite <- mat_prep |>
        tidyr::pivot_wider(
            names_from  = "CS",
            values_from = "Rel_Abund",
            values_fill = 0
        ) |>
        tibble::column_to_rownames("Taxon_Name") |>
        as.matrix()
    mat_elite[is.na(mat_elite)] <- 0
    row_means <- rowMeans(mat_elite)
    hc <- stats::hclust(stats::dist(mat_elite, method = "euclidean"),
        method = "complete"
    )
    dend <- stats::as.dendrogram(hc)
    dend_reorder <- stats::reorder(dend, row_means, agglo.FUN = mean)
    labels(dend_reorder)
}

#' @noRd
.hm_process_domain <- function(df_samp, dom, target_cs, top_n,
                               analysis_rank, all_cs, log_msg) {
    df_dom <- df_samp |>
        dplyr::filter(.data$Domain == dom) |>
        dplyr::group_by(.data$CS, .data$Taxon_Name) |>
        dplyr::summarise(Counts = sum(.data$Counts, na.rm = TRUE), .groups = "drop")
    if (nrow(df_dom) == 0) {
        return(NULL)
    }
    survivors <- df_dom |>
        dplyr::filter(.data$CS == target_cs) |>
        dplyr::arrange(dplyr::desc(.data$Counts))
    elite_names <- dplyr::slice_head(survivors, n = top_n) |>
        dplyr::pull(.data$Taxon_Name)
    if (length(elite_names) == 0) {
        log_msg("    Info: ", dom, " — no survivors at target CS. Plotting loss groups only.")
    }
    df_elite <- dplyr::filter(df_dom, .data$Taxon_Name %in% elite_names)
    df_rest <- dplyr::filter(df_dom, !.data$Taxon_Name %in% elite_names)
    df_agg <- .hm_aggregate_loss_groups(df_rest, all_cs, target_cs, analysis_rank)
    df_combined <- dplyr::bind_rows(df_elite, df_agg)
    if (nrow(df_combined) == 0) {
        return(NULL)
    }
    total_per_cs <- df_dom |>
        dplyr::group_by(.data$CS) |>
        dplyr::summarise(Total_CS_Reads = sum(.data$Counts), .groups = "drop")
    df_metrics <- df_combined |>
        dplyr::left_join(total_per_cs, by = "CS") |>
        dplyr::mutate(Rel_Abund = (.data$Counts / .data$Total_CS_Reads) * 100) |>
        tidyr::complete(Taxon_Name,
            CS = all_cs,
            fill = list(Rel_Abund = 0, Counts = 0)
        ) |>
        dplyr::mutate(Domain = dom, CS_Label = sprintf("%02d", .data$CS))
    elite_order <- .hm_cluster_elite(df_metrics, elite_names)
    agg_names_ord <- if (nrow(df_agg) > 0) unique(df_agg$Taxon_Name) else character(0)
    final_levels <- c(agg_names_ord, elite_order)
    dplyr::mutate(df_metrics, Order_Index = match(.data$Taxon_Name, final_levels))
}

#' @noRd
.hm_build_heatmap <- function(full_plot_data, analysis_rank, samp,
                              target_cs, barro_palette) {
    ggplot2::ggplot(
        full_plot_data,
        ggplot2::aes(x = .data$CS_Label, y = .data$Taxon_Name, fill = .data$Rel_Abund)
    ) +
        ggplot2::geom_tile(color = "white", linewidth = 0.2) +
        ggplot2::facet_grid(Domain ~ ., scales = "free_y", space = "free_y") +
        ggplot2::scale_fill_gradientn(
            colors = barro_palette, name = "Rel. Abund. (%)",
            values = c(0, 0.1, 0.25, 0.6, 1), limits = c(0, 100)
        ) +
        ggplot2::guides(fill = ggplot2::guide_colorbar(
            title.position = "right", title.hjust = 1, title.vjust = 1,
            barwidth = ggplot2::unit(7, "cm"),
            barheight = ggplot2::unit(0.3, "cm")
        )) +
        ggplot2::labs(
            title = paste(analysis_rank, "Relative Abundance -", samp),
            subtitle = paste(
                "Detailed aggregation (Threshold: CS",
                sprintf("%02d", target_cs), ")"
            ),
            x = "Kraken Confidence Score (%)", y = NULL
        ) +
        theme_kariocas() +
        ggplot2::theme(
            axis.text.x        = ggplot2::element_text(angle = 0, hjust = 0.5, size = 10),
            axis.text.y        = ggplot2::element_text(size = 9, face = "italic"),
            axis.ticks.y       = ggplot2::element_blank(),
            strip.text.y       = ggplot2::element_text(angle = 0, face = "bold", size = 10),
            strip.background   = ggplot2::element_rect(fill = "grey95", color = NA),
            panel.grid         = ggplot2::element_blank(),
            axis.line          = ggplot2::element_blank(),
            panel.spacing      = ggplot2::unit(0.5, "cm"),
            legend.position    = "bottom",
            legend.direction   = "horizontal",
            legend.title       = ggplot2::element_text(size = 8, face = "plain"),
            legend.box.spacing = ggplot2::unit(0.8, "cm")
        )
}

#' @noRd
.hm_process_sample <- function(df_proc, samp, DOMAINS, confidence_score, top_n,
                               analysis_rank, output_dir, log_msg, barro_palette) {
    log_msg("------------------------------------------------")
    log_msg("  Processing Sample: ", samp)
    df_samp <- dplyr::filter(df_proc, .data$sample == samp)
    max_cs <- .safe_max(df_samp$CS, default = 0)
    target_cs <- if (!is.null(confidence_score) && confidence_score <= max_cs) {
        confidence_score
    } else {
        max_cs
    }
    df_samp <- dplyr::filter(df_samp, .data$CS <= target_cs)
    all_cs <- sort(unique(df_samp$CS))
    domain_results <- lapply(DOMAINS, function(dom) {
        .hm_process_domain(
            df_samp, dom, target_cs, top_n,
            analysis_rank, all_cs, log_msg
        )
    })
    full_plot_data <- dplyr::bind_rows(Filter(Negate(is.null), domain_results))
    if (nrow(full_plot_data) == 0) {
        log_msg("    No data to plot for this sample.")
        return(invisible(NULL))
    }
    full_plot_data <- full_plot_data |>
        dplyr::arrange(.data$Domain, .data$Order_Index) |>
        dplyr::mutate(
            Taxon_Name = factor(.data$Taxon_Name, levels = unique(.data$Taxon_Name))
        )
    p <- .hm_build_heatmap(
        full_plot_data, analysis_rank, samp,
        target_cs, barro_palette
    )
    file_name <- paste0(
        samp, "_Heatmap_", analysis_rank,
        "_CS", sprintf("%02d", target_cs), ".pdf"
    )
    ggplot2::ggsave(
        file.path(output_dir, file_name), p,
        width = get_kariocas_dims()$height,
        height = get_kariocas_dims()$width
    )
    log_msg("    -> Generated: ", file_name)
}

# ==============================================================================
# EXPORTED FUNCTION
# ==============================================================================

#' Generate Heatmaps of Taxa Abundance with Extinction Patterns (Step 005)
#'
#' Creates Relative Abundance (%) heatmaps focused on survivors at a specific
#' Confidence Score. "Elite" survivors are clustered by similarity, while lost
#' taxa are aggregated into "Loss Groups" at the bottom. Handles domains with
#' zero survivors at the target CS (showing only loss groups).
#'
#' @param project_dir Path to the project root.
#' @param analysis_rank Taxonomic rank to analyze. Defaults to \code{"Genus"}.
#' @param confidence_score Target CS to define survivors (e.g., 90).
#'   Defaults to the highest available CS.
#' @param top_n Number of top survivors to display individually (default: 20).
#'
#' @return Invisibly returns \code{NULL}. PDF plots are saved to
#'   \code{<project_dir>/005_heatmaps/}.
#' @export
#' @importFrom dplyr filter mutate select group_by summarise arrange slice_head
#'   pull ungroup left_join bind_rows distinct rename all_of desc
#' @importFrom tidyr complete pivot_longer pivot_wider replace_na
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_gradientn labs theme
#'   element_text element_blank scale_x_discrete scale_y_discrete guides
#'   guide_colorbar facet_grid unit element_rect ggsave
#' @importFrom stats hclust dist as.dendrogram reorder
#' @importFrom forcats fct_reorder fct_inorder
#' @importFrom tibble column_to_rownames
#' @importFrom SummarizedExperiment rowData
#' @examples
#' toy_project <- system.file("extdata", "your_project_name", package = "karioCaS")
#'
#' # Basic usage (defaults to Genus, highest CS, top 20 taxa)
#' # heatmaps_karioCaS(project_dir = toy_project)
#'
#' # Advanced: Species level at CS 50, top 30 taxa
#' # heatmaps_karioCaS(
#' #   project_dir      = toy_project,
#' #   analysis_rank    = "Species",
#' #   confidence_score = 50,
#' #   top_n            = 30
#' # )
heatmaps_karioCaS <- function(project_dir,
                              analysis_rank = NULL,
                              confidence_score = NULL,
                              top_n = 20) {
    if (is.null(analysis_rank)) analysis_rank <- "Genus"
    if (!is.null(confidence_score)) {
        cs_pct <- .cs_arg_to_percent(confidence_score)
        if (is.na(cs_pct)) {
            stop(
                "Invalid 'confidence_score': ", confidence_score,
                ". Use a Kraken fraction (0-1) or a percentage (0-100)."
            )
        }
        confidence_score <- cs_pct
    }
    setup <- .hm_setup(project_dir, analysis_rank, confidence_score)
    df_proc <- .hm_load_and_enrich(project_dir, analysis_rank)
    barro_pal <- c("#FFFFFF", "#FFEDA0", "#FEB24C", "#F03B20", "#800026")
    SAMPLES <- unique(df_proc$sample)
    DOMAINS <- names(get_kariocas_colors("domains"))
    setup$log_msg(">>> Starting Heatmap Analysis for ", length(SAMPLES), " samples.")
    for (samp in SAMPLES) {
        .hm_process_sample(
            df_proc, samp, DOMAINS, confidence_score, top_n,
            analysis_rank, setup$output_dir, setup$log_msg, barro_pal
        )
    }
    setup$log_msg("SUCCESS: Heatmap analysis completed.")
    invisible(NULL)
}
