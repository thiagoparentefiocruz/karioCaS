# ==============================================================================
# PRIVATE HELPERS - taxa_resolution()
# ==============================================================================

#' @noRd
.txr_setup <- function(project_dir, parent_level, child_level) {
    output_dir <- file.path(project_dir, "004_taxa_resolution")
    log_dir <- file.path(project_dir, "logs")
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)
    log_file <- file.path(log_dir, "log_004_taxa_resolution.txt")
    writeLines(c(
        "====================================================",
        "LOG: 004_TAXA_RESOLUTION (Max-Based + Doc Fix)",
        paste0("PROJECT DIR: ", project_dir),
        paste0("ANALYSIS: ", parent_level, " vs ", child_level),
        "===================================================="
    ), con = log_file)
    log_msg <- function(...) {
        msg <- paste0(...)
        message(msg)
        write(
            paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", msg),
            file = log_file, append = TRUE
        )
    }
    list(output_dir = output_dir, log_msg = log_msg)
}

#' @noRd
.txr_load_and_enrich <- function(project_dir, parent_level, child_level, log_msg) {
    log_msg(">>> Loading Data...")
    df_long <- .get_tidy_data(project_dir)
    required_cols <- c(parent_level, child_level)
    if (!all(required_cols %in% colnames(df_long))) {
        log_msg(">>> Enriching data with full taxonomy from TSE...")
        tse_path <- file.path(
            project_dir, "000_karioCaS_input_matrix", "karioCaS_TSE.rds"
        )
        if (!file.exists(tse_path)) stop("TSE file missing for enrichment.")
        tse <- readRDS(tse_path)
        tax_df <- SummarizedExperiment::rowData(tse) |> as.data.frame()
        if (!"Taxonomy_Full" %in% colnames(df_long)) {
            if ("Taxonomy" %in% colnames(df_long)) {
                df_long <- dplyr::rename(df_long, Taxonomy_Full = "Taxonomy")
            } else {
                stop("Could not find Taxonomy key column.")
            }
        }
        cols_to_add <- setdiff(colnames(tax_df), colnames(df_long))
        tax_sub <- dplyr::select(
            tax_df, "Taxonomy_Full", dplyr::all_of(cols_to_add)
        )
        df_long <- dplyr::left_join(df_long, tax_sub, by = "Taxonomy_Full")
    }
    df_long |>
        dplyr::mutate(
            Parent_Name = .data[[parent_level]],
            Child_Name  = .data[[child_level]],
            Rank        = trimws(.data$Rank)
        )
}

#' @noRd
.txr_calc_resolution <- function(df_dom, parent_level, child_level, top_n) {
    parent_totals <- df_dom |>
        dplyr::filter(.data$Rank == parent_level) |>
        dplyr::group_by(.data$Parent_Name) |>
        dplyr::summarise(
            Total_Clade_Reads = .safe_max(.data$Counts),
            .groups = "drop"
        ) |>
        dplyr::filter(!is.na(.data$Parent_Name)) |>
        dplyr::arrange(dplyr::desc(.data$Total_Clade_Reads)) |>
        dplyr::slice_head(n = top_n)
    top_parents <- parent_totals$Parent_Name
    if (length(top_parents) == 0) {
        return(NULL)
    }
    df_parent <- df_dom |>
        dplyr::filter(
            .data$Parent_Name %in% top_parents,
            .data$Rank == parent_level
        ) |>
        dplyr::group_by(.data$Parent_Name) |>
        dplyr::summarise(
            Parent_Cumulative_Total = .safe_max(.data$Counts),
            .groups = "drop"
        )
    df_child <- df_dom |>
        dplyr::filter(
            .data$Parent_Name %in% top_parents,
            .data$Rank == child_level
        ) |>
        dplyr::group_by(.data$Parent_Name) |>
        dplyr::summarise(
            Child_Sum_Resolved = sum(.data$Counts, na.rm = TRUE),
            .groups = "drop"
        )
    df_parent |>
        dplyr::left_join(df_child, by = "Parent_Name") |>
        dplyr::mutate(
            Child_Sum_Resolved = tidyr::replace_na(.data$Child_Sum_Resolved, 0),
            Parent_Exclusive = pmax(
                0, .data$Parent_Cumulative_Total - .data$Child_Sum_Resolved
            )
        ) |>
        dplyr::arrange(dplyr::desc(.data$Parent_Cumulative_Total))
}

#' @noRd
.txr_audit_log <- function(calc_df, log_msg) {
    # Report the top parent clade (by cumulative reads) as a sanity check.
    audit <- calc_df[1, ]
    log_msg(sprintf(
        "    AUDIT [%s]: Parent(Max)=%d | Child(Sum)=%d | Exclusive(Diff)=%d",
        audit$Parent_Name[1],
        audit$Parent_Cumulative_Total[1],
        audit$Child_Sum_Resolved[1],
        audit$Parent_Exclusive[1]
    ))
}

#' @noRd
.txr_build_stack <- function(calc_df, top_n, parent_level, child_level) {
    label_excl <- paste0(parent_level, "-exclusive")
    label_chld <- child_level
    top_parents <- calc_df$Parent_Name
    stack <- calc_df |>
        dplyr::select("Parent_Name", "Child_Sum_Resolved", "Parent_Exclusive") |>
        tidyr::pivot_longer(
            cols      = c("Child_Sum_Resolved", "Parent_Exclusive"),
            names_to  = "Category_Raw",
            values_to = "Reads"
        ) |>
        dplyr::mutate(Category = dplyr::case_when(
            .data$Category_Raw == "Child_Sum_Resolved" ~ label_chld,
            .data$Category_Raw == "Parent_Exclusive" ~ label_excl
        ))
    missing <- top_n - length(top_parents)
    if (missing > 0) {
        ghost_df <- data.frame(
            Parent_Name  = paste0("Ghost_", seq_len(missing)),
            Category_Raw = label_excl,
            Reads        = 0,
            Category     = label_excl
        )
        stack <- dplyr::bind_rows(stack, ghost_df)
        top_parents <- c(top_parents, ghost_df$Parent_Name)
    }
    stack$Parent_Name <- factor(stack$Parent_Name, levels = rev(top_parents))
    stack$Category <- factor(stack$Category, levels = c(label_excl, label_chld))
    list(stack = stack, label_excl = label_excl, label_chld = label_chld)
}

#' @noRd
.txr_domain_plot <- function(df_dom, dom, parent_level, child_level,
                             top_n, log_msg) {
    calc_df <- .txr_calc_resolution(df_dom, parent_level, child_level, top_n)
    if (is.null(calc_df)) {
        return(plot_kariocas_empty(dom, "No Data"))
    }
    .txr_audit_log(calc_df, log_msg)
    res <- .txr_build_stack(calc_df, top_n, parent_level, child_level)
    spec_colors <- get_kariocas_colors("special")
    fill_values <- stats::setNames(
        c(spec_colors[["Parent"]], spec_colors[["Child"]]),
        c(res$label_excl, res$label_chld)
    )
    clean_labels <- function(x) ifelse(grepl("^Ghost_", x), "", x)
    ggplot2::ggplot(
        res$stack,
        ggplot2::aes(x = .data$Reads, y = .data$Parent_Name, fill = .data$Category)
    ) +
        ggplot2::geom_col(width = 0.7) +
        ggplot2::scale_fill_manual(values = fill_values) +
        ggplot2::scale_x_continuous(labels = label_k_number) +
        ggplot2::scale_y_discrete(labels = clean_labels) +
        ggplot2::labs(title = dom, x = "Reads", y = NULL) +
        theme_kariocas() +
        ggplot2::theme(
            axis.text.y        = ggplot2::element_text(face = "italic"),
            panel.grid.major.y = ggplot2::element_blank(),
            legend.position    = "bottom",
            legend.title       = ggplot2::element_blank()
        )
}

#' @noRd
.txr_save_panel <- function(plots, DOMAINS, samp, cs,
                            parent_level, child_level, output_dir, log_msg) {
    layout <- (plots[["Bacteria"]] | plots[["Archaea"]]) /
        (plots[["Eukaryota"]] | plots[["Viruses"]]) +
        patchwork::plot_annotation(
            title = paste(samp, "- CS", sprintf("%02d", cs)),
            subtitle = paste("Taxa Resolution:", parent_level, "vs", child_level),
            theme = ggplot2::theme(
                plot.title = ggplot2::element_text(
                    face = "bold", size = 16, hjust = 0.5
                ),
                plot.subtitle = ggplot2::element_text(
                    size = 12, hjust = 0.5, color = "grey30"
                )
            )
        ) +
        patchwork::plot_layout(guides = "collect") &
        ggplot2::theme(legend.position = "bottom")
    fname <- paste0(
        samp, "_CS", sprintf("%02d", cs),
        "_Resolution_", parent_level, "_vs_", child_level, ".pdf"
    )
    ggplot2::ggsave(
        file.path(output_dir, fname), layout,
        width = get_kariocas_dims()$width,
        height = get_kariocas_dims()$height
    )
    log_msg("    -> Generated: ", fname)
}

# ==============================================================================
# EXPORTED FUNCTION
# ==============================================================================

#' Generate Taxa Resolution Analysis (Step 004)
#'
#' Creates stacked bar plots showing resolution efficiency between a Parent Rank
#' and a Child Rank. Uses \code{max()} for parent aggregation to prevent
#' double-counting of children in cumulative data.
#'
#' @param project_dir Path to the project root.
#' @param parent_level Name of the parent rank (default: \code{"Genus"}).
#' @param child_level Name of the child rank (default: \code{"Species"}).
#' @param top_n Number of top taxa to display per domain (default: 10).
#'
#' @return Invisibly returns \code{NULL}. PDF plots are saved to
#'   \code{<project_dir>/004_taxa_resolution/}.
#' @export
#' @importFrom dplyr filter mutate select group_by summarise arrange slice_head
#'   left_join bind_rows distinct case_when pull rename all_of desc
#' @importFrom tidyr pivot_longer replace_na
#' @importFrom ggplot2 ggplot aes geom_col scale_fill_manual scale_x_continuous
#'   scale_y_discrete labs theme element_text guides guide_legend ggsave
#' @importFrom patchwork plot_layout plot_annotation
#' @importFrom scales label_number
#' @importFrom ggtext element_markdown
#' @importFrom SummarizedExperiment rowData
#' @examples
#' toy_project <- system.file("extdata", "your_project_name", package = "karioCaS")
#'
#' # Basic usage (Genus vs Species, top 10)
#' # taxa_resolution(project_dir = toy_project)
#'
#' # Family to Genus, top 5
#' # taxa_resolution(
#' #   project_dir  = toy_project,
#' #   parent_level = "Family",
#' #   child_level  = "Genus",
#' #   top_n        = 5
#' # )
taxa_resolution <- function(project_dir,
                            parent_level = "Genus",
                            child_level = "Species",
                            top_n = 10) {
    setup <- .txr_setup(project_dir, parent_level, child_level)
    df_proc <- .txr_load_and_enrich(
        project_dir, parent_level, child_level, setup$log_msg
    )
    SAMPLES <- unique(df_proc$sample)
    DOMAINS <- names(get_kariocas_colors("domains"))
    CS_LIST <- unique(df_proc$CS)
    setup$log_msg(">>> Starting Analysis Loop...")
    for (samp in SAMPLES) {
        setup$log_msg("------------------------------------------------")
        setup$log_msg("  Processing Sample: ", samp)
        df_samp <- dplyr::filter(df_proc, .data$sample == samp)
        for (cs in CS_LIST) {
            df_curr <- dplyr::filter(df_samp, .data$CS == cs)
            if (nrow(df_curr) == 0) next
            plots <- stats::setNames(
                lapply(DOMAINS, function(dom) {
                    df_dom <- dplyr::filter(df_curr, .data$Domain == dom)
                    .txr_domain_plot(
                        df_dom, dom, parent_level, child_level,
                        top_n, setup$log_msg
                    )
                }),
                DOMAINS
            )
            .txr_save_panel(
                plots, DOMAINS, samp, cs, parent_level,
                child_level, setup$output_dir, setup$log_msg
            )
        }
    }
    setup$log_msg("SUCCESS: Resolution analysis completed.")
    invisible(NULL)
}
