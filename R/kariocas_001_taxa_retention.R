# ==============================================================================
# PRIVATE HELPERS - taxa_retention()
# ==============================================================================

#' @noRd
.tr_setup <- function(project_dir) {
    output_dir <- file.path(project_dir, "001_taxa_retention")
    log_dir <- file.path(project_dir, "logs")
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)
    log_file <- file.path(log_dir, "log_001_cs_retention.txt")
    writeLines(c(
        "====================================================",
        "LOG: 001_CS_RETENTION_ANALYSIS",
        paste0("PROJECT DIR: ", project_dir),
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
    list(output_dir = output_dir, log_file = log_file, log_msg = log_msg)
}

#' @noRd
.tr_load_data <- function(project_dir, log_msg) {
    log_msg(">>> Loading Data (Auto-detected format)...")
    rank_levels <- c(
        "Domain", "Kingdom", "Phylum", "Class",
        "Order", "Family", "Genus", "Species"
    )
    df_long <- .get_tidy_data(project_dir)
    df_proc <- df_long |>
        dplyr::filter(.data$Rank %in% rank_levels) |>
        dplyr::mutate(Rank = factor(.data$Rank, levels = rank_levels))
    if (nrow(df_proc) == 0) {
        log_msg("CRITICAL ERROR: No data found for specified ranks.")
        stop("No data found.")
    }
    df_proc
}

#' @noRd
.tr_baseline_stats <- function(df_samp, log_msg, fmt_num) {
    domain_totals <- df_samp |>
        dplyr::filter(.data$Rank == "Domain", .data$Lowest_Rank == "Domain") |>
        dplyr::group_by(.data$Domain, .data$CS) |>
        dplyr::summarise(Global_Reads = max(.data$Counts), .groups = "drop")
    if (nrow(domain_totals) == 0) {
        domain_totals <- df_samp |>
            dplyr::filter(.data$Rank == "Domain") |>
            dplyr::group_by(.data$Domain, .data$CS) |>
            dplyr::summarise(Global_Reads = max(.data$Counts), .groups = "drop")
    }
    stats_df <- df_samp |>
        dplyr::group_by(.data$Domain, .data$Rank, .data$CS) |>
        dplyr::summarise(
            Rank_Reads = sum(.data$Counts),
            Rank_Taxa  = dplyr::n_distinct(.data$Taxon_Name),
            .groups    = "drop"
        ) |>
        dplyr::left_join(domain_totals, by = c("Domain", "CS"))
    baseline_df <- stats_df |>
        dplyr::filter(.data$CS == 0) |>
        dplyr::rename(
            Base_Global = "Global_Reads",
            Base_Reads  = "Rank_Reads",
            Base_Taxa   = "Rank_Taxa"
        ) |>
        dplyr::select("Domain", "Rank", "Base_Global", "Base_Reads", "Base_Taxa")
    qc <- dplyr::filter(baseline_df, .data$Rank == "Species")
    if (nrow(qc) > 0) {
        log_msg("    > QC STATS (CS00 Baseline):")
        for (i in seq_len(nrow(qc))) {
            log_msg(
                "      - ", qc$Domain[i], ": ",
                fmt_num(qc$Base_Global[i]), " Total Reads | ",
                fmt_num(qc$Base_Taxa[i]), " Distinct Species"
            )
        }
    }
    list(stats_df = stats_df, baseline_df = baseline_df)
}

#' @noRd
.tr_compute_pct <- function(stats_df, baseline_df) {
    stats_df |>
        dplyr::left_join(baseline_df, by = c("Domain", "Rank")) |>
        dplyr::mutate(
            Pct_Domain_Retained = ifelse(
                .data$Base_Global > 0,
                (.data$Global_Reads / .data$Base_Global) * 100, 0
            ),
            Pct_Rank_Reads = ifelse(
                .data$Base_Reads > 0,
                (.data$Rank_Reads / .data$Base_Reads) * 100, 0
            ),
            Pct_Taxa = ifelse(
                .data$Base_Taxa > 0,
                (.data$Rank_Taxa / .data$Base_Taxa) * 100, 0
            )
        )
}

#' @noRd
.tr_domain_plot_a <- function(df_calc, dom, baseline_df, fmt_num) {
    plot_ranks <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
    df_dom <- dplyr::filter(
        df_calc, .data$Domain == dom,
        .data$Rank %in% plot_ranks
    )
    if (nrow(df_dom) == 0) {
        return(plot_kariocas_empty(dom, "No Data"))
    }
    base <- dplyr::filter(baseline_df, .data$Domain == dom)
    dom_row <- dplyr::filter(base, .data$Rank == "Domain")
    total <- if (nrow(dom_row) > 0) dom_row$Base_Global[1] else 0
    cnt <- function(r) {
        v <- base$Base_Taxa[base$Rank == r]
        if (length(v) == 0) 0 else v
    }
    sub_str <- paste0(
        "Reads: ", fmt_num(total),
        " | P: ", cnt("Phylum"), " | C: ", cnt("Class"), " | O: ", cnt("Order"),
        " | F: ", cnt("Family"), " | G: ", cnt("Genus"), " | S: ", cnt("Species")
    )
    lbls <- get_kariocas_labels()
    ggplot2::ggplot(
        df_dom,
        ggplot2::aes(
            x = .data$CS, y = .data$Pct_Taxa,
            color = .data$Rank, group = .data$Rank, shape = .data$Rank
        )
    ) +
        ggplot2::geom_line(linewidth = 1) +
        ggplot2::geom_point(size = 2.5) +
        ggplot2::scale_color_manual(values = get_kariocas_colors("ranks")) +
        ggplot2::scale_shape_manual(values = get_kariocas_shapes("ranks")) +
        scale_y_kariocas_log10(limits = c(0.01, 105)) +
        ggplot2::scale_x_continuous(breaks = seq(0, 100, 20), limits = c(0, 100)) +
        ggplot2::labs(
            title = dom, subtitle = sub_str,
            x = lbls$y_confidence, y = lbls$y_log10_retained
        ) +
        theme_kariocas() +
        ggplot2::guides(color = ggplot2::guide_legend(nrow = 1))
}

#' @noRd
.tr_save_plot_a <- function(df_calc, baseline_df, DOMAINS, samp,
                            fmt_num, output_dir, log_msg) {
    plots <- stats::setNames(
        lapply(DOMAINS, function(d) {
            .tr_domain_plot_a(df_calc, d, baseline_df, fmt_num)
        }),
        DOMAINS
    )
    layout <- (plots[["Bacteria"]] | plots[["Archaea"]]) /
        (plots[["Eukaryota"]] | plots[["Viruses"]]) +
        patchwork::plot_annotation(
            title = paste(samp, "- Retention (All Levels)"),
            subtitle = "Comparison of taxa loss across ranks",
            theme = ggplot2::theme(
                plot.title = ggplot2::element_text(face = "bold", size = 16, hjust = 0.5)
            )
        ) +
        patchwork::plot_layout(guides = "collect") &
        ggplot2::theme(legend.position = "bottom")
    fname <- paste0(samp, "_CS_Retention_All_Levels.pdf")
    ggplot2::ggsave(
        file.path(output_dir, fname), layout,
        width = get_kariocas_dims()$width, height = get_kariocas_dims()$height
    )
    log_msg("    -> Generated: ", fname)
}

#' @noRd
.tr_prep_b_data <- function(df_calc, dom, r, leg_taxa, leg_reads, leg_total) {
    df_viz <- dplyr::filter(df_calc, .data$Domain == dom, .data$Rank == r)
    if (nrow(df_viz) == 0) {
        return(NULL)
    }
    df_viz |>
        dplyr::select("CS", "Pct_Taxa", "Pct_Rank_Reads", "Pct_Domain_Retained") |>
        tidyr::pivot_longer(
            cols      = c("Pct_Taxa", "Pct_Rank_Reads", "Pct_Domain_Retained"),
            names_to  = "Metric_Type",
            values_to = "Pct_Value"
        ) |>
        dplyr::mutate(
            Metric_Label = dplyr::case_when(
                .data$Metric_Type == "Pct_Taxa" ~ leg_taxa,
                .data$Metric_Type == "Pct_Rank_Reads" ~ leg_reads,
                .data$Metric_Type == "Pct_Domain_Retained" ~ leg_total
            ),
            Metric_Label = factor(
                .data$Metric_Label,
                levels = c(leg_taxa, leg_total, leg_reads)
            )
        )
}

#' @noRd
.tr_domain_plot_b <- function(df_calc, dom, r,
                              leg_taxa, leg_reads, leg_total, fmt_num) {
    df_long <- .tr_prep_b_data(df_calc, dom, r, leg_taxa, leg_reads, leg_total)
    if (is.null(df_long)) {
        return(plot_kariocas_empty(dom, "No Data"))
    }
    df_viz <- dplyr::filter(df_calc, .data$Domain == dom, .data$Rank == r)
    sub_str <- paste0(
        fmt_num(df_viz$Base_Global[1]), " Total Reads; ",
        fmt_num(df_viz$Base_Reads[1]), " assigned to ", r,
        " | ", fmt_num(df_viz$Base_Taxa[1]), " ", r
    )
    spec <- get_kariocas_colors("special")
    shps <- get_kariocas_shapes("ranks")
    lts <- get_kariocas_linetypes()
    lbls <- get_kariocas_labels()
    nms <- c(leg_taxa, leg_total, leg_reads)
    col_v <- setNames(c(spec[["Level Taxa"]], spec[["Total Reads"]], spec[["Level Reads"]]), nms)
    shp_v <- setNames(c(shps[["Level Taxa"]], shps[["Total Reads"]], shps[["Level Reads"]]), nms)
    lt_v <- setNames(c(lts[["Level Taxa"]], lts[["Total Reads"]], lts[["Level Reads"]]), nms)
    ggplot2::ggplot(
        df_long,
        ggplot2::aes(
            x = .data$CS, y = .data$Pct_Value,
            color = .data$Metric_Label, linetype = .data$Metric_Label,
            shape = .data$Metric_Label
        )
    ) +
        ggplot2::geom_line(linewidth = 1) +
        ggplot2::geom_point(size = 3) +
        ggplot2::scale_color_manual(values = col_v) +
        ggplot2::scale_linetype_manual(values = lt_v) +
        ggplot2::scale_shape_manual(values = shp_v) +
        scale_y_kariocas_log10(limits = c(0.01, 105)) +
        ggplot2::scale_x_continuous(breaks = seq(0, 100, 20), limits = c(0, 100)) +
        ggplot2::labs(
            title = dom, subtitle = sub_str,
            x = lbls$y_confidence, y = lbls$y_log10_retained,
            color = NULL, linetype = NULL, shape = NULL
        ) +
        theme_kariocas() +
        ggplot2::theme(legend.position = "bottom")
}

#' @noRd
.tr_save_plots_b <- function(df_calc, DOMAINS, samp, fmt_num, output_dir, log_msg) {
    rank_map <- c(
        "Phylum" = "Phyla", "Class" = "Classes", "Order" = "Orders",
        "Family" = "Families", "Genus" = "Genera", "Species" = "Species"
    )
    for (r in names(rank_map)) {
        leg_taxa <- r
        leg_reads <- paste0(r, "-exclusive Reads")
        leg_total <- "Total Reads"
        plots <- stats::setNames(
            lapply(DOMAINS, function(d) {
                .tr_domain_plot_b(df_calc, d, r, leg_taxa, leg_reads, leg_total, fmt_num)
            }),
            DOMAINS
        )
        layout <- (plots[["Bacteria"]] | plots[["Archaea"]]) /
            (plots[["Eukaryota"]] | plots[["Viruses"]]) +
            patchwork::plot_annotation(
                title = paste(samp, "- Retention:", rank_map[[r]]),
                subtitle = "Comparison of Taxa vs Reads Retention",
                theme = ggplot2::theme(
                    plot.title = ggplot2::element_text(face = "bold", size = 16, hjust = 0.5)
                )
            ) +
            patchwork::plot_layout(guides = "collect") &
            ggplot2::theme(legend.position = "bottom")
        fname <- paste0(samp, "_CS_Retention_", rank_map[[r]], ".pdf")
        ggplot2::ggsave(
            file.path(output_dir, fname), layout,
            width = get_kariocas_dims()$width, height = get_kariocas_dims()$height
        )
        log_msg("    -> Generated: ", fname)
    }
}

#' @noRd
.tr_process_sample <- function(df_proc, samp, DOMAINS, output_dir, log_msg) {
    fmt_num <- function(x) format(x, big.mark = ",", scientific = FALSE)
    log_msg("------------------------------------------------")
    log_msg("  Processing Sample: ", samp)
    df_samp <- dplyr::filter(df_proc, .data$sample == samp)
    baseline <- .tr_baseline_stats(df_samp, log_msg, fmt_num)
    df_calc <- .tr_compute_pct(baseline$stats_df, baseline$baseline_df)
    .tr_save_plot_a(
        df_calc, baseline$baseline_df, DOMAINS, samp,
        fmt_num, output_dir, log_msg
    )
    .tr_save_plots_b(df_calc, DOMAINS, samp, fmt_num, output_dir, log_msg)
}

#' @noRd
.tr_group_overlay <- function(full_audit, DOMAINS, tax_level, method,
                              output_dir, log_msg) {
    df <- full_audit
    df$Group <- .grp_parse_group(df$Sample)
    df$sample <- df$Sample
    df$x <- df$CS
    df$y <- df$Pct_Retained
    prim <- df |>
        dplyr::filter(.data$SI_Type == "Primary_SI") |>
        dplyr::group_by(.data$Domain) |>
        dplyr::summarise(
            vline = stats::median(.data$CS, na.rm = TRUE), .groups = "drop"
        )
    lbls <- get_kariocas_labels()
    apply_scales <- function(p) {
        p +
            ggplot2::scale_x_continuous(
                breaks = seq(0, 100, 20), limits = c(0, 100)
            ) +
            ggplot2::scale_y_continuous(limits = c(0, 105))
    }
    for (grp in unique(df$Group)) {
        df_g <- dplyr::filter(df, .data$Group == grp)
        n_samples <- dplyr::n_distinct(df_g$sample)
        prim_g <- dplyr::filter(prim, .data$Domain %in% unique(df_g$Domain))
        vlines <- stats::setNames(prim_g$vline, prim_g$Domain)
        log_msg("  Group: ", grp, " (", n_samples, " samples)")
        plots <- .grp_overlay_plots(
            df_g, DOMAINS, lbls$y_confidence, "**% Retained**",
            apply_scales, vlines = vlines
        )
        .grp_assemble_2x2(
            plots,
            paste0(grp, " - Group Retention (", tax_level, ")"),
            paste0(
                "n = ", n_samples,
                " samples | bold = group mean, dashed = median optimal CS (",
                method, ")"
            ),
            paste0(grp, "_Group_Retention_", tax_level, ".pdf"),
            output_dir, log_msg
        )
    }
}

# ==============================================================================
# EXPORTED FUNCTION
# ==============================================================================

#' Run Confidence Score Retention & Optimization (Step 001)
#'
#' Executes taxa retention analysis based on Confidence Score (Kraken/Bracken)
#' and, in the same step, computes the objective optimal CS (Stability Index, SI)
#' for each domain. By default it produces a single, low-clutter \strong{group
#' overlay} per biological group: every sample of the group is drawn as a faint
#' line with the group mean (\eqn{\pm}SD) highlighted and each domain's median
#' optimal CS marked with a dashed line, faceted by Domain. It also writes the
#' machine-readable SI audit (\code{SI_Audit_<rank>.tsv}/\code{.rds}) used by
#' \code{\link{retrieve_selected_taxa}}. Detailed per-sample panels (all ranks,
#' taxa vs reads) are written only on request.
#'
#' Groups are inferred from sample names by stripping trailing digits
#' (e.g. \code{SAMPLE33}, \code{SAMPLE34} both belong to group \code{SAMPLE}).
#'
#' The optimal CS is found with a multi-strategy engine: \code{"kneedle"}
#' (default, parameter-free elbow detection), \code{"postcliff"} (a more
#' conservative threshold past the steepest drop), \code{"segmented"}
#' (broken-stick regression), \code{"dynamic"}, or \code{"manual"}.
#'
#' @param project_dir Root path of the project.
#' @param tax_level Taxonomic rank used for the group overlay and SI
#'   (default: \code{"Species"}).
#' @param method Optimal-CS strategy. One of \code{"kneedle"} (default),
#'   \code{"postcliff"}, \code{"segmented"}, \code{"dynamic"}, or \code{"manual"}.
#' @param manual_toll Numeric or named list. Acceptable step-wise loss percentage,
#'   used only when \code{method = "manual"} (default: 1.0).
#' @param detail_samples Which samples to also render as detailed per-sample
#'   panels. \code{NULL} (default) writes only the group overlay; \code{"all"}
#'   renders every sample; a comma-separated string such as
#'   \code{"SAMPLE33, SAMPLE45"} (or a character vector) renders just those.
#'   Detailed PDFs are saved to a \code{per_sample/} subfolder.
#'
#' @return Invisibly returns a \code{data.frame} with the full SI audit trail.
#'   PDF plots and \code{SI_Audit_<rank>} files are saved to
#'   \code{<project_dir>/001_taxa_retention/}.
#' @export
#' @importFrom dplyr filter select group_by summarise mutate left_join arrange
#'   rename bind_rows pull distinct n_distinct case_when all_of lag
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_vline geom_ribbon
#'   scale_color_manual scale_linetype_manual scale_shape_manual labs ggsave
#'   scale_y_continuous scale_x_continuous theme guides guide_legend element_text
#' @importFrom patchwork plot_layout plot_annotation
#' @importFrom scales label_number
#' @importFrom stats median setNames
#' @importFrom ggtext element_markdown
#' @examples
#' toy_project <- system.file("extdata", "your_project_name", package = "karioCaS")
#'
#' # Group retention overlay + optimal CS (default Kneedle method)
#' # taxa_retention(project_dir = toy_project)
taxa_retention <- function(project_dir,
                           tax_level = "Species",
                           method = c(
                               "kneedle", "postcliff", "segmented",
                               "dynamic", "manual"
                           ),
                           manual_toll = 1.0,
                           detail_samples = NULL) {
    method <- match.arg(method)
    setup <- .tr_setup(project_dir)
    df_proc <- .tr_load_data(project_dir, setup$log_msg)
    DOMAINS <- names(get_kariocas_colors("domains"))
    SAMPLES <- unique(df_proc$sample)

    setup$log_msg(
        ">>> Computing Stability Index (method: ", method,
        ") at rank: ", tax_level
    )
    df_rank <- dplyr::filter(df_proc, .data$Rank == tax_level)
    audit_list <- list()
    for (samp in SAMPLES) {
        df_samp <- dplyr::filter(df_rank, .data$sample == samp)
        for (dom in DOMAINS) {
            a <- .si_domain_audit(
                df_samp, dom, method, manual_toll, samp, setup$log_msg
            )
            if (!is.null(a)) audit_list[[length(audit_list) + 1]] <- a
        }
    }
    full_audit <- .si_export_audit(
        audit_list, tax_level, setup$output_dir, setup$log_msg
    )

    if (!is.null(full_audit) && nrow(full_audit) > 0) {
        setup$log_msg(">>> Building group overlay(s)...")
        .tr_group_overlay(
            full_audit, DOMAINS, tax_level, method,
            setup$output_dir, setup$log_msg
        )
    } else {
        setup$log_msg("    [WARNING] No SI audit produced; skipping overlay.")
    }

    detail <- .grp_resolve_detail(detail_samples, SAMPLES, setup$log_msg)
    if (length(detail) > 0) {
        detail_dir <- file.path(setup$output_dir, "per_sample")
        if (!dir.exists(detail_dir)) dir.create(detail_dir, recursive = TRUE)
        setup$log_msg(
            ">>> Rendering detailed panels for ", length(detail), " sample(s)."
        )
        for (samp in detail) {
            .tr_process_sample(df_proc, samp, DOMAINS, detail_dir, setup$log_msg)
        }
    }
    setup$log_msg("SUCCESS: Retention analysis completed.")
    invisible(full_audit)
}
