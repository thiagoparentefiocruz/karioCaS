# ==============================================================================
# PRIVATE HELPERS - taxa_retention()
# ==============================================================================

#' @noRd
.tr_setup <- function(project_dir, export = TRUE) {
    .kcs_setup_step(
        project_dir, "002_taxa_retention", "log_002_taxa_retention.txt",
        export = export
    )
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
    # Rank_Reads: reads classified to the rank, i.e. the sum of the rank's own
    # MPA rows (already cumulative), never adding their descendants again.
    stats_df <- df_samp |>
        dplyr::mutate(
            Own_Row = .data$Lowest_Rank == as.character(.data$Rank)
        ) |>
        dplyr::group_by(.data$Domain, .data$Rank, .data$CS) |>
        dplyr::summarise(
            Rank_Reads = sum(.data$Counts[.data$Own_Row]),
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
                            fmt_num, setup, dir) {
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
    path <- .kcs_save_plot(layout, fname, setup, dir = dir)
    list(plots = stats::setNames(list(layout), "All_Levels"), paths = path)
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
.tr_save_plots_b <- function(df_calc, DOMAINS, samp, fmt_num, setup, dir) {
    rank_map <- c(
        "Phylum" = "Phyla", "Class" = "Classes", "Order" = "Orders",
        "Family" = "Families", "Genus" = "Genera", "Species" = "Species"
    )
    plots <- list()
    paths <- character(0)
    for (r in names(rank_map)) {
        leg_taxa <- r
        leg_reads <- paste0("Reads classified to ", r)
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
        plots[[rank_map[[r]]]] <- layout
        paths <- c(paths, .kcs_save_plot(layout, fname, setup, dir = dir))
    }
    list(plots = plots, paths = paths)
}

#' @noRd
.tr_process_sample <- function(df_proc, samp, DOMAINS, setup, dir) {
    fmt_num <- function(x) format(x, big.mark = ",", scientific = FALSE)
    log_msg <- setup$log_msg
    log_msg("------------------------------------------------")
    log_msg("  Processing Sample: ", samp)
    df_samp <- dplyr::filter(df_proc, .data$sample == samp)
    baseline <- .tr_baseline_stats(df_samp, log_msg, fmt_num)
    df_calc <- .tr_compute_pct(baseline$stats_df, baseline$baseline_df)
    a <- .tr_save_plot_a(
        df_calc, baseline$baseline_df, DOMAINS, samp, fmt_num, setup, dir
    )
    b <- .tr_save_plots_b(df_calc, DOMAINS, samp, fmt_num, setup, dir)
    plots <- c(a$plots, b$plots)
    names(plots) <- paste0(samp, "_", names(plots))
    list(plots = plots, paths = c(a$paths, b$paths))
}

#' @noRd
.tr_group_overlay <- function(full_audit, DOMAINS, tax_level, method, setup) {
    log_msg <- setup$log_msg
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
    plots <- list()
    paths <- character(0)
    for (grp in unique(df$Group)) {
        df_g <- dplyr::filter(df, .data$Group == grp)
        n_samples <- dplyr::n_distinct(df_g$sample)
        prim_g <- dplyr::filter(prim, .data$Domain %in% unique(df_g$Domain))
        vlines <- stats::setNames(prim_g$vline, prim_g$Domain)
        log_msg("  Group: ", grp, " (", n_samples, " samples)")
        plots_g <- .grp_overlay_plots(
            df_g, DOMAINS, lbls$y_confidence, "**% Retained**",
            apply_scales,
            vlines = vlines
        )
        out <- .grp_assemble_2x2(
            plots_g,
            paste0(grp, " - Group Retention (", tax_level, ")"),
            paste0(
                "n = ", n_samples,
                " samples | bold = group mean, dashed = median optimal CS (",
                method, ")"
            ),
            paste0(grp, "_Group_Retention_", tax_level, ".pdf"),
            setup
        )
        plots[[paste0(grp, "_Group_Retention")]] <- out$plot
        paths <- c(paths, out$path)
    }
    list(plots = plots, paths = paths)
}

# ==============================================================================
# EXPORTED FUNCTION
# ==============================================================================

#' Run Confidence Score Retention & Optimization (Step 002)
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
#' @param export Logical (default \code{TRUE}). When \code{TRUE}, PDF plots,
#'   the \code{SI_Audit_<rank>.tsv}/\code{.rds} files and a log are written to
#'   \code{<project_dir>/002_taxa_retention/}. When \code{FALSE}, nothing is
#'   written to disk and all results are only returned. Note that
#'   \code{\link{retrieve_selected_taxa}} with \code{CS_* = "auto"} reads the
#'   exported audit file, so run this step with \code{export = TRUE} before it.
#'
#' @section How the quantities are calculated:
#' For each sample and domain, at each CS:
#' \itemize{
#'   \item \code{Taxa_Count}: number of distinct taxa at \code{tax_level}
#'     detected with at least one read.
#'   \item \code{Pct_Retained}: \code{Taxa_Count} as a percentage of its
#'     maximum across CS values (normally the value at CS 0).
#'   \item \code{Step_Loss_Pct}: percentage points of \code{Pct_Retained} lost
#'     since the previous CS.
#' }
#' The optimal CS (Primary SI) is then located on the \code{Pct_Retained}
#' curve. \code{"kneedle"}: the CS at which the curve lies furthest below the
#' straight line joining its first and last points (falling back to the
#' steepest single drop if the curve is not convex). \code{"postcliff"}: the
#' first CS, at or after the steepest drop, whose \code{Step_Loss_Pct} is at
#' most a tolerance \eqn{\max(\bar{t} + 1.5 s_t, 0.5)}, where \eqn{\bar{t}} and
#' \eqn{s_t} are the mean and SD of the step losses at CS \eqn{\ge} 50.
#' \code{"segmented"}: the breakpoint minimising the total residual sum of
#' squares of two straight-line fits. \code{"dynamic"}: the first CS above the
#' lowest one whose step loss is within the same tolerance. \code{"manual"}:
#' the first CS above the lowest one whose step loss is within
#' \code{manual_toll} (a single value, or a named list per domain). When no CS
#' qualifies, the highest CS is used. A Secondary SI, a stricter alternative,
#' is also reported by \code{"kneedle"} (the post-cliff value, when stricter
#' than the elbow), \code{"dynamic"} (the next CS whose step loss is at most
#' \eqn{\max(\bar{t}, 0.1)}) and \code{"manual"} (the next CS whose step loss
#' is at most 0.2). A domain needs at least three CS values and some loss of
#' taxa to be assessed.
#'
#' In the detailed per-sample panels, curves are expressed relative to CS 0:
#' taxa at each rank (distinct taxa detected), reads classified to that rank
#' (the sum of the rank's own MPA rows, which already include the reads
#' classified to lower ranks) and total reads of the domain (the domain row).
#'
#' @return A \code{\link{kariocas_result}} object with \code{$data}, the SI
#'   audit table (one row per sample, domain and CS, with the retention
#'   percentages and the Stability Index tags), \code{$plots}, a named list of
#'   \code{patchwork} figures (one group overlay per biological group, plus the
#'   detailed per-sample panels when requested) and \code{$paths}, the files
#'   written when \code{export = TRUE}.
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
#' # Copy the bundled toy project to a temporary folder and import it
#' toy_project <- file.path(tempdir(), "toy_karioCaS")
#' dir.create(toy_project, showWarnings = FALSE)
#' file.copy(
#'     system.file("extdata", "your_project_name", "000_mpa_original",
#'         package = "karioCaS"
#'     ),
#'     toy_project, recursive = TRUE
#' )
#' import_karioCaS(toy_project)
#'
#' # Retention curves + optimal CS (Kneedle), results kept in memory only
#' res <- taxa_retention(toy_project, export = FALSE)
#' res
#' head(res$data)
#' names(res$plots)
#' res$plots[[1]]
#'
#' unlink(toy_project, recursive = TRUE)
taxa_retention <- function(project_dir,
                           tax_level = "Species",
                           method = c(
                               "kneedle", "postcliff", "segmented",
                               "dynamic", "manual"
                           ),
                           manual_toll = 1.0,
                           detail_samples = NULL,
                           export = TRUE) {
    method <- match.arg(method)
    setup <- .tr_setup(project_dir, export)
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
    audit <- .si_export_audit(audit_list, tax_level, setup)
    full_audit <- audit$data
    plots <- list()
    paths <- audit$paths

    if (!is.null(full_audit) && nrow(full_audit) > 0) {
        setup$log_msg(">>> Building group overlay(s)...")
        ov <- .tr_group_overlay(full_audit, DOMAINS, tax_level, method, setup)
        plots <- c(plots, ov$plots)
        paths <- c(paths, ov$paths)
    } else {
        setup$log_msg("    [WARNING] No SI audit produced; skipping overlay.")
    }

    detail <- .grp_resolve_detail(detail_samples, SAMPLES)
    if (length(detail) > 0) {
        detail_dir <- if (export) file.path(setup$output_dir, "per_sample") else NULL
        setup$log_msg(
            ">>> Rendering detailed panels for ", length(detail), " sample(s)."
        )
        for (samp in detail) {
            ps <- .tr_process_sample(df_proc, samp, DOMAINS, setup, detail_dir)
            plots <- c(plots, ps$plots)
            paths <- c(paths, ps$paths)
        }
    }
    .kcs_finish_step(
        "002_taxa_retention", full_audit, plots, paths, setup, setup$log_msg,
        what = "retention results"
    )
}
