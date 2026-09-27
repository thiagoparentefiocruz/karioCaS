# ==============================================================================
# PRIVATE HELPERS - reads_per_taxa()
# ==============================================================================

#' @noRd
.rpt_setup <- function(project_dir, analysis_level, method, export = TRUE) {
    .kcs_setup_step(
        project_dir, "003_reads_saturation", "log_003_reads_saturation.txt",
        export = export,
        extra_header = c(
            paste0("ANALYSIS LEVEL: ", analysis_level),
            paste0("METHOD: ", toupper(method))
        )
    )
}

#' Adaptive read-count cutoffs for a saturation curve.
#'
#' Always includes the low-end anchors \code{1, 2, 3, 4, 5, 7, 10} (even when the
#' maximum read count is smaller), then adds log-spaced points
#' (\code{1, 2, 3, 5, 7} x 10^k) up to and just past \code{max_count}. This gives
#' fine resolution at the low end (where rare/background taxa are shed) while
#' adapting the upper range to the actual data instead of a hard-coded ceiling.
#' @param max_count Maximum read count in the sample/domain.
#' @return A sorted numeric vector of cutoffs.
#' @noRd
.rpt_cutoffs_for <- function(max_count) {
    anchors <- c(1, 2, 3, 4, 5, 7, 10)
    if (!is.finite(max_count) || max_count <= 10) {
        return(anchors)
    }
    mult <- c(1, 2, 3, 5, 7)
    decades <- 10^seq_len(ceiling(log10(max_count)))
    hi <- sort(unique(as.vector(outer(mult, decades))))
    hi <- hi[hi > 10]
    keep <- hi[hi <= max_count]
    nxt <- hi[hi > max_count][1]
    if (!is.na(nxt)) keep <- c(keep, nxt)
    sort(unique(c(anchors, keep)))
}

# ------------------------------------------------------------------------------
# Group overlay + optimal-reads audit
# ------------------------------------------------------------------------------

#' Per-sample saturation curves + per-sample optimal-reads audit, one domain.
#' @return list(overlay = long df for plotting, audit = tagged SI rows).
#' @noRd
.rpt_cs_domain <- function(df_cs_dom, method, cs, dom) {
    samples <- unique(df_cs_dom$sample)
    overlay <- list()
    audit <- list()
    for (s in samples) {
        df_s <- dplyr::filter(df_cs_dom, .data$sample == s)
        if (nrow(df_s) == 0) next
        cutoffs <- .rpt_cutoffs_for(max(df_s$Counts))
        if (length(cutoffs) == 0) next
        total <- nrow(df_s)
        surv <- vapply(cutoffs, function(k) sum(df_s$Counts >= k), numeric(1))
        overlay[[s]] <- data.frame(
            sample = s, Domain = dom, x = cutoffs, y = surv / total * 100
        )
        el <- .si_reads_elbow(data.frame(Cutoff = cutoffs, Taxa_Count = surv), method)
        if (!is.null(el)) {
            a <- el$calc
            a$Sample <- s
            a$CS <- cs
            a$Domain <- dom
            a$SI_Type <- dplyr::case_when(
                a$Cutoff == el$opt ~ "Primary_SI",
                !is.na(el$sec) & a$Cutoff == el$sec ~ "Secondary_SI_1",
                TRUE ~ NA_character_
            )
            audit[[s]] <- a[, c(
                "Sample", "CS", "Domain", "Cutoff",
                "Taxa_Count", "Pct_Retained", "Step_Loss_Pct", "SI_Type"
            )]
        }
    }
    list(
        overlay = dplyr::bind_rows(overlay),
        audit = dplyr::bind_rows(audit)
    )
}

#' @noRd
.rpt_save_group_overlay <- function(overlay, vlines, grp, cs, n_samples,
                                    analysis_level, DOMAINS, setup) {
    apply_scales <- function(p) {
        p +
            scale_x_kariocas_log10(labels = label_k_number) +
            ggplot2::coord_cartesian(ylim = c(0, 105))
    }
    plots <- .grp_overlay_plots(
        overlay, DOMAINS, get_kariocas_labels()$x_log10_reads,
        "**% Retained**", apply_scales,
        vlines = vlines
    )
    fname <- paste0(
        grp, "_Group_CS", sprintf("%02d", cs),
        "_", analysis_level, "_Saturation.pdf"
    )
    .grp_assemble_2x2(
        plots,
        paste0(grp, " - CS", sprintf("%02d", cs), " | Saturation"),
        paste0(
            "n = ", n_samples, " samples (", analysis_level,
            ") | bold = group mean, dashed = median optimal reads"
        ),
        fname, setup
    )
}

#' @noRd
.rpt_group_analysis <- function(df_proc, CS_LIST, DOMAINS,
                                analysis_level, method, setup) {
    log_msg <- setup$log_msg
    audit_all <- list()
    plots <- list()
    paths <- character(0)
    for (grp in unique(df_proc$Group)) {
        df_grp <- dplyr::filter(df_proc, .data$Group == grp)
        n_samples <- dplyr::n_distinct(df_grp$sample)
        log_msg("  Group: ", grp, " (", n_samples, " samples)")
        for (cs in CS_LIST) {
            res <- lapply(DOMAINS, function(dom) {
                df_cs_dom <- dplyr::filter(
                    df_grp, .data$CS == cs, .data$Domain == dom, .data$Counts > 0
                )
                .rpt_cs_domain(df_cs_dom, method, cs, dom)
            })
            overlay <- dplyr::bind_rows(lapply(res, `[[`, "overlay"))
            audit <- dplyr::bind_rows(lapply(res, `[[`, "audit"))
            if (nrow(audit) > 0) audit_all[[length(audit_all) + 1]] <- audit
            if (nrow(overlay) == 0) next
            vlines <- NULL
            if (nrow(audit) > 0) {
                prim <- audit |>
                    dplyr::filter(.data$SI_Type == "Primary_SI") |>
                    dplyr::group_by(.data$Domain) |>
                    dplyr::summarise(
                        v = stats::median(.data$Cutoff), .groups = "drop"
                    )
                vlines <- stats::setNames(prim$v, prim$Domain)
            }
            out <- .rpt_save_group_overlay(
                overlay, vlines, grp, cs, n_samples,
                analysis_level, DOMAINS, setup
            )
            plots[[paste0(grp, "_CS", sprintf("%02d", cs), "_Saturation")]] <- out$plot
            paths <- c(paths, out$path)
        }
    }
    list(data = dplyr::bind_rows(audit_all), plots = plots, paths = paths)
}

# ------------------------------------------------------------------------------
# Per-sample detail (saturation curve: taxa + reads)
# ------------------------------------------------------------------------------

#' @noRd
.rpt_calc_stats <- function(df_dom, relevant_cutoffs, total_reads, total_taxa) {
    stats_list <- lapply(relevant_cutoffs, function(k) {
        survivors <- dplyr::filter(df_dom, .data$Counts >= k)
        data.frame(
            Cutoff        = k,
            Ret_Reads_Pct = sum(survivors$Counts) / total_reads,
            Ret_Taxa_Pct  = nrow(survivors) / total_taxa
        )
    })
    do.call(rbind, stats_list)
}

#' @noRd
.rpt_detail_domain_plot <- function(df_curr, dom, analysis_level) {
    df_dom <- dplyr::filter(df_curr, .data$Domain == dom, .data$Counts > 0)
    x_lab <- get_kariocas_labels()$x_log10_reads
    if (nrow(df_dom) == 0) {
        return(plot_kariocas_empty(dom, "No reads detected", x_lab, "**% Retained**"))
    }
    total_reads <- sum(df_dom$Counts)
    total_taxa <- nrow(df_dom)
    max_val <- max(df_dom$Counts)
    cutoffs <- .rpt_cutoffs_for(max_val)
    df_stats <- .rpt_calc_stats(df_dom, cutoffs, total_reads, total_taxa)
    df_plot <- df_stats |>
        tidyr::pivot_longer(
            cols = c("Ret_Reads_Pct", "Ret_Taxa_Pct"),
            names_to = "Metric", values_to = "Pct"
        ) |>
        dplyr::mutate(
            Metric_Key = factor(
                dplyr::case_when(
                    .data$Metric == "Ret_Reads_Pct" ~ "Level Reads",
                    .data$Metric == "Ret_Taxa_Pct" ~ "Level Taxa"
                ),
                levels = c("Level Taxa", "Level Reads")
            )
        )
    legend_labels <- c(analysis_level, "Reads")
    ggplot2::ggplot(
        df_plot,
        ggplot2::aes(x = .data$Cutoff, y = .data$Pct, group = .data$Metric_Key)
    ) +
        ggplot2::geom_line(
            ggplot2::aes(color = .data$Metric_Key, linetype = .data$Metric_Key),
            linewidth = 1
        ) +
        ggplot2::geom_point(
            ggplot2::aes(color = .data$Metric_Key, shape = .data$Metric_Key),
            size = 3
        ) +
        ggplot2::scale_y_continuous(
            breaks = seq(0, 1, 0.25), labels = function(x) x * 100,
            limits = c(0, 1.05)
        ) +
        ggplot2::scale_color_manual(
            values = get_kariocas_colors("special"), labels = legend_labels
        ) +
        ggplot2::scale_shape_manual(
            values = get_kariocas_shapes("ranks"), labels = legend_labels
        ) +
        ggplot2::scale_linetype_manual(
            values = get_kariocas_linetypes(), labels = legend_labels
        ) +
        scale_x_kariocas_log10(labels = label_k_number) +
        ggplot2::labs(
            title = dom,
            subtitle = paste0(
                label_kariocas_auto(total_reads), " Reads | ",
                label_kariocas_auto(total_taxa), " ", analysis_level
            ),
            x = x_lab, y = "**% Retained**"
        ) +
        theme_kariocas() +
        ggplot2::guides(
            color = ggplot2::guide_legend(nrow = 1, title = NULL),
            shape = ggplot2::guide_legend(nrow = 1, title = NULL),
            linetype = ggplot2::guide_legend(nrow = 1, title = NULL)
        )
}

#' @noRd
.rpt_detail_cs <- function(df_proc, samp, cs, DOMAINS, analysis_level,
                           setup, dir) {
    df_curr <- dplyr::filter(df_proc, .data$sample == samp, .data$CS == cs)
    if (nrow(df_curr) == 0) {
        return(NULL)
    }
    plots <- stats::setNames(
        lapply(DOMAINS, function(dom) {
            .rpt_detail_domain_plot(df_curr, dom, analysis_level)
        }),
        DOMAINS
    )
    layout <- (plots[["Bacteria"]] | plots[["Archaea"]]) /
        (plots[["Eukaryota"]] | plots[["Viruses"]]) +
        patchwork::plot_annotation(
            title = paste(samp, "- CS", sprintf("%02d", cs), "| Saturation"),
            subtitle = paste("Retention of Reads vs", analysis_level),
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
        "_Cutoff_", analysis_level, "_Saturation.pdf"
    )
    list(plot = layout, path = .kcs_save_plot(layout, fname, setup, dir = dir))
}

# ==============================================================================
# EXPORTED FUNCTION
# ==============================================================================

#' Read Cutoff Saturation Analysis & Optimal Minimum Reads (Step 003)
#'
#' Saturation analysis: progressively raises a per-taxon read-count cutoff and
#' tracks how many taxa survive, on a log read axis. By default it draws one
#' \strong{group overlay} per biological group (every sample a faint line, group
#' mean highlighted) and marks each domain's \strong{median optimal minimum
#' reads} - the elbow of the saturation curve, found with the same engine used
#' for the optimal CS - as a dashed line. The per-sample optimal-reads values are
#' written to \code{Reads_Audit_<rank>.tsv}/\code{.rds}, giving a quantitative
#' threshold for excluding low-abundance background/false-positive taxa.
#'
#' The optimal reads is computed \emph{per Confidence Score}, since the
#' saturation curve changes with CS; read it off at your chosen optimal CS
#' (Step 002).
#'
#' @param project_dir Path to the project root.
#' @param analysis_level Taxonomic rank to analyze (default: \code{"Species"}).
#' @param method Elbow strategy for the optimal reads. One of \code{"kneedle"}
#'   (default), \code{"postcliff"} or \code{"segmented"}.
#' @param detail_samples Which samples to also render as detailed per-sample
#'   saturation panels. \code{NULL} (default) writes only the group overlays;
#'   \code{"all"} renders every sample; a comma-separated string such as
#'   \code{"SAMPLE33, SAMPLE45"} (or a character vector) renders just those.
#'   Detailed PDFs are saved to a \code{per_sample/} subfolder.
#' @param export Logical (default \code{TRUE}). When \code{TRUE}, PDF plots,
#'   the \code{Reads_Audit_<rank>.tsv}/\code{.rds} files and a log are written
#'   to \code{<project_dir>/003_reads_saturation/}. When \code{FALSE}, nothing
#'   is written and all results are only returned. Note that
#'   \code{\link{retrieve_selected_taxa}} with \code{reads_min_* = "auto"}
#'   reads the exported audit file, so run this step with \code{export = TRUE}
#'   before it.
#'
#' @section How the quantities are calculated:
#' For each sample, CS and domain, the taxa at \code{analysis_level} are
#' counted at increasing minimum-read cutoffs (1, 2, 3, 4, 5, 7, 10, then
#' 1, 2, 3, 5 and 7 \eqn{\times 10^k} up to just past the largest count).
#' \code{Taxa_Count} is the number of taxa with at least \code{Cutoff} reads,
#' where a taxon's reads are the count of its own MPA row (which already
#' includes reads classified to lower ranks). \code{Pct_Retained} is
#' \code{Taxa_Count} as a percentage of the count at the lowest cutoff, and
#' \code{Step_Loss_Pct} the loss since the previous cutoff. The optimal
#' minimum reads (Primary SI) is the elbow of this curve, found with the same
#' engine as in \code{\link{taxa_retention}} on a \eqn{\log_{10}} read axis;
#' with \code{"kneedle"}, the more conservative post-cliff cutoff is reported
#' as Secondary SI when it is larger.
#'
#' @return A \code{\link{kariocas_result}} object with \code{$data}, the
#'   optimal-reads audit table (one row per sample, CS, domain and read cutoff),
#'   \code{$plots}, a named list of \code{patchwork} figures (one saturation
#'   overlay per group and CS, plus detailed per-sample panels when requested)
#'   and \code{$paths}, the files written when \code{export = TRUE}.
#' @export
#' @importFrom dplyr filter mutate group_by summarise bind_rows case_when
#'   n_distinct
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_line geom_point scale_color_manual
#'   scale_shape_manual scale_linetype_manual scale_y_continuous labs
#'   coord_cartesian guides guide_legend element_text ggsave theme
#' @importFrom patchwork plot_layout plot_annotation
#' @importFrom readr write_tsv write_rds
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
#' # Saturation overlays + optimal minimum reads (Kneedle), in memory only
#' res <- reads_per_taxa(toy_project, export = FALSE)
#' res
#' head(res$data)
#' res$plots[[1]]
#'
#' unlink(toy_project, recursive = TRUE)
reads_per_taxa <- function(project_dir,
                           analysis_level = "Species",
                           method = c("kneedle", "postcliff", "segmented"),
                           detail_samples = NULL,
                           export = TRUE) {
    method <- match.arg(method)
    setup <- .rpt_setup(project_dir, analysis_level, method, export)
    setup$log_msg(">>> Loading Data (Auto-detected format)...")
    df_long <- .get_tidy_data(project_dir)
    # One row per taxon: its own (cumulative) MPA row, not its descendants'.
    df_proc <- dplyr::filter(
        df_long,
        .data$Rank == analysis_level, .data$Lowest_Rank == analysis_level
    )
    if (nrow(df_proc) == 0) {
        setup$log_msg("CRITICAL ERROR: No data found for Rank: ", analysis_level)
        stop("No data for specified rank.")
    }
    df_proc$Group <- .grp_parse_group(df_proc$sample)
    SAMPLES <- unique(df_proc$sample)
    DOMAINS <- names(get_kariocas_colors("domains"))
    CS_LIST <- sort(unique(df_proc$CS))

    setup$log_msg(
        ">>> Saturation + optimal reads (method: ", method, ")..."
    )
    ga <- .rpt_group_analysis(
        df_proc, CS_LIST, DOMAINS, analysis_level, method, setup
    )
    full_audit <- ga$data
    plots <- ga$plots
    paths <- ga$paths
    if (!is.null(full_audit) && nrow(full_audit) > 0) {
        paths <- c(paths, .kcs_save_table(
            full_audit, paste0("Reads_Audit_", analysis_level), setup
        ))
    }

    detail <- .grp_resolve_detail(detail_samples, SAMPLES)
    if (length(detail) > 0) {
        detail_dir <- if (export) file.path(setup$output_dir, "per_sample") else NULL
        setup$log_msg(
            ">>> Rendering detailed panels for ", length(detail), " sample(s)."
        )
        for (samp in detail) {
            for (cs in CS_LIST) {
                out <- .rpt_detail_cs(
                    df_proc, samp, cs, DOMAINS, analysis_level, setup, detail_dir
                )
                if (is.null(out)) next
                plots[[paste0(samp, "_CS", sprintf("%02d", cs), "_Saturation")]] <- out$plot
                paths <- c(paths, out$path)
            }
        }
    }
    .kcs_finish_step(
        "003_reads_saturation", full_audit, plots, paths, setup, setup$log_msg,
        what = "saturation results"
    )
}
