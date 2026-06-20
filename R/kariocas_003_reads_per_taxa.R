# ==============================================================================
# PRIVATE HELPERS - reads_per_taxa()
# ==============================================================================

#' @noRd
.rpt_setup <- function(project_dir, analysis_level, method) {
    output_dir <- file.path(project_dir, "003_cutoffs")
    log_dir <- file.path(project_dir, "logs")
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)
    log_file <- file.path(log_dir, "log_003_cutoff_analysis.txt")
    writeLines(c(
        "====================================================",
        "LOG: 003_CUTOFF_ANALYSIS (Saturation + Optimal Reads)",
        paste0("PROJECT DIR: ", project_dir),
        paste0("ANALYSIS LEVEL: ", analysis_level),
        paste0("METHOD: ", toupper(method)),
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
.rpt_cutoff_template <- function() {
    base_steps <- c(1, 3, 5)
    multipliers <- c(1, 10, 100, 1000, 10000, 100000, 1000000)
    sort(unique(as.vector(outer(base_steps, multipliers, "*"))))
}

#' Read-count cutoffs relevant to one sample/domain (truncated to its range).
#' @noRd
.rpt_sample_cutoffs <- function(max_count, cutoffs_sat) {
    ct <- cutoffs_sat[cutoffs_sat <= max_count]
    nxt <- cutoffs_sat[cutoffs_sat > max_count][1]
    if (!is.na(nxt)) c(ct, nxt) else ct
}

# ------------------------------------------------------------------------------
# Group overlay + optimal-reads audit
# ------------------------------------------------------------------------------

#' Per-sample saturation curves + per-sample optimal-reads audit, one domain.
#' @return list(overlay = long df for plotting, audit = tagged SI rows).
#' @noRd
.rpt_cs_domain <- function(df_cs_dom, cutoffs_sat, method, cs, dom) {
    samples <- unique(df_cs_dom$sample)
    overlay <- list()
    audit <- list()
    for (s in samples) {
        df_s <- dplyr::filter(df_cs_dom, .data$sample == s)
        if (nrow(df_s) == 0) next
        cutoffs <- .rpt_sample_cutoffs(max(df_s$Counts), cutoffs_sat)
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
                                    analysis_level, DOMAINS, output_dir, log_msg) {
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
        fname, output_dir, log_msg
    )
}

#' @noRd
.rpt_group_analysis <- function(df_proc, CS_LIST, DOMAINS, cutoffs_sat,
                                analysis_level, method, output_dir, log_msg) {
    audit_all <- list()
    for (grp in unique(df_proc$Group)) {
        df_grp <- dplyr::filter(df_proc, .data$Group == grp)
        n_samples <- dplyr::n_distinct(df_grp$sample)
        log_msg("  Group: ", grp, " (", n_samples, " samples)")
        for (cs in CS_LIST) {
            res <- lapply(DOMAINS, function(dom) {
                df_cs_dom <- dplyr::filter(
                    df_grp, .data$CS == cs, .data$Domain == dom, .data$Counts > 0
                )
                .rpt_cs_domain(df_cs_dom, cutoffs_sat, method, cs, dom)
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
            .rpt_save_group_overlay(
                overlay, vlines, grp, cs, n_samples,
                analysis_level, DOMAINS, output_dir, log_msg
            )
        }
    }
    dplyr::bind_rows(audit_all)
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
.rpt_detail_domain_plot <- function(df_curr, dom, analysis_level, cutoffs_sat) {
    df_dom <- dplyr::filter(df_curr, .data$Domain == dom, .data$Counts > 0)
    x_lab <- get_kariocas_labels()$x_log10_reads
    if (nrow(df_dom) == 0) {
        return(plot_kariocas_empty(dom, "No reads detected", x_lab, "**% Retained**"))
    }
    total_reads <- sum(df_dom$Counts)
    total_taxa <- nrow(df_dom)
    max_val <- max(df_dom$Counts)
    cutoffs <- .rpt_sample_cutoffs(max_val, cutoffs_sat)
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
                           cutoffs_sat, output_dir, log_msg) {
    df_curr <- dplyr::filter(df_proc, .data$sample == samp, .data$CS == cs)
    if (nrow(df_curr) == 0) {
        return(invisible(NULL))
    }
    plots <- stats::setNames(
        lapply(DOMAINS, function(dom) {
            .rpt_detail_domain_plot(df_curr, dom, analysis_level, cutoffs_sat)
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
    ggplot2::ggsave(
        file.path(output_dir, fname), layout,
        width = get_kariocas_dims()$width, height = get_kariocas_dims()$height
    )
    log_msg("    -> Saved: ", fname)
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
#' (Step 001).
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
#'
#' @return Invisibly returns a \code{data.frame} with the optimal-reads audit.
#'   PDF plots and \code{Reads_Audit_<rank>} files are saved to
#'   \code{<project_dir>/003_cutoffs/}.
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
#' toy_project <- system.file("extdata", "your_project_name", package = "karioCaS")
#'
#' # Group saturation overlays + optimal minimum reads (default Kneedle)
#' # reads_per_taxa(project_dir = toy_project)
reads_per_taxa <- function(project_dir,
                           analysis_level = "Species",
                           method = c("kneedle", "postcliff", "segmented"),
                           detail_samples = NULL) {
    method <- match.arg(method)
    setup <- .rpt_setup(project_dir, analysis_level, method)
    setup$log_msg(">>> Loading Data (Auto-detected format)...")
    df_long <- .get_tidy_data(project_dir)
    df_proc <- dplyr::filter(df_long, .data$Rank == analysis_level)
    if (nrow(df_proc) == 0) {
        setup$log_msg("CRITICAL ERROR: No data found for Rank: ", analysis_level)
        stop("No data for specified rank.")
    }
    df_proc$Group <- .grp_parse_group(df_proc$sample)
    SAMPLES <- unique(df_proc$sample)
    DOMAINS <- names(get_kariocas_colors("domains"))
    CS_LIST <- sort(unique(df_proc$CS))
    cutoffs_sat <- .rpt_cutoff_template()

    setup$log_msg(
        ">>> Saturation + optimal reads (method: ", method, ")..."
    )
    full_audit <- .rpt_group_analysis(
        df_proc, CS_LIST, DOMAINS, cutoffs_sat,
        analysis_level, method, setup$output_dir, setup$log_msg
    )
    if (!is.null(full_audit) && nrow(full_audit) > 0) {
        tsv_path <- file.path(
            setup$output_dir, paste0("Reads_Audit_", analysis_level, ".tsv")
        )
        rds_path <- file.path(
            setup$output_dir, paste0("Reads_Audit_", analysis_level, ".rds")
        )
        readr::write_tsv(full_audit, tsv_path)
        readr::write_rds(full_audit, rds_path)
        setup$log_msg("SAVED READS AUDIT TSV: ", tsv_path)
        setup$log_msg("SAVED READS AUDIT RDS: ", rds_path)
    }

    detail <- .grp_resolve_detail(detail_samples, SAMPLES, setup$log_msg)
    if (length(detail) > 0) {
        detail_dir <- file.path(setup$output_dir, "per_sample")
        if (!dir.exists(detail_dir)) dir.create(detail_dir, recursive = TRUE)
        setup$log_msg(
            ">>> Rendering detailed panels for ", length(detail), " sample(s)."
        )
        for (samp in detail) {
            for (cs in CS_LIST) {
                .rpt_detail_cs(
                    df_proc, samp, cs, DOMAINS, analysis_level,
                    cutoffs_sat, detail_dir, setup$log_msg
                )
            }
        }
    }
    setup$log_msg("SUCCESS: Cutoff analysis completed.")
    invisible(full_audit)
}
