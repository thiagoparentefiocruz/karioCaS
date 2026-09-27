# ==============================================================================
# PRIVATE HELPERS - upset_kariocas()
# ==============================================================================

#' @noRd
.ups_setup <- function(project_dir, export = TRUE) {
    .kcs_setup_step(
        project_dir, "005_taxa_intersections_across_CS",
        "log_005_taxa_intersections_cs.txt",
        export = export
    )
}

#' @noRd
.ups_binary_matrix <- function(df_long, samp, dom, lvl) {
    df_sub <- df_long |>
        dplyr::filter(
            .data$sample == samp,
            .data$Domain == dom,
            .data$Rank == lvl,
            .data$Counts > 0
        ) |>
        dplyr::select("CS", "Taxon_Name") |>
        dplyr::distinct()
    if (nrow(df_sub) == 0) {
        return(NULL)
    }
    df_sub |>
        dplyr::mutate(
            Present  = 1L,
            CS_Label = sprintf("CS%02d", as.numeric(.data$CS))
        ) |>
        dplyr::select(-"CS") |>
        tidyr::pivot_wider(
            names_from  = "CS_Label",
            values_from = "Present",
            values_fill = 0L
        ) |>
        as.data.frame()
}

#' Build the UpSet object for one sample/domain (does not draw it)
#' @noRd
.ups_build <- function(binary_matrix, upset_cols, lvl) {
    UpSetR::upset(
        binary_matrix,
        sets                = rev(upset_cols),
        keep.order          = TRUE,
        order.by            = "freq",
        empty.intersections = NULL,
        mainbar.y.label     = paste(lvl, "Intersections"),
        sets.x.label        = paste("Total", lvl, "per CS"),
        text.scale          = c(1.5, 1.5, 1.2, 1.2, 1.5, 1.3),
        mb.ratio            = c(0.6, 0.4),
        main.bar.color      = get_kariocas_colors("upset")$main,
        sets.bar.color      = get_kariocas_colors("upset")$sets,
        matrix.color        = "grey20",
        shade.color         = "grey90"
    )
}

#' @noRd
.ups_process_sample <- function(df_long, samp, DOMAINS, lvl, setup) {
    log_msg <- setup$log_msg
    log_msg("------------------------------------------------")
    log_msg("  Processing Sample: ", samp)
    samp_dir <- if (isTRUE(setup$export)) {
        file.path(setup$output_dir, samp)
    } else {
        NULL
    }
    data <- list()
    plots <- list()
    paths <- character(0)
    failed <- character(0)
    for (dom in DOMAINS) {
        mat <- .ups_binary_matrix(df_long, samp, dom, lvl)
        if (is.null(mat)) {
            log_msg("    Skipping ", dom, "-", lvl, ": No data found.")
            next
        }
        upset_cols <- setdiff(colnames(mat), "Taxon_Name")
        if (length(upset_cols) < 2) {
            log_msg(
                "    Skipping ", dom, "-", lvl,
                ": Not enough intersection levels (Found: ",
                length(upset_cols), ")"
            )
            next
        }
        key <- paste0(samp, "_", dom, "_", lvl)
        up <- tryCatch(
            .ups_build(mat, upset_cols, lvl),
            error = function(e) {
                failed <<- c(failed, paste0(key, ": ", conditionMessage(e)))
                NULL
            }
        )
        if (is.null(up)) next
        data[[dom]] <- data.frame(
            sample = samp, Domain = dom, Rank = lvl, mat,
            N_CS = rowSums(mat[, upset_cols, drop = FALSE]),
            check.names = FALSE, stringsAsFactors = FALSE
        )
        if (!is.null(samp_dir)) {
            path <- .kcs_draw_upset_pdf(
                up, file.path(samp_dir, paste0(key, "_UpSet.pdf")),
                paste(samp, "-", dom, "|", lvl, "Intersection Analysis"),
                log_msg
            )
            if (is.character(path) && !file.exists(path)) {
                failed <- c(failed, paste0(key, ": ", attr(path, "error")))
                next
            }
            paths <- c(paths, path)
        }
        plots[[key]] <- up
    }
    list(
        data = dplyr::bind_rows(data), plots = plots,
        paths = paths, failed = failed
    )
}

# ==============================================================================
# EXPORTED FUNCTION
# ==============================================================================

#' Generate UpSet Plots per Sample and Domain (Step 005)
#'
#' "karioCaS never are upset!" Generates UpSet plots showing taxon persistence
#' across Confidence Score levels, for a single taxonomic rank, with detailed
#' logging. One plot per sample and domain is written to a per-sample subfolder.
#'
#' @param project_dir Path to the project root.
#' @param tax_level Taxonomic rank to analyze (default: \code{"Species"}).
#' @param export Logical (default \code{TRUE}). When \code{TRUE}, one PDF per
#'   sample and domain is written to a per-sample subfolder of
#'   \code{<project_dir>/005_taxa_intersections_across_CS/}, together with a
#'   log. When \code{FALSE}, nothing is written and the UpSet plots are only
#'   returned.
#'
#' @return A \code{\link{kariocas_result}} object. \code{$data} is the
#'   presence/absence table behind the plots: one row per sample, domain and
#'   taxon, one \code{CSxx} column per Confidence Score (1 = detected with
#'   at least one read at that CS) and \code{N_CS}, the number of CS levels
#'   at which the taxon is detected. \code{$plots} holds one \code{UpSetR}
#'   plot per sample and domain; printing an element draws it on the current
#'   device. \code{$paths} lists the PDFs written when \code{export = TRUE}.
#'   If a plot cannot be drawn, a warning names it and it is not counted as a
#'   success.
#' @export
#' @importFrom dplyr filter mutate select distinct case_when pull
#' @importFrom tidyr pivot_wider
#' @importFrom UpSetR upset
#' @importFrom grDevices pdf dev.off
#' @importFrom grid grid.text gpar
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
#' # Genus persistence across Confidence Scores, kept in memory
#' res <- upset_kariocas(toy_project, tax_level = "Genus", export = FALSE)
#' res
#' head(res$data)
#' res$plots[["SAMPLE01_Bacteria_Genus"]]
#'
#' unlink(toy_project, recursive = TRUE)
upset_kariocas <- function(project_dir, tax_level = "Species", export = TRUE) {
    setup <- .ups_setup(project_dir, export)
    setup$log_msg(">>> Loading Data...")
    df_long <- .get_tidy_data(project_dir)
    if (!tax_level %in% unique(df_long$Rank)) {
        stop(
            "Rank '", tax_level, "' not found. Available: ",
            paste(sort(unique(df_long$Rank)), collapse = ", ")
        )
    }
    SAMPLES <- unique(df_long$sample)
    DOMAINS <- names(get_kariocas_colors("domains"))
    setup$log_msg(
        ">>> Starting UpSet Analysis (", tax_level, ") for ",
        length(SAMPLES), " samples."
    )
    results <- lapply(SAMPLES, function(samp) {
        .ups_process_sample(df_long, samp, DOMAINS, tax_level, setup)
    })
    .kcs_warn_failed_plots(unlist(lapply(results, `[[`, "failed")), setup$log_msg)
    plots <- do.call(c, unname(lapply(results, `[[`, "plots")))
    .kcs_finish_step(
        "005_taxa_intersections_across_CS",
        dplyr::bind_rows(lapply(results, `[[`, "data")),
        if (is.null(plots)) list() else plots,
        unlist(lapply(results, `[[`, "paths")),
        setup, setup$log_msg,
        what = "UpSet plots"
    )
}
