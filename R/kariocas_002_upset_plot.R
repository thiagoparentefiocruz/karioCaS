# ==============================================================================
# PRIVATE HELPERS - upset_kariocas()
# ==============================================================================

#' @noRd
.ups_setup <- function(project_dir) {
    output_dir <- file.path(project_dir, "002_UpSetComparison_Plots")
    log_dir <- file.path(project_dir, "logs")
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)
    log_file <- file.path(log_dir, "log_002_upsetplots.txt")
    writeLines(c(
        "====================================================",
        "LOG: 002_UPSET_PLOTS (karioCaS never are upset)",
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
    list(output_dir = output_dir, log_msg = log_msg)
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

#' @noRd
.ups_plot_one <- function(binary_matrix, upset_cols, samp, dom,
                          lvl, full_path, log_msg) {
    grDevices::pdf(
        file    = full_path,
        width   = get_kariocas_dims()$width,
        height  = get_kariocas_dims()$height,
        onefile = FALSE
    )
    tryCatch(
        {
            # UpSetR requires explicit method dispatch inside pdf() device.
            # show() is used here as the Bioconductor-preferred dispatcher.
            show(
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
            )
            grid::grid.text(
                label = paste(samp, "-", dom, "|", lvl, "Intersection Analysis"),
                x = 0.5, y = 0.98,
                gp = grid::gpar(fontsize = 16, fontface = "bold")
            )
            log_msg("    -> Generated: ", basename(full_path))
        },
        error = function(e) {
            log_msg("    ERROR plotting ", basename(full_path), ": ", e$message)
        }
    )
    grDevices::dev.off()
}

#' @noRd
.ups_process_sample <- function(df_long, samp, DOMAINS, LEVELS,
                                output_dir, log_msg) {
    log_msg("------------------------------------------------")
    log_msg("  Processing Sample: ", samp)
    samp_dir <- file.path(output_dir, samp)
    if (!dir.exists(samp_dir)) dir.create(samp_dir)
    for (dom in DOMAINS) {
        for (lvl in LEVELS) {
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
            fname <- paste0(samp, "_", dom, "_", lvl, "_UpSet.pdf")
            full_path <- file.path(samp_dir, fname)
            .ups_plot_one(mat, upset_cols, samp, dom, lvl, full_path, log_msg)
        }
    }
}

# ==============================================================================
# EXPORTED FUNCTION
# ==============================================================================

#' Generate UpSet Plots per Sample/Domain/Level (Step 002)
#'
#' "karioCaS never are upset!" Generates UpSet plots showing taxon persistence
#' across Confidence Score levels, with detailed logging.
#'
#' @param project_dir Path to the project root.
#'
#' @return Invisibly returns \code{NULL}. PDF plots are saved to
#'   \code{<project_dir>/002_UpSetComparison_Plots/}.
#' @export
#' @importFrom dplyr filter mutate select distinct case_when pull
#' @importFrom tidyr pivot_wider
#' @importFrom UpSetR upset
#' @importFrom grDevices pdf dev.off
#' @importFrom grid grid.text gpar
#' @examples
#' toy_project <- system.file("extdata", "your_project_name", package = "karioCaS")
#'
#' # upset_kariocas(project_dir = toy_project)
upset_kariocas <- function(project_dir) {
    if (!requireNamespace("UpSetR", quietly = TRUE)) {
        stop("Package 'UpSetR' is required.")
    }
    setup <- .ups_setup(project_dir)
    setup$log_msg(">>> Loading Data...")
    df_long <- .get_tidy_data(project_dir)
    SAMPLES <- unique(df_long$sample)
    DOMAINS <- names(get_kariocas_colors("domains"))
    LEVELS <- c("Species", "Genus", "Family")
    setup$log_msg(">>> Starting UpSet Analysis for ", length(SAMPLES), " samples.")
    for (samp in SAMPLES) {
        .ups_process_sample(
            df_long, samp, DOMAINS, LEVELS,
            setup$output_dir, setup$log_msg
        )
    }
    setup$log_msg("SUCCESS: UpSet analysis completed.")
    invisible(NULL)
}
