# ==============================================================================
# SHARED EXPORT / RESULT INFRASTRUCTURE
# ==============================================================================
# Every analysis step can either write its outputs to <project_dir>/<step>/
# (export = TRUE, the default) or keep everything in memory (export = FALSE).
# In both cases the step returns a `kariocas_result` object.
# ==============================================================================

#' Build the logging closure used by every step
#'
#' Messages always go to the console via \code{message()}. When
#' \code{log_file} is not \code{NULL} they are also appended to that file.
#' @param log_file Path of the log file, or \code{NULL} for console only.
#' @param header Character vector written once at the top of the log file.
#' @return A function \code{log_msg(...)}.
#' @noRd
.kcs_make_logger <- function(log_file = NULL, header = character(0)) {
    if (!is.null(log_file)) {
        writeLines(header, con = log_file)
    }
    function(...) {
        msg <- paste0(...)
        message(msg)
        if (!is.null(log_file)) {
            write(
                paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", msg),
                file = log_file, append = TRUE
            )
        }
        invisible(msg)
    }
}

#' Prepare the output folder and logger of a step
#'
#' @param project_dir Root path of the project (must exist).
#' @param step Sub-folder name of the step, e.g. \code{"002_taxa_retention"}.
#' @param log_name File name of the log, e.g. \code{"log_002_taxa_retention.txt"}.
#' @param export Logical. When \code{FALSE}, nothing is created on disk and the
#'   logger writes to the console only.
#' @param extra_header Optional extra lines for the log header.
#' @return A list with \code{output_dir} (\code{NULL} if not exporting),
#'   \code{log_file}, \code{log_msg} and \code{export}.
#' @noRd
.kcs_setup_step <- function(project_dir, step, log_name, export = TRUE,
                            extra_header = character(0)) {
    if (!dir.exists(project_dir)) {
        stop("Project directory not found: ", project_dir)
    }
    if (!isTRUE(export) && !isFALSE(export)) {
        stop("'export' must be TRUE or FALSE.")
    }
    output_dir <- NULL
    log_file <- NULL
    if (export) {
        output_dir <- file.path(project_dir, step)
        log_dir <- file.path(project_dir, "logs")
        if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
        if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)
        log_file <- file.path(log_dir, log_name)
    }
    header <- c(
        "====================================================",
        paste0("LOG: ", toupper(step)),
        paste0("PROJECT DIR: ", project_dir),
        extra_header,
        "===================================================="
    )
    list(
        output_dir = output_dir,
        log_file = log_file,
        log_msg = .kcs_make_logger(log_file, header),
        export = export
    )
}

#' Save a ggplot/patchwork object as PDF when exporting
#'
#' @param plot A ggplot or patchwork object.
#' @param fname File name (PDF).
#' @param setup The list returned by \code{.kcs_setup_step()}; \code{output_dir}
#'   may be overridden with \code{dir}.
#' @param dir Optional target directory (defaults to \code{setup$output_dir}).
#' @param width,height Page size in inches.
#' @return The full path of the written file, or \code{NULL} if not exporting.
#' @noRd
.kcs_save_plot <- function(plot, fname, setup, dir = NULL,
                           width = get_kariocas_dims()$width,
                           height = get_kariocas_dims()$height) {
    if (!isTRUE(setup$export)) {
        return(NULL)
    }
    dir <- if (is.null(dir)) setup$output_dir else dir
    if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
    path <- file.path(dir, fname)
    ggplot2::ggsave(path, plot, width = width, height = height)
    setup$log_msg("    -> Generated: ", fname)
    path
}

#' Save a data frame as TSV + RDS when exporting
#'
#' @param df Data frame.
#' @param stem File name without extension.
#' @param setup The list returned by \code{.kcs_setup_step()}.
#' @return Character vector of written paths (empty if not exporting).
#' @noRd
.kcs_save_table <- function(df, stem, setup) {
    if (!isTRUE(setup$export)) {
        return(character(0))
    }
    tsv_path <- file.path(setup$output_dir, paste0(stem, ".tsv"))
    rds_path <- file.path(setup$output_dir, paste0(stem, ".rds"))
    readr::write_tsv(df, tsv_path)
    readr::write_rds(df, rds_path)
    setup$log_msg("    -> Saved: ", basename(tsv_path), " / ", basename(rds_path))
    c(tsv_path, rds_path)
}

#' Finish a step: honest completion message and result object
#'
#' Emits \code{SUCCESS} only when at least one result (table row or plot) was
#' produced; otherwise raises a warning explaining that nothing was generated.
#' @inheritParams kariocas_result
#' @param log_msg Logging closure.
#' @return A \code{kariocas_result} object (visibly).
#' @noRd
.kcs_finish_step <- function(step, data, plots, paths, setup, log_msg,
                             what = "result(s)") {
    n_rows <- if (is.null(data)) 0L else nrow(data)
    n_plots <- length(plots)
    if (n_rows == 0L && n_plots == 0L) {
        log_msg("FAILED: no ", what, " were generated for ", step, ".")
        warning(
            "karioCaS ", step, ": no ", what, " generated. ",
            "Check the input data and the arguments used.",
            call. = FALSE
        )
    } else {
        log_msg(
            "SUCCESS: ", step, " completed (",
            n_rows, " table row(s), ", n_plots, " plot(s)",
            if (isTRUE(setup$export)) paste0(", ", length(paths), " file(s) written") else "",
            ")."
        )
    }
    kariocas_result(
        step = step, data = data, plots = plots, paths = paths,
        output_dir = setup$output_dir
    )
}

#' Result of a karioCaS analysis step
#'
#' Every analysis function of karioCaS returns an object of class
#' \code{kariocas_result}. It bundles the computed table(s), the generated
#' plots and, when \code{export = TRUE}, the paths of the files written to
#' the project folder. Elements are accessed with \code{$}.
#'
#' @param step Character. Name of the step, e.g. \code{"002_taxa_retention"}.
#' @param data A data frame with the computed results, or \code{NULL}.
#' @param plots A named list of \code{ggplot}/\code{patchwork} objects (may be
#'   empty).
#' @param paths Character vector of files written (empty if nothing exported).
#' @param output_dir Output folder, or \code{NULL} when not exporting.
#'
#' @return A list of class \code{kariocas_result} with elements \code{step},
#'   \code{data}, \code{plots}, \code{paths} and \code{output_dir}.
#' @examples
#' res <- kariocas_result(
#'     step = "demo", data = data.frame(x = 1:3),
#'     plots = list(), paths = character(0), output_dir = NULL
#' )
#' res
#' res$data
#' @export
kariocas_result <- function(step, data = NULL, plots = list(),
                            paths = character(0), output_dir = NULL) {
    structure(
        list(
            step = step, data = data, plots = plots,
            paths = paths, output_dir = output_dir
        ),
        class = "kariocas_result"
    )
}

#' @describeIn kariocas_result Compact summary of a result object.
#' @param x A \code{kariocas_result} object.
#' @param ... Ignored.
#' @export
print.kariocas_result <- function(x, ...) {
    n_rows <- if (is.null(x$data)) 0L else nrow(x$data)
    cat("<kariocas_result> step: ", x$step, "\n", sep = "")
    cat("  data  : ", n_rows, " row(s)",
        if (n_rows > 0) paste0(" x ", ncol(x$data), " column(s)") else "",
        "\n",
        sep = ""
    )
    cat("  plots : ", length(x$plots),
        if (length(x$plots) > 0) {
            paste0(
                " [", paste(utils::head(names(x$plots), 3), collapse = ", "),
                if (length(x$plots) > 3) ", ..." else "", "]"
            )
        } else {
            ""
        },
        "\n",
        sep = ""
    )
    if (is.null(x$output_dir)) {
        cat("  export: none (export = FALSE)\n")
    } else {
        cat("  export: ", length(x$paths), " file(s) in ", x$output_dir, "\n", sep = "")
    }
    invisible(x)
}

#' Draw an UpSetR object into a PDF file
#'
#' UpSetR draws directly on a graphics device, so it cannot be saved with
#' \code{ggsave()}. Drawing errors are caught and reported to the caller
#' (never swallowed): on failure the returned path carries an \code{"error"}
#' attribute and the partial file is removed.
#' @param up An object returned by \code{UpSetR::upset()}.
#' @param path Output PDF path.
#' @param title Title drawn at the top of the page.
#' @param log_msg Logging closure.
#' @param draw Function used to draw \code{up} (injectable for testing).
#' @return \code{path}, invisibly.
#' @noRd
.kcs_draw_upset_pdf <- function(up, path, title, log_msg,
                                draw = methods::show) {
    if (!dir.exists(dirname(path))) dir.create(dirname(path), recursive = TRUE)
    grDevices::pdf(
        file = path, width = get_kariocas_dims()$width,
        height = get_kariocas_dims()$height, onefile = FALSE
    )
    err <- tryCatch(
        {
            # UpSetR objects are drawn by their print method; show() dispatches
            # to it (BiocCheck discourages print() outside show methods).
            draw(up)
            grid::grid.text(
                label = title, x = 0.5, y = 0.98,
                gp = grid::gpar(fontsize = 14, fontface = "bold")
            )
            NULL
        },
        error = function(e) conditionMessage(e),
        finally = grDevices::dev.off()
    )
    if (!is.null(err)) {
        log_msg("    ERROR plotting ", basename(path), ": ", err)
        unlink(path)
        attr(path, "error") <- err
        return(invisible(path))
    }
    log_msg("    -> Generated: ", basename(path))
    invisible(path)
}

#' Turn failed UpSet drawings into a visible warning
#' @param failed Character vector of "plot: error message" entries.
#' @param log_msg Logging closure.
#' @noRd
.kcs_warn_failed_plots <- function(failed, log_msg) {
    if (length(failed) == 0) {
        return(invisible(NULL))
    }
    log_msg("  [WARNING] ", length(failed), " plot(s) could not be drawn.")
    warning(
        length(failed), " UpSet plot(s) could not be drawn and were not ",
        "written:\n  ", paste(failed, collapse = "\n  "),
        call. = FALSE
    )
    invisible(NULL)
}
