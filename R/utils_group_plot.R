# ==============================================================================
# SHARED GROUP-OVERLAY ENGINE
# ==============================================================================
# Reused by taxa_retention() (001) and reads_per_taxa() (003) to draw a single
# per-group figure in which every sample is a faint line and the group mean is a
# bold line with a +/-SD band, faceted 2x2 by Domain.
# ==============================================================================

#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_ribbon geom_vline
#'   labs scale_x_continuous scale_y_continuous ggsave theme element_text
#' @importFrom dplyr filter group_by summarise mutate arrange ungroup n_distinct
#' @importFrom tidyr replace_na
#' @importFrom stats sd setNames
#' @importFrom patchwork plot_annotation
NULL

#' Derive the biological group of a sample from its name
#'
#' The group is the sample name with any trailing run of digits removed, so
#' \code{SAMPLE33}, \code{SAMPLE34}, ... all collapse to \code{SAMPLE}, while
#' \code{CONTROL01} and \code{TREATED01} become \code{CONTROL} and \code{TREATED}.
#' Names without a trailing number are used verbatim.
#' @param sample Character vector of sample names.
#' @return Character vector of group labels.
#' @noRd
.grp_parse_group <- function(sample) {
    grp <- sub("[0-9]+$", "", sample)
    ifelse(grp == "", sample, grp)
}

#' Resolve the `detail_samples` argument to a vector of sample names
#'
#' Accepts \code{NULL} (no detail), the string \code{"all"}, a comma-separated
#' string such as \code{"SAMPLE33, SAMPLE45"}, or a character vector. Unknown
#' names are dropped with a warning.
#' @param detail_samples User input.
#' @param all_samples Character vector of valid sample names.
#' @param log_msg Logging closure.
#' @return Character vector of sample names to render in detail (possibly empty).
#' @noRd
.grp_resolve_detail <- function(detail_samples, all_samples, log_msg) {
    if (is.null(detail_samples)) {
        return(character(0))
    }
    if (length(detail_samples) == 1 && tolower(detail_samples) == "all") {
        return(all_samples)
    }
    if (length(detail_samples) == 1 && grepl(",", detail_samples)) {
        detail_samples <- trimws(strsplit(detail_samples, ",")[[1]])
    } else {
        detail_samples <- trimws(as.character(detail_samples))
    }
    detail_samples <- detail_samples[nzchar(detail_samples)]
    unknown <- setdiff(detail_samples, all_samples)
    if (length(unknown) > 0) {
        log_msg(
            "    [WARNING] Unknown sample(s) ignored: ",
            paste(unknown, collapse = ", ")
        )
    }
    intersect(detail_samples, all_samples)
}

#' Draw one domain panel of the group overlay
#'
#' @param df_dom Data frame for a single domain with columns \code{sample},
#'   \code{x}, \code{y}.
#' @param dom,dom_color Domain label and its colour.
#' @param x_lab,y_lab Axis titles (ggtext markdown allowed).
#' @param apply_scales Function taking a ggplot and returning it with scales.
#' @param vline Optional numeric x position for a dashed reference line (NA = none).
#' @noRd
.grp_overlay_domain <- function(df_dom, dom, dom_color, x_lab, y_lab,
                                apply_scales, vline = NA) {
    if (nrow(df_dom) == 0 || dplyr::n_distinct(df_dom$x) < 2) {
        return(plot_kariocas_empty(dom, x_label = x_lab, y_label = y_lab))
    }
    summ <- df_dom |>
        dplyr::group_by(.data$x) |>
        dplyr::summarise(
            Mean_Y = mean(.data$y, na.rm = TRUE),
            SD_Y = stats::sd(.data$y, na.rm = TRUE),
            .groups = "drop"
        ) |>
        dplyr::mutate(
            SD_Y = tidyr::replace_na(.data$SD_Y, 0),
            Lower = pmax(0, .data$Mean_Y - .data$SD_Y),
            Upper = .data$Mean_Y + .data$SD_Y
        )
    p <- ggplot2::ggplot()
    if (!is.na(vline)) {
        p <- p + ggplot2::geom_vline(
            xintercept = vline, linetype = "dashed",
            color = "grey40", linewidth = 0.5
        )
    }
    p <- p +
        ggplot2::geom_line(
            data = df_dom,
            ggplot2::aes(x = .data$x, y = .data$y, group = .data$sample),
            color = "grey70", linewidth = 0.4, alpha = 0.6
        ) +
        ggplot2::geom_ribbon(
            data = summ,
            ggplot2::aes(x = .data$x, ymin = .data$Lower, ymax = .data$Upper),
            fill = dom_color, alpha = 0.18
        ) +
        ggplot2::geom_line(
            data = summ,
            ggplot2::aes(x = .data$x, y = .data$Mean_Y),
            color = dom_color, linewidth = 1.2
        ) +
        ggplot2::geom_point(
            data = summ,
            ggplot2::aes(x = .data$x, y = .data$Mean_Y),
            color = dom_color, size = 2
        ) +
        ggplot2::labs(title = dom, x = x_lab, y = y_lab) +
        theme_kariocas()
    apply_scales(p)
}

#' Build the four domain panels of a group overlay
#'
#' @param df Data frame with columns \code{sample}, \code{Domain}, \code{x},
#'   \code{y}.
#' @param domains Character vector of domains (panel order by name).
#' @param x_lab,y_lab Axis titles.
#' @param apply_scales Scale-adding function.
#' @param vlines Optional named numeric vector of reference x per domain.
#' @noRd
.grp_overlay_plots <- function(df, domains, x_lab, y_lab,
                               apply_scales, vlines = NULL) {
    dom_colors <- get_kariocas_colors("domains")
    stats::setNames(
        lapply(domains, function(dom) {
            df_dom <- dplyr::filter(df, .data$Domain == dom)
            vl <- if (!is.null(vlines) && dom %in% names(vlines)) {
                vlines[[dom]]
            } else {
                NA
            }
            .grp_overlay_domain(
                df_dom, dom, dom_colors[[dom]],
                x_lab, y_lab, apply_scales, vline = vl
            )
        }),
        domains
    )
}

#' Assemble a 2x2 domain panel and save it
#'
#' @param plots Named list of ggplots (Bacteria/Archaea/Eukaryota/Viruses).
#' @param title,subtitle Annotation strings.
#' @param fname Output file name.
#' @param output_dir,log_msg Output dir and logging closure.
#' @noRd
.grp_assemble_2x2 <- function(plots, title, subtitle, fname, output_dir, log_msg) {
    layout <- (plots[["Bacteria"]] | plots[["Archaea"]]) /
        (plots[["Eukaryota"]] | plots[["Viruses"]]) +
        patchwork::plot_annotation(
            title = title,
            subtitle = subtitle,
            theme = ggplot2::theme(
                plot.title = ggplot2::element_text(
                    face = "bold", size = 16, hjust = 0.5
                ),
                plot.subtitle = ggplot2::element_text(
                    size = 12, hjust = 0.5, color = "grey30"
                )
            )
        )
    ggplot2::ggsave(
        file.path(output_dir, fname), layout,
        width = get_kariocas_dims()$width,
        height = get_kariocas_dims()$height
    )
    log_msg("    -> Generated: ", fname)
}
