#' karioCaS Visualization Styles
#'
#' Single Source of Truth for visual elements in karioCaS.
#' Follows Nature/Science publication standards.
#'
#' @export
#' @importFrom ggplot2 theme_classic theme element_text element_blank element_line element_rect margin unit guide_legend ggplot annotate labs scale_y_continuous scale_x_continuous expansion
#' @importFrom scales label_number cut_short_scale log_trans
kariocas_colors <- list(
  ranks = c(
    "Phylum"  = "#000000", "Class"   = "#E69F00", "Order"   = "#56B4E9",
    "Family"  = "#009E73", "Genus"   = "#F0E442", "Species" = "#D55E00"
  ),
  domains = c(
    "Bacteria"  = "#3C5488", "Archaea"   = "#E64B35",
    "Viruses"   = "#00A087", "Eukaryota" = "#4DBBD5"
  ),
  special = c(
    "Total Reads" = "#000000", "Level Reads" = "#3C5488",
    "Level Taxa"  = "#E64B35", "Parent"      = "#8491B4", "Child" = "#E64B35"
  ),
  upset = list(main = "#3C5488", sets = "#7E6148"),
  heatmap = c("#4DBBD5", "#FFFFFF", "#E64B35")
)

#' Official karioCaS Shapes
#' @export
kariocas_shapes <- c(
  "Phylum" = 15, "Class" = 16, "Order" = 17,
  "Family" = 18, "Genus" = 25, "Species" = 8
)

#' Official karioCaS Line Types
#' @export
kariocas_linetypes <- c(
  "Total Reads" = "longdash", "Level Reads" = "longdash", "Level Taxa" = "solid"
)

# ==============================================================================
# 2. FORMATTERS & LABELS
# ==============================================================================

#' Format numbers with K/M suffixes (e.g., 10000 -> 10K)
#' @export
label_k_number <- function(x) {
  scales::label_number(scale_cut = scales::cut_short_scale())(x)
}

#' Smart Number Formatter for Logs
#' Displays integers as integers (100) and small fractions with decimals (0.01).
#' @export
label_kariocas_auto <- function(x) {
  format(x, scientific = FALSE, big.mark = ".", drop0trailing = TRUE)
}

#' Standard Axis Labels with Mathematical Notation
#' Uses atop() for line breaks in expressions.
#' @export
kariocas_labels <- list(
  y_log10_retained = expression(atop("% Retained", paste("(axis scaled to ", log[10], ")"))),
  x_log10_reads    = expression(atop("Reads", paste("(axis scaled to ", log[10], ")"))),
  y_confidence     = "Confidence Score (%)"
)

# ==============================================================================
# 3. CUSTOM SCALES
# ==============================================================================

#' Log10 Scale with Smart Integer Formatting
#' @param ... Arguments passed to scale_y_continuous
#' @export
scale_y_kariocas_log10 <- function(...) {
  ggplot2::scale_y_continuous(
    trans = "log10",
    labels = label_kariocas_auto,
    breaks = c(0.01, 0.1, 1, 10, 100),
    ...
  )
}

#' @rdname scale_y_kariocas_log10
#' @export
scale_x_kariocas_log10 <- function(...) {
  ggplot2::scale_x_continuous(
    trans = "log10",
    labels = label_kariocas_auto,
    ...
  )
}

# ==============================================================================
# 4. HELPER PLOTS
# ==============================================================================

#' Generate a Standard "No Data" Plot (Skeleton Style)
#' @export
plot_kariocas_empty <- function(title_text = "No Data",
                                subtitle_text = NULL,
                                x_label = NULL,
                                y_label = NULL) {
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No Data",
                      size = 5, fontface = "italic", color = "grey60") +
    ggplot2::labs(
      title = title_text,
      subtitle = subtitle_text,
      x = x_label,
      y = y_label
    ) +
    theme_kariocas() +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank()
    )
}

# ==============================================================================
# 5. DIMENSIONS & THEME
# ==============================================================================

#' Standard Output Dimensions (A4 Landscape)
#' @export
kariocas_dims <- list(width = 11.69, height = 8.27)

#' karioCaS Publication Theme
#' @export
theme_kariocas <- function(base_size = 12, base_family = "sans") {
  ggplot2::theme_classic(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      # Typography
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 4, hjust = 0.5, color = "black", margin = ggplot2::margin(b = 10)),
      plot.subtitle = ggplot2::element_text(size = base_size - 2, hjust = 0.5, color = "grey30", margin = ggplot2::margin(b = 15)),

      # Axis
      axis.title = ggplot2::element_text(face = "bold", size = base_size - 1, color = "black"),
      axis.text = ggplot2::element_text(size = base_size - 2, color = "black"),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 10)),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 10)),

      # Legend
      legend.position = "bottom",
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = base_size - 2),
      legend.margin = ggplot2::margin(t = -5),
      legend.box = "horizontal",

      # Grid & Lines
      panel.grid.major.y = ggplot2::element_line(color = "grey90", linetype = "dotted"),
      panel.grid.major.x = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(linewidth = 0.5, color = "black"),

      # Facets
      strip.background = ggplot2::element_rect(fill = "grey95", color = NA),
      strip.text = ggplot2::element_text(face = "bold", size = base_size - 2, color = "black"),

      # Margins
      plot.margin = ggplot2::margin(10, 15, 10, 10)
    )
}
