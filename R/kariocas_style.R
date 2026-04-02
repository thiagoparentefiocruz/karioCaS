#' karioCaS Visualization Styles & Helpers
#'
#' Single Source of Truth for visual elements in karioCaS.
#' Follows Nature/Science publication standards.
#'
#' @export
#' @importFrom ggplot2 theme_classic theme element_text element_blank element_line element_rect margin unit guide_legend ggplot annotate labs scale_y_continuous scale_x_continuous
#' @importFrom scales label_number cut_short_scale log_trans
#' @importFrom ggtext element_markdown

# ==============================================================================
# 1. CORE VISUAL ASSETS (Hidden from direct user manipulation)
# ==============================================================================

.kariocas_internal_colors <- list(
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

.kariocas_internal_shapes <- list(
  ranks = c(
    "Phylum" = 15, "Class" = 16, "Order" = 17,
    "Family" = 18, "Genus" = 25, "Species" = 8,
    "Total Reads" = 15, "Level Reads" = 17, "Level Taxa" = 16
  )
)

.kariocas_internal_linetypes <- c(
  "Total Reads" = "longdash", "Level Reads" = "longdash", "Level Taxa" = "solid"
)

# ==============================================================================
# 2. ACCESSOR FUNCTIONS (The Bioconductor Way)
# ==============================================================================

#' Get karioCaS Official Colors
#'
#' Securely retrieves the official color palettes used across the karioCaS package.
#'
#' @param type Character. The palette to retrieve: "ranks", "domains", "special", "upset", or "heatmap".
#' @return A named vector or list of colors.
#' @export
get_kariocas_colors <- function(type = c("ranks", "domains", "special", "upset", "heatmap")) {
  type <- match.arg(type)
  return(.kariocas_internal_colors[[type]])
}

#' Get karioCaS Official Shapes
#'
#' Securely retrieves the official ggplot2 shapes for taxonomic ranks and special metrics.
#'
#' @param type Character. Currently supports "ranks".
#' @return A named vector of shape integers.
#' @export
get_kariocas_shapes <- function(type = "ranks") {
  if (type == "ranks") return(.kariocas_internal_shapes$ranks)
  stop("Shape type not recognized.")
}

#' Get karioCaS Official Linetypes
#'
#' Retrieves standard linetypes for read and taxa retention plots.
#'
#' @return A named character vector of linetypes.
#' @export
get_kariocas_linetypes <- function() {
  return(.kariocas_internal_linetypes)
}

# ==============================================================================
# 3. FORMATTERS & LABELS
# ==============================================================================

#' Format numbers with K/M suffixes
#' @export
label_k_number <- function(x) {
  scales::label_number(scale_cut = scales::cut_short_scale())(x)
}

#' Smart Number Formatter for Logs
#' Default uses comma for thousands to avoid conflict with decimal dot.
#' @export
label_kariocas_auto <- function(x) {
  format(x, scientific = FALSE, big.mark = ",", drop0trailing = TRUE)
}

#' Standard Axis Labels with HTML/Markdown Styling
#' Uses ggtext syntax for mixing fonts and sizes.
#' @export
get_kariocas_labels <- function() {
  list(
    y_log10_retained = "**% Retained**<br><span style='font-size:8pt;color:grey40'>(axis scaled to log<sub>10</sub>)</span>",
    x_log10_reads    = "**Reads**<br><span style='font-size:8pt;color:grey40'>(axis scaled to log<sub>10</sub>)</span>",
    y_confidence     = "**Confidence Score (%)**"
  )
}

# ==============================================================================
# 4. CUSTOM SCALES
# ==============================================================================

#' Log10 Scale with Smart Integer Formatting
#'
#' @param labels A formatting function or vector. Defaults to label_kariocas_auto.
#' @param ... Other arguments passed to scale_y_continuous
#' @export
scale_y_kariocas_log10 <- function(labels = label_kariocas_auto, ...) {
  ggplot2::scale_y_continuous(
    trans = "log10",
    labels = labels,
    breaks = c(0.01, 0.1, 1, 10, 100),
    ...
  )
}

#' @rdname scale_y_kariocas_log10
#' @export
scale_x_kariocas_log10 <- function(labels = label_kariocas_auto, ...) {
  ggplot2::scale_x_continuous(
    trans = "log10",
    labels = labels,
    ...
  )
}

# ==============================================================================
# 5. HELPER PLOTS
# ==============================================================================

#' Generate a Standard "No Data" Plot
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
# 6. DIMENSIONS & THEME
# ==============================================================================

#' Standard Output Dimensions (A4 Landscape)
#' @export
get_kariocas_dims <- function() {
  list(width = 11.69, height = 8.27)
}

#' karioCaS Publication Theme
#' @export
theme_kariocas <- function(base_size = 12, base_family = "sans") {
  ggplot2::theme_classic(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      # Typography
      plot.title = ggplot2::element_text(
        face = "bold", size = 14, hjust = 0.5, color = "black",
        margin = ggplot2::margin(b = 5)
      ),
      plot.subtitle = ggplot2::element_text(
        size = 10, hjust = 0.5, color = "grey30",
        margin = ggplot2::margin(b = 10)
      ),
      
      # Axis Titles (Compact Margins)
      axis.title.y = ggtext::element_markdown(
        hjust = 0.5, size = 11, color = "black", lineheight = 1.2,
        margin = ggplot2::margin(r = 5)
      ),
      axis.title.x = ggtext::element_markdown(
        hjust = 0.5, size = 11, color = "black", lineheight = 1.2,
        margin = ggplot2::margin(t = 5)
      ),
      
      axis.text = ggplot2::element_text(size = 10, color = "black"),
      
      # Legend
      legend.position = "bottom",
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 10),
      legend.margin = ggplot2::margin(t = -5),
      legend.box = "horizontal",
      
      # Grid & Lines
      panel.grid.major.y = ggplot2::element_line(color = "grey90", linetype = "dotted"),
      panel.grid.major.x = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(linewidth = 0.5, color = "black"),
      
      # Facets
      strip.background = ggplot2::element_rect(fill = "grey95", color = NA),
      strip.text = ggplot2::element_text(face = "bold", size = 10, color = "black"),
      
      # Plot Margins
      plot.margin = ggplot2::margin(10, 15, 10, 10)
    )
}