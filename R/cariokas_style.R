#' karioCaS Visualization Style Guide
#'
#' Centralized definitions for themes, palettes, and visual constants
#' to ensure "Nature-style" consistency across all plots.
#'
#' @keywords internal
#' @importFrom ggplot2 theme_classic theme element_text element_blank element_line unit margin
#' @noRd

# --- 1. Paletas de Cores (The Nature Look) ---
karioCaS_cols <- list(
  # Cores principais por Domínio
  domain = c(
    "Bacteria"  = "#3C5488", # Navy Blue
    "Archaea"   = "#E64B35", # Red/Orange
    "Viruses"   = "#00A087", # Teal
    "Eukaryota" = "#4DBBD5"  # Light Blue
  ),
  # Cores para Resolução Taxonômica (004)
  resolution = c(
    "Child"  = "#3C5488", # Matches Bacteria
    "Parent" = "#E64B35"  # Matches Archaea
  ),
  # Cores para Heatmap (005)
  heatmap = c("white", "#3C5488"),
  # Cor única para gráficos monocromáticos (ex: UpSet)
  main = "#3C5488"
)

# --- 2. Formas e Linhas (Consistência entre 001 e 003) ---
karioCaS_shapes <- c("Bacteria"=16, "Archaea"=17, "Viruses"=15, "Eukaryota"=18)
karioCaS_linetypes <- c("Bacteria"="solid", "Archaea"="dashed", "Viruses"="dotdash", "Eukaryota"="dotted")

# --- 3. Tema Central (karioCaS Theme) ---
#' Theme karioCaS
#'
#' A unified theme based on theme_classic with specific adjustments for
#' bottom legends, bold titles, and clean standardized text.
#'
#' @param base_size Base font size (default: 12)
#' @return A ggplot2 theme object
#' @export
theme_karioCaS <- function(base_size = 12) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      # Legenda
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 9, face = "bold"),
      legend.title = ggplot2::element_blank(), # Cleaner look

      # Títulos
      plot.title = ggplot2::element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 10, color = "grey30"),

      # Eixos
      axis.text = ggplot2::element_text(color = "black"),

      # Margens e Grids Específicos
      plot.margin = ggplot2::margin(t = 10, r = 15, b = 10, l = 10)
    )
}
