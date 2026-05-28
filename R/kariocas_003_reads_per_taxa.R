# ==============================================================================
# PRIVATE HELPERS — reads_per_taxa()
# ==============================================================================

#' @noRd
.rpt_setup <- function(project_dir, analysis_level,
                       x_max_bac, x_max_arc, x_max_euk, x_max_vir) {
  output_dir <- file.path(project_dir, "003_cutoffs")
  log_dir    <- file.path(project_dir, "logs")
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  if (!dir.exists(log_dir))    dir.create(log_dir,    recursive = TRUE)
  log_file <- file.path(log_dir, "log_003_cutoff_analysis.txt")
  writeLines(c(
    "====================================================",
    "LOG: 003_CUTOFF_ANALYSIS (Saturation + Rare Taxa)",
    paste0("PROJECT DIR: ",      project_dir),
    paste0("ANALYSIS LEVEL: ",   analysis_level),
    paste0("RARE TAXA MAX: B=",  x_max_bac, "|A=", x_max_arc,
           "|E=", x_max_euk, "|V=", x_max_vir),
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
  rare_limits <- list(
    Bacteria  = x_max_bac, Archaea   = x_max_arc,
    Eukaryota = x_max_euk, Viruses   = x_max_vir
  )
  list(output_dir = output_dir, log_msg = log_msg, rare_limits = rare_limits)
}

#' @noRd
.rpt_cutoff_template <- function() {
  base_steps  <- c(1, 3, 5)
  multipliers <- c(1, 10, 100, 1000, 10000, 100000, 1000000)
  sort(unique(c(0, as.vector(outer(base_steps, multipliers, "*")))))
}

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
.rpt_build_plot_data <- function(df_stats, analysis_level) {
  df_stats |>
    tidyr::pivot_longer(
      cols      = c("Ret_Reads_Pct", "Ret_Taxa_Pct"),
      names_to  = "Metric",
      values_to = "Pct"
    ) |>
    dplyr::mutate(
      Metric_Key = dplyr::case_when(
        .data$Metric == "Ret_Reads_Pct" ~ "Level Reads",
        .data$Metric == "Ret_Taxa_Pct"  ~ "Level Taxa"
      ),
      Metric_Key = factor(.data$Metric_Key,
                          levels = c("Level Taxa", "Level Reads"))
    )
}

#' @noRd
.rpt_base_plot <- function(df_plot, analysis_level, dom,
                           total_reads, total_taxa, legend_labels) {
  subtitle_stats <- paste0(
    label_kariocas_auto(total_reads), " Reads | ",
    label_kariocas_auto(total_taxa),  " ", analysis_level
  )
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
      breaks = seq(0, 1, 0.25),
      labels = function(x) x * 100
    ) +
    ggplot2::scale_color_manual(
      values = get_kariocas_colors("special"), labels = legend_labels
    ) +
    ggplot2::scale_shape_manual(
      values = get_kariocas_shapes("ranks"),   labels = legend_labels
    ) +
    ggplot2::scale_linetype_manual(
      values = get_kariocas_linetypes(),       labels = legend_labels
    ) +
    ggplot2::labs(
      title = dom, subtitle = subtitle_stats,
      y = "**% Retained**"
    ) +
    theme_kariocas() +
    ggplot2::guides(
      color    = ggplot2::guide_legend(nrow = 1, title = NULL),
      shape    = ggplot2::guide_legend(nrow = 1, title = NULL),
      linetype = ggplot2::guide_legend(nrow = 1, title = NULL)
    ) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5))
}

#' @noRd
.rpt_apply_saturation_scales <- function(p, max_val) {
  display_limit <- min(max_val, 10000)
  all_powers    <- c(1, 10, 100, 1000, 10000, 100000, 1000000)
  clean_breaks  <- all_powers[all_powers <= display_limit]
  p +
    scale_x_kariocas_log10(breaks = clean_breaks, labels = label_k_number) +
    ggplot2::coord_cartesian(
      xlim = c(NA, display_limit), ylim = c(0, 1.05), clip = "on"
    ) +
    ggplot2::labs(x = get_kariocas_labels()$x_log10_reads)
}

#' @noRd
.rpt_apply_raretaxa_scales <- function(p, dom_limit) {
  x_breaks <- if (dom_limit <= 20) {
    seq(1, dom_limit)
  } else {
    b <- unique(round(scales::breaks_pretty(n = 10)(c(1, dom_limit))))
    b[b >= 1 & b <= dom_limit]
  }
  p +
    ggplot2::scale_x_continuous(breaks = x_breaks) +
    ggplot2::coord_cartesian(
      xlim = c(1, dom_limit), ylim = c(0, 1.05), clip = "on"
    ) +
    ggplot2::labs(x = "**Reads**")
}

#' @noRd
.rpt_domain_plot <- function(df_curr, dom, mode, analysis_level,
                             rare_limits, cutoff_template_log) {
  df_dom      <- dplyr::filter(df_curr, .data$Domain == dom, .data$Counts > 0)
  total_reads <- sum(df_dom$Counts)
  total_taxa  <- nrow(df_dom)
  empty_x     <- if (mode == "Saturation") {
    get_kariocas_labels()$x_log10_reads
  } else {
    "**Reads**"
  }
  if (total_reads == 0) {
    return(plot_kariocas_empty(
      title_text    = dom,
      subtitle_text = "No reads detected",
      x_label       = empty_x,
      y_label       = "**% Retained**"
    ))
  }
  max_val <- max(df_dom$Counts)
  cutoffs <- if (mode == "Saturation") {
    ct <- cutoff_template_log[cutoff_template_log <= max_val]
    nxt <- cutoff_template_log[which(cutoff_template_log > max_val)[1]]
    if (!is.na(nxt)) c(ct, nxt) else ct
  } else {
    seq(1, rare_limits[[dom]], by = 1)
  }
  df_stats <- .rpt_calc_stats(df_dom, cutoffs, total_reads, total_taxa)
  df_plot  <- .rpt_build_plot_data(df_stats, analysis_level)
  p        <- .rpt_base_plot(df_plot, analysis_level, dom,
                             total_reads, total_taxa,
                             c(analysis_level, "Reads"))
  if (mode == "Saturation") {
    .rpt_apply_saturation_scales(p, max_val)
  } else {
    .rpt_apply_raretaxa_scales(p, rare_limits[[dom]])
  }
}

#' @noRd
.rpt_save_panel <- function(plots, DOMAINS, samp, cs,
                            mode, analysis_level, output_dir, log_msg) {
  layout <- (plots[["Bacteria"]] | plots[["Archaea"]]) /
    (plots[["Eukaryota"]] | plots[["Viruses"]]) +
    patchwork::plot_annotation(
      title    = paste(samp, "- CS", sprintf("%02d", cs), "| Analysis:", mode),
      subtitle = paste("Retention of Reads vs", analysis_level),
      theme    = ggplot2::theme(
        plot.title    = ggplot2::element_text(
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
    "_Cutoff_", analysis_level, "_", mode, ".pdf"
  )
  ggplot2::ggsave(
    file.path(output_dir, fname), layout,
    width  = get_kariocas_dims()$width,
    height = get_kariocas_dims()$height
  )
  log_msg("    -> Saved [", mode, "]: ", fname)
}

#' @noRd
.rpt_process_cs <- function(df_proc, samp, cs, DOMAINS, modes,
                            analysis_level, rare_limits,
                            cutoff_template_log, output_dir, log_msg) {
  df_curr <- dplyr::filter(df_proc, .data$sample == samp, .data$CS == cs)
  if (nrow(df_curr) == 0) {
    log_msg("    Skipping CS", cs, ": No data.")
    return(invisible(NULL))
  }
  for (mode in modes) {
    plots <- stats::setNames(
      lapply(DOMAINS, function(dom) {
        .rpt_domain_plot(df_curr, dom, mode, analysis_level,
                         rare_limits, cutoff_template_log)
      }),
      DOMAINS
    )
    .rpt_save_panel(plots, DOMAINS, samp, cs,
                    mode, analysis_level, output_dir, log_msg)
  }
}

# ==============================================================================
# EXPORTED FUNCTION
# ==============================================================================

#' Generate Read Cutoff Saturation Analysis (Step 003)
#'
#' Performs a saturation analysis by progressively filtering taxa with low read
#' counts. Generates two views: "Saturation" (log scale overview) and
#' "Rare_Taxa" (linear scale focus on low counts).
#'
#' @param project_dir Path to the project root.
#' @param analysis_level Taxonomic rank to analyze (default: \code{"Species"}).
#' @param x_max_bac Integer. Max X-axis for Bacteria Rare Taxa plot (default: 10).
#' @param x_max_arc Integer. Max X-axis for Archaea Rare Taxa plot (default: 10).
#' @param x_max_euk Integer. Max X-axis for Eukaryota Rare Taxa plot (default: 10).
#' @param x_max_vir Integer. Max X-axis for Viruses Rare Taxa plot (default: 10).
#'
#' @return Invisibly returns \code{NULL}. PDF plots are saved to
#'   \code{<project_dir>/003_cutoffs/}.
#' @export
#' @importFrom dplyr filter mutate select group_by summarise arrange pull
#'   n distinct case_when left_join bind_rows
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_line geom_point scale_color_manual
#'   scale_shape_manual scale_linetype_manual scale_y_continuous
#'   scale_x_continuous labs coord_cartesian guides guide_legend
#'   element_text ggsave theme
#' @importFrom patchwork plot_layout plot_annotation
#' @importFrom scales percent label_number breaks_pretty
#' @importFrom ggtext element_markdown
#' @examples
#' toy_project <- system.file("extdata", "your_project_name", package = "karioCaS")
#'
#' # Basic usage (defaults to Species, x_max = 10 for all domains)
#' # reads_per_taxa(project_dir = toy_project)
#'
#' # Genus level with custom Rare Taxa zoom
#' # reads_per_taxa(
#' #   project_dir    = toy_project,
#' #   analysis_level = "Genus",
#' #   x_max_bac      = 50,
#' #   x_max_vir      = 5
#' # )
reads_per_taxa <- function(project_dir,
                           analysis_level = "Species",
                           x_max_bac      = 10,
                           x_max_arc      = 10,
                           x_max_euk      = 10,
                           x_max_vir      = 10) {
  setup <- .rpt_setup(project_dir, analysis_level,
                      x_max_bac, x_max_arc, x_max_euk, x_max_vir)
  setup$log_msg(">>> Loading Data (Auto-detected format)...")
  df_long <- .get_tidy_data(project_dir)
  df_proc <- dplyr::filter(df_long, .data$Rank == analysis_level)
  if (nrow(df_proc) == 0) {
    setup$log_msg("CRITICAL ERROR: No data found for Rank: ", analysis_level)
    stop("No data for specified rank.")
  }
  SAMPLES            <- unique(df_proc$sample)
  DOMAINS            <- names(get_kariocas_colors("domains"))
  CS_LIST            <- unique(df_proc$CS)
  cutoff_template    <- .rpt_cutoff_template()
  modes              <- c("Saturation", "Rare_Taxa")
  setup$log_msg(">>> Starting Cutoff Analysis for ", length(SAMPLES), " samples.")
  for (samp in SAMPLES) {
    setup$log_msg("------------------------------------------------")
    setup$log_msg("  Processing Sample: ", samp)
    for (cs in CS_LIST) {
      .rpt_process_cs(
        df_proc, samp, cs, DOMAINS, modes,
        analysis_level, setup$rare_limits,
        cutoff_template, setup$output_dir, setup$log_msg
      )
    }
  }
  setup$log_msg("SUCCESS: Cutoff analysis completed.")
  invisible(NULL)
}