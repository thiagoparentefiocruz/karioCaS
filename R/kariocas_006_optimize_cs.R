# ==============================================================================
# PRIVATE HELPERS — optimize_CS()
# ==============================================================================

#' @noRd
.ocs_setup <- function(project_dir, tax_level, method) {
  output_dir <- file.path(project_dir, "006_optimize_CS")
  log_dir    <- file.path(project_dir, "logs")
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  if (!dir.exists(log_dir))    dir.create(log_dir,    recursive = TRUE)
  log_file <- file.path(log_dir, "log_006_optimize_cs.txt")
  writeLines(c(
    "====================================================",
    "LOG: 006_OPTIMIZE_CS (Stability Index - SI)",
    paste0("PROJECT DIR: ", project_dir),
    paste0("TAXA LEVEL:  ", tax_level),
    paste0("METHOD:      ", toupper(method)),
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
.ocs_load_data <- function(project_dir, tax_level, log_msg) {
  log_msg(">>> Loading Data (Auto-detected format)...")
  df_long <- .get_tidy_data(project_dir)
  df_proc <- dplyr::filter(df_long, .data$Rank == tax_level)
  if (nrow(df_proc) == 0) {
    log_msg("CRITICAL ERROR: No data found for Rank: ", tax_level)
    stop("No data for specified rank.")
  }
  df_proc
}

#' @noRd
.ocs_method_segmented <- function(calc_df, n_pts, dom, log_msg) {
  if (n_pts >= 5) {
    best_rss    <- Inf
    primary_idx <- 2
    for (i in 3:(n_pts - 2)) {
      fit_l     <- stats::lm(Pct_Retained ~ CS, data = calc_df[seq_len(i), ])
      fit_r     <- stats::lm(Pct_Retained ~ CS, data = calc_df[i:n_pts, ])
      rss_total <- sum(stats::resid(fit_l)^2) + sum(stats::resid(fit_r)^2)
      if (rss_total < best_rss) {
        best_rss    <- rss_total
        primary_idx <- i
      }
    }
    primary_cs <- calc_df$CS[primary_idx]
  } else {
    primary_cs <- calc_df$CS[min(2, n_pts)]
  }
  log_msg(sprintf("    %s -> Regime Shift: CS %02d", dom, primary_cs))
  list(
    primary_cs = primary_cs,
    sec1_cs    = NA,
    sub_txt    = sprintf(
      "Method: Segmented | Shift Breakpoint: CS %02d", primary_cs
    )
  )
}

#' @noRd
.ocs_method_dynamic <- function(calc_df, min_cs, max_cs, dom, log_msg) {
  tail_df <- dplyr::filter(calc_df, .data$CS >= 50)
  if (nrow(tail_df) >= 2) {
    bg_mean      <- mean(tail_df$Step_Loss_Pct, na.rm = TRUE)
    bg_sd        <- stats::sd(tail_df$Step_Loss_Pct, na.rm = TRUE)
    dynamic_toll <- max(bg_mean + (1.5 * bg_sd), 0.5)
    ultra_toll   <- max(bg_mean, 0.1)
  } else {
    dynamic_toll <- 1.0
    ultra_toll   <- 0.2
  }
  stable    <- dplyr::filter(
    calc_df, .data$CS > min_cs, .data$Step_Loss_Pct <= dynamic_toll
  )
  primary_cs <- if (nrow(stable) > 0) stable$CS[1] else max_cs
  ultra      <- dplyr::filter(
    calc_df, .data$CS > primary_cs, .data$Step_Loss_Pct <= ultra_toll
  )
  sec1_cs <- if (nrow(ultra) >= 1) ultra$CS[1] else NA
  log_msg(sprintf(
    "    %s -> Dynamic Toll: %.2f%% | Primary: CS %02d", dom, dynamic_toll, primary_cs
  ))
  list(
    primary_cs = primary_cs,
    sec1_cs    = sec1_cs,
    sub_txt    = sprintf(
      "Method: Dynamic (Toll: %.2f%%) | Primary: CS %02d | Sec: %s",
      dynamic_toll, primary_cs,
      ifelse(is.na(sec1_cs), "NA", sprintf("%02d", sec1_cs))
    )
  )
}

#' @noRd
.ocs_method_manual <- function(calc_df, min_cs, max_cs,
                               dom, manual_toll, log_msg) {
  current_toll <- if (is.list(manual_toll) && !is.null(manual_toll[[dom]])) {
    manual_toll[[dom]]
  } else if (is.numeric(manual_toll) && length(manual_toll) == 1) {
    manual_toll
  } else {
    1.0
  }
  stable     <- dplyr::filter(
    calc_df, .data$CS > min_cs, .data$Step_Loss_Pct <= current_toll
  )
  primary_cs <- if (nrow(stable) > 0) stable$CS[1] else max_cs
  ultra      <- dplyr::filter(
    calc_df, .data$CS > primary_cs, .data$Step_Loss_Pct <= 0.2
  )
  sec1_cs <- if (nrow(ultra) >= 1) ultra$CS[1] else NA
  log_msg(sprintf(
    "    %s -> Manual Toll: %.2f%% | Primary: CS %02d", dom, current_toll, primary_cs
  ))
  list(
    primary_cs = primary_cs,
    sec1_cs    = sec1_cs,
    sub_txt    = sprintf(
      "Method: Manual (Toll: %.2f%%) | Primary: CS %02d | Sec: %s",
      current_toll, primary_cs,
      ifelse(is.na(sec1_cs), "NA", sprintf("%02d", sec1_cs))
    )
  )
}

#' @noRd
.ocs_tag_audit <- function(calc_df, dom, samp, primary_cs, sec1_cs) {
  calc_df |>
    dplyr::mutate(
      Domain  = dom,
      Sample  = samp,
      SI_Type = dplyr::case_when(
        .data$CS == primary_cs                    ~ "Primary_SI",
        !is.na(sec1_cs) & .data$CS == sec1_cs    ~ "Secondary_SI_1",
        TRUE                                       ~ NA_character_
      )
    ) |>
    dplyr::select(
      "Sample", "Domain", "CS", "Taxa_Count",
      "Pct_Retained", "Step_Loss_Pct", "SI_Type"
    )
}

#' @noRd
.ocs_domain_plot <- function(calc_df, dom, sub_txt) {
  spec_colors <- get_kariocas_colors("special")
  labels_vec  <- get_kariocas_labels()
  si_pts      <- dplyr::filter(calc_df, !is.na(.data$SI_Type))
  ggplot2::ggplot(calc_df, ggplot2::aes(x = .data$CS, y = .data$Pct_Retained)) +
    ggplot2::geom_line(color = "black", linewidth = 1) +
    ggplot2::geom_point(
      color = "black", size = 2.5,
      shape = get_kariocas_shapes("ranks")[["Level Taxa"]]
    ) +
    ggplot2::geom_point(
      data  = si_pts,
      ggplot2::aes(x = .data$CS, y = .data$Pct_Retained),
      color = spec_colors[["Level Taxa"]], size = 4, shape = 4, stroke = 1.5
    ) +
    ggplot2::scale_x_continuous(breaks = seq(0, 100, 20), limits = c(0, 100)) +
    ggplot2::scale_y_continuous(limits = c(0, 105)) +
    ggplot2::labs(
      title    = dom,
      subtitle = sub_txt,
      x        = labels_vec$y_confidence,
      y        = "**% Retained**"
    ) +
    theme_kariocas()
}

#' @noRd
.ocs_process_domain <- function(df_samp, dom, method,
                                manual_toll, samp, log_msg) {
  df_dom   <- dplyr::filter(df_samp, .data$Domain == dom)
  stats_df <- df_dom |>
    dplyr::group_by(.data$CS) |>
    dplyr::summarise(
      Taxa_Count = dplyr::n_distinct(.data$Taxon_Name), .groups = "drop"
    ) |>
    dplyr::arrange(.data$CS)
  if (nrow(stats_df) < 3) {
    log_msg("    Skipping ", dom, ": Not enough CS data points for curve fitting.")
    return(list(plot = plot_kariocas_empty(dom, "Insufficient CS points"),
                audit = NULL))
  }
  min_cs   <- min(stats_df$CS)
  max_cs   <- max(stats_df$CS)
  max_taxa <- max(stats_df$Taxa_Count)
  min_taxa <- min(stats_df$Taxa_Count)
  if (max_taxa == min_taxa) {
    log_msg("    Skipping ", dom, ": Flat curve (no taxa loss).")
    return(list(plot = plot_kariocas_empty(dom, "Flat Curve (No elbow)"),
                audit = NULL))
  }
  calc_df <- stats_df |>
    dplyr::arrange(.data$CS) |>
    dplyr::mutate(
      Pct_Retained  = (.data$Taxa_Count / max_taxa) * 100,
      Step_Loss_Pct = dplyr::lag(.data$Pct_Retained, default = 100) -
        .data$Pct_Retained
    )
  n_pts  <- nrow(calc_df)
  result <- switch(method,
                   segmented = .ocs_method_segmented(calc_df, n_pts, dom, log_msg),
                   dynamic   = .ocs_method_dynamic(calc_df, min_cs, max_cs, dom, log_msg),
                   manual    = .ocs_method_manual(calc_df, min_cs, max_cs,
                                                  dom, manual_toll, log_msg)
  )
  audit_df <- .ocs_tag_audit(
    calc_df, dom, samp, result$primary_cs, result$sec1_cs
  )
  list(
    plot  = .ocs_domain_plot(audit_df, dom, result$sub_txt),
    audit = audit_df
  )
}

#' @noRd
.ocs_save_panel <- function(plots, samp, tax_level, output_dir, log_msg) {
  layout <- (plots[["Bacteria"]] | plots[["Archaea"]]) /
    (plots[["Eukaryota"]] | plots[["Viruses"]]) +
    patchwork::plot_annotation(
      title    = paste(samp, "- CS Optimization (SI)"),
      subtitle = paste(
        "Level:", tax_level,
        "| Red 'X': Stability Indices (Primary & Secondary)"
      ),
      theme = ggplot2::theme(
        plot.title    = ggplot2::element_text(
          face = "bold", size = 16, hjust = 0.5
        ),
        plot.subtitle = ggplot2::element_text(
          size = 12, hjust = 0.5, color = "grey30"
        )
      )
    )
  fname <- paste0(samp, "_Optimize_CS_", tax_level, ".pdf")
  dims  <- get_kariocas_dims()
  ggplot2::ggsave(
    file.path(output_dir, fname), layout,
    width = dims$width, height = dims$height
  )
  log_msg("    -> Generated Plot: ", fname)
}

#' @noRd
.ocs_export_audit <- function(audit_list, tax_level, output_dir, log_msg) {
  log_msg(">>> Consolidating Audit Data...")
  if (length(audit_list) == 0) {
    log_msg("WARNING: No audit data generated. Check inputs.")
    return(invisible(NULL))
  }
  full_audit <- dplyr::bind_rows(audit_list)
  tsv_path   <- file.path(output_dir, paste0("SI_Audit_", tax_level, ".tsv"))
  rds_path   <- file.path(output_dir, paste0("SI_Audit_", tax_level, ".rds"))
  readr::write_tsv(full_audit, tsv_path)
  readr::write_rds(full_audit, rds_path)
  log_msg("SAVED AUDIT TSV: ", tsv_path)
  log_msg("SAVED AUDIT RDS: ", rds_path)
  full_audit
}

# ==============================================================================
# EXPORTED FUNCTION
# ==============================================================================

#' Optimize Confidence Score using Multi-Strategy Stability Index (Step 006)
#'
#' Objectively identifies the optimal Confidence Score (CS) threshold to separate
#' statistical noise from true biological signal using a multi-strategy engine:
#' \itemize{
#'   \item \code{"dynamic"}: Adapts to baseline noise (ideal for pathogens).
#'   \item \code{"segmented"}: Broken-stick regression for regime shifts
#'     (ideal for ecology/dark matter).
#'   \item \code{"manual"}: Expert-defined acceptable loss tolls.
#' }
#'
#' @param project_dir Character. Path to the project root.
#' @param tax_level Character. Taxonomic rank to analyze (default: \code{"Species"}).
#' @param method Character. One of \code{"dynamic"} (default), \code{"segmented"},
#'   or \code{"manual"}.
#' @param manual_toll Numeric or named list. Acceptable step-wise loss percentage.
#'   Used only when \code{method = "manual"} (default: 1.0).
#'
#' @return Invisibly returns a \code{data.frame} with the full SI audit trail.
#'   TSV and RDS copies are saved to \code{<project_dir>/006_optimize_CS/}.
#' @export
#' @importFrom dplyr filter select group_by summarise arrange mutate case_when
#'   bind_rows lag n_distinct
#' @importFrom ggplot2 ggplot aes geom_line geom_point scale_x_continuous
#'   scale_y_continuous labs theme element_text ggsave
#' @importFrom patchwork plot_layout plot_annotation
#' @importFrom readr write_tsv write_rds
#' @importFrom stats lm resid sd
#' @examples
#' toy_project <- system.file("extdata", "your_project_name", package = "karioCaS")
#'
#' # Dynamic method (default — ideal for clinical/pathogen focus)
#' # optimize_CS(project_dir = toy_project, tax_level = "Species")
#'
#' # Segmented regression (ideal for ecology / biological dark matter)
#' # optimize_CS(
#' #   project_dir = toy_project,
#' #   tax_level   = "Species",
#' #   method      = "segmented"
#' # )
optimize_CS <- function(project_dir,
                        tax_level    = "Species",
                        method       = c("dynamic", "segmented", "manual"),
                        manual_toll  = 1.0) {
  method  <- match.arg(method)
  setup   <- .ocs_setup(project_dir, tax_level, method)
  df_proc <- .ocs_load_data(project_dir, tax_level, setup$log_msg)
  SAMPLES <- unique(df_proc$sample)
  DOMAINS <- names(get_kariocas_colors("domains"))
  audit_list <- list()
  setup$log_msg(">>> Starting SI Calculation for ", length(SAMPLES), " samples.")
  for (samp in SAMPLES) {
    setup$log_msg("------------------------------------------------")
    setup$log_msg("  Processing Sample: ", samp)
    df_samp  <- dplyr::filter(df_proc, .data$sample == samp)
    results  <- lapply(
      stats::setNames(DOMAINS, DOMAINS),
      function(dom) .ocs_process_domain(
        df_samp, dom, method, manual_toll, samp, setup$log_msg
      )
    )
    plots      <- lapply(results, `[[`, "plot")
    audit_list <- c(audit_list, Filter(Negate(is.null), lapply(results, `[[`, "audit")))
    if (length(plots) > 0) {
      .ocs_save_panel(plots, samp, tax_level, setup$output_dir, setup$log_msg)
    }
  }
  full_audit <- .ocs_export_audit(
    audit_list, tax_level, setup$output_dir, setup$log_msg
  )
  setup$log_msg("SUCCESS: CS Optimization completed.")
  invisible(full_audit)
}