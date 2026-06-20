# ==============================================================================
# STABILITY INDEX (SI) ENGINE
# ==============================================================================
# Computes the optimal Confidence Score per sample/domain from a taxa-retention
# decay curve, and writes the machine-readable audit consumed by
# retrieve_selected_taxa(). Used by taxa_retention() (Step 001).
# ==============================================================================

#' @importFrom dplyr filter group_by summarise arrange mutate lag n_distinct
#'   case_when select bind_rows
#' @importFrom stats lm resid sd
#' @importFrom readr write_tsv write_rds
NULL

# ------------------------------------------------------------------------------
# Method strategies (each returns primary_cs, sec1_cs, sub_txt)
# ------------------------------------------------------------------------------

#' Tail-derived step-loss tolerance (shared by dynamic & post-cliff)
#' @noRd
.ocs_dyn_toll <- function(calc_df) {
    t <- calc_df$Step_Loss_Pct[calc_df$CS >= 50]
    t <- t[!is.na(t)]
    if (length(t) >= 2) {
        max(mean(t) + 1.5 * stats::sd(t), 0.5)
    } else {
        1.0
    }
}

#' Kneedle elbow: CS of maximum distance below the start-end chord
#' @noRd
.ocs_knee_cs <- function(calc_df) {
    x <- calc_df$CS
    y <- calc_df$Pct_Retained
    n <- length(x)
    if (n < 3 || x[n] == x[1]) {
        return(x[min(2, n)])
    }
    line <- y[1] + (y[n] - y[1]) * (x - x[1]) / (x[n] - x[1])
    d <- line - y
    if (all(is.na(d)) || max(d, na.rm = TRUE) <= 0) {
        # Non-convex curve: fall back to the steepest single drop.
        return(x[which.max(calc_df$Step_Loss_Pct)])
    }
    x[which.max(d)]
}

#' Post-cliff stabilisation: first CS at/after the biggest drop within tolerance
#' @noRd
.ocs_postcliff_cs <- function(calc_df, toll, max_cs) {
    pk <- which.max(calc_df$Step_Loss_Pct)
    idx <- which(
        seq_len(nrow(calc_df)) >= pk & calc_df$Step_Loss_Pct <= toll
    )
    if (length(idx) > 0) calc_df$CS[idx[1]] else max_cs
}

#' @noRd
.ocs_method_kneedle <- function(calc_df, max_cs, dom, log_msg) {
    primary_cs <- .ocs_knee_cs(calc_df)
    # Secondary SI: the more conservative post-cliff floor, if it is stricter.
    floor_cs <- .ocs_postcliff_cs(calc_df, .ocs_dyn_toll(calc_df), max_cs)
    sec1_cs <- if (!is.na(floor_cs) && floor_cs > primary_cs) floor_cs else NA
    log_msg(sprintf("    %s -> Kneedle Elbow: CS %02d", dom, primary_cs))
    list(
        primary_cs = primary_cs,
        sec1_cs = sec1_cs,
        sub_txt = sprintf(
            "Method: Kneedle | Elbow (Primary): CS %02d | Floor (Sec): %s",
            primary_cs,
            ifelse(is.na(sec1_cs), "NA", sprintf("%02d", sec1_cs))
        )
    )
}

#' @noRd
.ocs_method_postcliff <- function(calc_df, max_cs, dom, log_msg) {
    toll <- .ocs_dyn_toll(calc_df)
    primary_cs <- .ocs_postcliff_cs(calc_df, toll, max_cs)
    log_msg(sprintf(
        "    %s -> Post-cliff (Toll: %.2f%%): CS %02d", dom, toll, primary_cs
    ))
    list(
        primary_cs = primary_cs,
        sec1_cs = NA,
        sub_txt = sprintf(
            "Method: Post-cliff (Toll: %.2f%%) | Primary: CS %02d",
            toll, primary_cs
        )
    )
}

#' @noRd
.ocs_method_segmented <- function(calc_df, n_pts, dom, log_msg) {
    if (n_pts >= 5) {
        best_rss <- Inf
        primary_idx <- 2
        for (i in 3:(n_pts - 2)) {
            fit_l <- stats::lm(Pct_Retained ~ CS, data = calc_df[seq_len(i), ])
            fit_r <- stats::lm(Pct_Retained ~ CS, data = calc_df[i:n_pts, ])
            rss_total <- sum(stats::resid(fit_l)^2) + sum(stats::resid(fit_r)^2)
            if (rss_total < best_rss) {
                best_rss <- rss_total
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
        sec1_cs = NA,
        sub_txt = sprintf(
            "Method: Segmented | Shift Breakpoint: CS %02d", primary_cs
        )
    )
}

#' @noRd
.ocs_method_dynamic <- function(calc_df, min_cs, max_cs, dom, log_msg) {
    tail_df <- dplyr::filter(calc_df, .data$CS >= 50)
    if (nrow(tail_df) >= 2) {
        bg_mean <- mean(tail_df$Step_Loss_Pct, na.rm = TRUE)
        bg_sd <- stats::sd(tail_df$Step_Loss_Pct, na.rm = TRUE)
        dynamic_toll <- max(bg_mean + (1.5 * bg_sd), 0.5)
        ultra_toll <- max(bg_mean, 0.1)
    } else {
        dynamic_toll <- 1.0
        ultra_toll <- 0.2
    }
    stable <- dplyr::filter(
        calc_df, .data$CS > min_cs, .data$Step_Loss_Pct <= dynamic_toll
    )
    primary_cs <- if (nrow(stable) > 0) stable$CS[1] else max_cs
    ultra <- dplyr::filter(
        calc_df, .data$CS > primary_cs, .data$Step_Loss_Pct <= ultra_toll
    )
    sec1_cs <- if (nrow(ultra) >= 1) ultra$CS[1] else NA
    log_msg(sprintf(
        "    %s -> Dynamic Toll: %.2f%% | Primary: CS %02d",
        dom, dynamic_toll, primary_cs
    ))
    list(
        primary_cs = primary_cs,
        sec1_cs = sec1_cs,
        sub_txt = sprintf(
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
    stable <- dplyr::filter(
        calc_df, .data$CS > min_cs, .data$Step_Loss_Pct <= current_toll
    )
    primary_cs <- if (nrow(stable) > 0) stable$CS[1] else max_cs
    ultra <- dplyr::filter(
        calc_df, .data$CS > primary_cs, .data$Step_Loss_Pct <= 0.2
    )
    sec1_cs <- if (nrow(ultra) >= 1) ultra$CS[1] else NA
    log_msg(sprintf(
        "    %s -> Manual Toll: %.2f%% | Primary: CS %02d",
        dom, current_toll, primary_cs
    ))
    list(
        primary_cs = primary_cs,
        sec1_cs = sec1_cs,
        sub_txt = sprintf(
            "Method: Manual (Toll: %.2f%%) | Primary: CS %02d | Sec: %s",
            current_toll, primary_cs,
            ifelse(is.na(sec1_cs), "NA", sprintf("%02d", sec1_cs))
        )
    )
}

# ------------------------------------------------------------------------------
# Audit construction
# ------------------------------------------------------------------------------

#' @noRd
.ocs_tag_audit <- function(calc_df, dom, samp, primary_cs, sec1_cs) {
    calc_df |>
        dplyr::mutate(
            Domain = dom,
            Sample = samp,
            SI_Type = dplyr::case_when(
                .data$CS == primary_cs ~ "Primary_SI",
                !is.na(sec1_cs) & .data$CS == sec1_cs ~ "Secondary_SI_1",
                TRUE ~ NA_character_
            )
        ) |>
        dplyr::select(
            "Sample", "Domain", "CS", "Taxa_Count",
            "Pct_Retained", "Step_Loss_Pct", "SI_Type"
        )
}

#' Compute the SI audit for one domain of one sample (no plotting).
#' @return A tagged audit data frame, or NULL if the curve is unusable.
#' @noRd
.si_domain_audit <- function(df_samp, dom, method, manual_toll, samp, log_msg) {
    df_dom <- dplyr::filter(df_samp, .data$Domain == dom)
    stats_df <- df_dom |>
        dplyr::group_by(.data$CS) |>
        dplyr::summarise(
            Taxa_Count = dplyr::n_distinct(.data$Taxon_Name), .groups = "drop"
        ) |>
        dplyr::arrange(.data$CS)
    if (nrow(stats_df) < 3) {
        log_msg("    Skipping ", dom, ": Not enough CS data points.")
        return(NULL)
    }
    min_cs <- min(stats_df$CS)
    max_cs <- max(stats_df$CS)
    max_taxa <- max(stats_df$Taxa_Count)
    if (max_taxa == min(stats_df$Taxa_Count)) {
        log_msg("    Skipping ", dom, ": Flat curve (no taxa loss).")
        return(NULL)
    }
    calc_df <- stats_df |>
        dplyr::arrange(.data$CS) |>
        dplyr::mutate(
            Pct_Retained = (.data$Taxa_Count / max_taxa) * 100,
            Step_Loss_Pct = dplyr::lag(.data$Pct_Retained, default = 100) -
                .data$Pct_Retained
        )
    n_pts <- nrow(calc_df)
    result <- switch(method,
        kneedle = .ocs_method_kneedle(calc_df, max_cs, dom, log_msg),
        postcliff = .ocs_method_postcliff(calc_df, max_cs, dom, log_msg),
        segmented = .ocs_method_segmented(calc_df, n_pts, dom, log_msg),
        dynamic = .ocs_method_dynamic(calc_df, min_cs, max_cs, dom, log_msg),
        manual = .ocs_method_manual(
            calc_df, min_cs, max_cs, dom, manual_toll, log_msg
        )
    )
    .ocs_tag_audit(calc_df, dom, samp, result$primary_cs, result$sec1_cs)
}

#' Optimal minimum-reads (saturation-curve elbow) for one sample/domain.
#'
#' Operates on the read-count axis in log10 space so the elbow matches the
#' log-scaled saturation plot. Reuses the generic Kneedle / segmented /
#' post-cliff detectors. The CS-tail methods ("dynamic"/"manual") are not
#' meaningful for read counts, so anything other than "postcliff"/"segmented"
#' falls back to Kneedle.
#'
#' @param curve Data frame with columns \code{Cutoff} (>= 1) and
#'   \code{Taxa_Count} (taxa surviving at that cutoff).
#' @param method One of \code{"kneedle"}, \code{"postcliff"}, \code{"segmented"}.
#' @return A list with \code{calc} (the curve with Pct/Step columns),
#'   \code{opt} (optimal cutoff) and \code{sec} (conservative cutoff or NA),
#'   or \code{NULL} if the curve is unusable.
#' @noRd
.si_reads_elbow <- function(curve, method = "kneedle") {
    curve <- curve[curve$Cutoff >= 1, , drop = FALSE]
    curve <- curve[order(curve$Cutoff), , drop = FALSE]
    if (nrow(curve) < 3 || length(unique(curve$Taxa_Count)) < 2) {
        return(NULL)
    }
    max_taxa <- max(curve$Taxa_Count)
    calc_df <- data.frame(
        CS = log10(curve$Cutoff), # x in log space for the elbow geometry
        Cutoff = curve$Cutoff,
        Taxa_Count = curve$Taxa_Count,
        Pct_Retained = curve$Taxa_Count / max_taxa * 100
    )
    calc_df$Step_Loss_Pct <- c(0, -diff(calc_df$Pct_Retained))
    max_x <- max(calc_df$CS)
    n_pts <- nrow(calc_df)
    silent <- function(...) invisible(NULL)
    primary_x <- switch(method,
        postcliff = .ocs_postcliff_cs(calc_df, .ocs_dyn_toll(calc_df), max_x),
        segmented = .ocs_method_segmented(calc_df, n_pts, "", silent)$primary_cs,
        .ocs_knee_cs(calc_df)
    )
    opt_cutoff <- calc_df$Cutoff[which.min(abs(calc_df$CS - primary_x))]
    # Secondary: the conservative post-cliff floor, if it is stricter.
    floor_x <- .ocs_postcliff_cs(calc_df, .ocs_dyn_toll(calc_df), max_x)
    floor_cutoff <- calc_df$Cutoff[which.min(abs(calc_df$CS - floor_x))]
    sec_cutoff <- if (floor_cutoff > opt_cutoff) floor_cutoff else NA_real_
    list(calc = calc_df, opt = opt_cutoff, sec = sec_cutoff)
}

#' Bind per-domain audits, write TSV + RDS, and return the full audit.
#' @noRd
.si_export_audit <- function(audit_list, tax_level, output_dir, log_msg) {
    log_msg(">>> Consolidating Stability Index audit...")
    if (length(audit_list) == 0) {
        log_msg("WARNING: No SI audit data generated. Check inputs.")
        return(invisible(NULL))
    }
    full_audit <- dplyr::bind_rows(audit_list)
    tsv_path <- file.path(output_dir, paste0("SI_Audit_", tax_level, ".tsv"))
    rds_path <- file.path(output_dir, paste0("SI_Audit_", tax_level, ".rds"))
    readr::write_tsv(full_audit, tsv_path)
    readr::write_rds(full_audit, rds_path)
    log_msg("SAVED SI AUDIT TSV: ", tsv_path)
    log_msg("SAVED SI AUDIT RDS: ", rds_path)
    full_audit
}
