#' Optimize Confidence Score using Multi-Strategy Stability Index (SI) (Step 006)
#'
#' Objectively identifies the optimal Confidence Score (CS) threshold to separate 
#' statistical noise from true biological signal. It replaces legacy geometric approaches 
#' with a versatile multi-strategy mathematical engine:
#' - "dynamic": Adapts to the baseline noise of the sample (stricter, ideal for pathogens).
#' - "segmented": Uses a Broken-Stick regression model to find regime shifts (ideal for ecology/dark matter).
#' - "manual": Allows expert-defined acceptable loss tolls globally or per domain.
#'
#' @param project_dir Character. Path to the project root.
#' @param tax_level Character. Taxonomic rank to analyze and optimize (default: "Species").
#' @param method Character. Strategy to calculate the SI: "dynamic" (default), "segmented", or "manual".
#' @param manual_toll Numeric or named list. Acceptable step-wise loss percentage. Used only if method = "manual" (default: 1.0).
#'
#' @return A data.frame containing the full SI audit trail (invisibly).
#' @export
#' @importFrom dplyr filter select group_by summarise arrange mutate case_when bind_rows lag n_distinct
#' @importFrom ggplot2 ggplot aes geom_line geom_point scale_x_continuous scale_y_continuous labs theme element_text
#' @importFrom patchwork plot_layout plot_annotation
#' @importFrom readr write_tsv write_rds
#' @importFrom stats lm resid sd
#' @examples
#' # Get the path to the included toy dataset
#' toy_project <- system.file("extdata", "your_project_name", package = "karioCaS")
#'
#' # Default usage (Dynamic Toll - Ideal for clinical/pathogen focus)
#' optimal_dynamic <- optimize_CS(
#'   project_dir = toy_project,
#'   tax_level = "Species"
#' )
#' 
#' # Segmented Regression (Ideal for ecology and biological dark matter)
#' optimal_segmented <- optimize_CS(
#'   project_dir = toy_project,
#'   tax_level = "Species",
#'   method = "segmented"
#' )

optimize_CS <- function(project_dir, tax_level = "Species", method = c("dynamic", "segmented", "manual"), manual_toll = 1.0) {
  
  method <- match.arg(method) # Garante que o usuário digitou uma das 3 opções válidas
  
  # ==============================================================================
  # 1. SETUP & LOGGING
  # ==============================================================================
  output_dir <- file.path(project_dir, "006_optimize_CS")
  log_dir    <- file.path(project_dir, "logs")
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)
  
  log_file <- file.path(log_dir, "log_006_optimize_cs.txt")
  
  log_msg <- function(...) {
    msg <- paste0(...)
    message(msg)
    cat(paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", msg, "\n"),
        file = log_file, append = TRUE)
  }
  
  cat("====================================================\n", file = log_file)
  cat("LOG: 006_OPTIMIZE_CS (Stability Index - SI)\n", file = log_file, append = TRUE)
  cat("PROJECT DIR: ", project_dir, "\n", file = log_file, append = TRUE)
  cat("TAXA LEVEL:  ", tax_level, "\n", file = log_file, append = TRUE)
  cat("METHOD:      ", toupper(method), "\n", file = log_file, append = TRUE)
  cat("====================================================\n", file = log_file, append = TRUE)
  
  # ==============================================================================
  # 2. DATA LOADING & PREP
  # ==============================================================================
  log_msg(">>> Loading Data (Auto-detected format)...")
  df_long <- .get_tidy_data(project_dir)
  
  df_proc <- df_long %>% dplyr::filter(Rank == tax_level)
  
  if (nrow(df_proc) == 0) {
    log_msg("CRITICAL ERROR: No data found for Rank: ", tax_level)
    stop("No data for specified rank.")
  }
  
  SAMPLES <- unique(df_proc$sample)
  DOMAINS <- names(get_kariocas_colors("domains"))
  
  audit_list <- list()
  
  log_msg(">>> Starting SI Calculation for ", length(SAMPLES), " samples.")
  
  # ==============================================================================
  # 3. ANALYSIS LOOP
  # ==============================================================================
  for (samp in SAMPLES) {
    log_msg("------------------------------------------------")
    log_msg("  Processing Sample: ", samp)
    
    plot_list <- list()
    df_samp <- df_proc %>% dplyr::filter(sample == samp)
    
    for (dom in DOMAINS) {
      
      df_dom <- df_samp %>% dplyr::filter(Domain == dom)
    
  # ==========================================================================
  # 3.1 RAW COUNTS AGGREGATION
  # ==========================================================================
      
      stats_df <- df_dom %>%
        dplyr::group_by(CS) %>%
        dplyr::summarise(Taxa_Count = dplyr::n_distinct(Taxon_Name), .groups = "drop") %>%
        dplyr::arrange(CS)
      
      if (nrow(stats_df) < 3) {
        log_msg("    Skipping ", dom, ": Not enough CS data points for curve fitting.")
        plot_list[[dom]] <- plot_kariocas_empty(dom, "Insufficient CS points")
        next
      }
      
      min_cs <- min(stats_df$CS)
      max_cs <- max(stats_df$CS)
      max_taxa <- max(stats_df$Taxa_Count)
      min_taxa <- min(stats_df$Taxa_Count)
      
      if (max_taxa == min_taxa) {
        log_msg("    Skipping ", dom, ": Flat curve (no taxa loss).")
        plot_list[[dom]] <- plot_kariocas_empty(dom, "Flat Curve (No elbow)")
        next
      }

  # ==========================================================================
  # 3.2 THE MATH: MULTI-STRATEGY ENGINE
  # ==========================================================================
      
      calc_df <- stats_df %>%
        dplyr::arrange(CS) %>%
        dplyr::mutate(
          Pct_Retained = (Taxa_Count / max_taxa) * 100,
          Step_Loss_Pct = dplyr::lag(Pct_Retained, default = 100) - Pct_Retained
        )
      
      n_pts <- nrow(calc_df)
      sec1_cs <- NA # Inicializa vazio
      
      if (method == "segmented") {
        # --- ESTRATÉGIA 1: REGRESSÃO SEGMENTADA (ECOLOGIA / DARK MATTER) ---
        if (n_pts >= 5) {
          best_rss <- Inf
          primary_idx <- 2
          for (i in 3:(n_pts - 2)) {
            fit_left  <- lm(Pct_Retained ~ CS, data = calc_df[1:i, ])
            fit_right <- lm(Pct_Retained ~ CS, data = calc_df[i:n_pts, ])
            rss_total <- sum(resid(fit_left)^2) + sum(resid(fit_right)^2)
            
            if (rss_total < best_rss) {
              best_rss <- rss_total
              primary_idx <- i
            }
          }
          primary_cs <- calc_df$CS[primary_idx]
        } else {
          primary_cs <- calc_df$CS[min(2, n_pts)]
        }
        sub_txt <- sprintf("Method: Segmented | Shift Breakpoint: CS %02d", primary_cs)
        log_msg(sprintf("    %s -> Regime Shift: CS %02d", dom, primary_cs))
        
      } else if (method == "dynamic") {
        # --- ESTRATÉGIA 2: DYNAMIC TOLL (PATÓGENOS / RIGOR OBJETIVO) ---
        tail_df <- calc_df %>% dplyr::filter(CS >= 50)
        
        if (nrow(tail_df) >= 2) {
          bg_mean <- mean(tail_df$Step_Loss_Pct, na.rm = TRUE)
          bg_sd   <- sd(tail_df$Step_Loss_Pct, na.rm = TRUE)
          dynamic_toll <- max(bg_mean + (1.5 * bg_sd), 0.5)
        } else {
          dynamic_toll <- 1.0
        }
        
        stable_points <- calc_df %>% dplyr::filter(CS > min_cs & Step_Loss_Pct <= dynamic_toll)
        primary_cs <- if(nrow(stable_points) > 0) stable_points$CS[1] else max_cs
        
        ultra_toll <- if(nrow(tail_df) >= 2) max(bg_mean, 0.1) else 0.2
        ultra_stable <- calc_df %>% dplyr::filter(CS > primary_cs & Step_Loss_Pct <= ultra_toll)
        sec1_cs <- if(nrow(ultra_stable) >= 1) ultra_stable$CS[1] else NA
        
        sub_txt <- sprintf("Method: Dynamic (Toll: %.2f%%) | Primary: CS %02d | Sec: %s", 
                           dynamic_toll, primary_cs, ifelse(is.na(sec1_cs), "NA", sprintf("%02d", sec1_cs)))
        log_msg(sprintf("    %s -> Dynamic Toll: %.2f%% | Primary: CS %02d", dom, dynamic_toll, primary_cs))
        
      } else if (method == "manual") {
        # --- ESTRATÉGIA 3: MANUAL TOLL (CONTROLE DO ESPECIALISTA) ---
        # Verifica se manual_toll é uma lista com valor específico para o domínio, senão usa o valor numérico
        current_toll <- if(is.list(manual_toll) && !is.null(manual_toll[[dom]])) {
          manual_toll[[dom]]
        } else if (is.numeric(manual_toll) && length(manual_toll) == 1) {
          manual_toll
        } else {
          1.0 # Fallback de segurança
        }
        
        stable_points <- calc_df %>% dplyr::filter(CS > min_cs & Step_Loss_Pct <= current_toll)
        primary_cs <- if(nrow(stable_points) > 0) stable_points$CS[1] else max_cs
        
        ultra_stable <- calc_df %>% dplyr::filter(CS > primary_cs & Step_Loss_Pct <= 0.2)
        sec1_cs <- if(nrow(ultra_stable) >= 1) ultra_stable$CS[1] else NA
        
        sub_txt <- sprintf("Method: Manual (Toll: %.2f%%) | Primary: CS %02d | Sec: %s", 
                           current_toll, primary_cs, ifelse(is.na(sec1_cs), "NA", sprintf("%02d", sec1_cs)))
        log_msg(sprintf("    %s -> Manual Toll: %.2f%% | Primary: CS %02d", dom, current_toll, primary_cs))
      }
      
      # Tag the Audit DataFrame
      calc_df <- calc_df %>%
        dplyr::mutate(
          Domain = dom,
          Sample = samp,
          SI_Type = dplyr::case_when(
            CS == primary_cs ~ "Primary_SI",
            !is.na(sec1_cs) & CS == sec1_cs ~ "Secondary_SI_1",
            TRUE ~ NA_character_
          )
        ) %>%
        dplyr::select(Sample, Domain, CS, Taxa_Count, Pct_Retained, Step_Loss_Pct, SI_Type)
      
      audit_list[[length(audit_list) + 1]] <- calc_df
      
  # ==========================================================================
  # 3.3 PLOTTING THE ELBOW/SHIFT
  # ==========================================================================
      
      spec_colors <- get_kariocas_colors("special")
      col_primary <- spec_colors[["Level Taxa"]]
      labels_vec  <- get_kariocas_labels()
      
      p <- ggplot2::ggplot(calc_df, ggplot2::aes(x = CS, y = Pct_Retained)) +
        
        # Actual Data Curve
        ggplot2::geom_line(color = "black", linewidth = 1) +
        ggplot2::geom_point(color = "black", size = 2.5, shape = get_kariocas_shapes("ranks")[["Level Taxa"]]) +        
        
        # SI Points Highlight
        ggplot2::geom_point(
          data = calc_df %>% dplyr::filter(!is.na(SI_Type)),
          ggplot2::aes(x = CS, y = Pct_Retained),
          color = col_primary, size = 4, shape = 4, stroke = 1.5
        ) +
        
        # Scales & Labels
        ggplot2::scale_x_continuous(breaks = seq(0, 100, 20), limits = c(0, 100)) +
        ggplot2::scale_y_continuous(limits = c(0, 105)) +
        ggplot2::labs(
          title = dom,
          subtitle = sub_txt, # Agora varia conforme a estratégia escolhida!
          x = labels_vec$y_confidence,
          y = "**% Retained**"
        ) +
        theme_kariocas()
      
      plot_list[[dom]] <- p
    }      

  # ==========================================================================
  # 3.4 ASSEMBLE PANEL
  # ==========================================================================

    if (length(plot_list) > 0) {
      final_layout <- (plot_list[["Bacteria"]] | plot_list[["Archaea"]]) /
        (plot_list[["Eukaryota"]] | plot_list[["Viruses"]]) +
        patchwork::plot_annotation(
          title = paste(samp, "- CS Optimization (SI)"),
          subtitle = paste("Level:", tax_level, "| Red 'X': Stability Indices (Primary & Secondary)"),
          theme = ggplot2::theme(
            plot.title = ggplot2::element_text(face = "bold", size = 16, hjust = 0.5),
            plot.subtitle = ggplot2::element_text(size = 12, hjust = 0.5, color = "grey30")
          )
        )
      
      file_name <- paste0(samp, "_Optimize_CS_", tax_level, ".pdf")
      save_path <- file.path(output_dir, file_name)
      
      dims <- get_kariocas_dims()
      ggplot2::ggsave(save_path, final_layout, width = dims$width, height = dims$height)
      log_msg("    -> Generated Plot: ", file_name)
    }
  }

  # ==============================================================================
  # 4. EXPORT AUDIT DATAFRAME
  # ==============================================================================
  log_msg(">>> Consolidating Audit Data...")
  
  if (length(audit_list) > 0) {
    full_audit_df <- dplyr::bind_rows(audit_list)
    
    tsv_path <- file.path(output_dir, paste0("SI_Audit_", tax_level, ".tsv"))
    rds_path <- file.path(output_dir, paste0("SI_Audit_", tax_level, ".rds"))
    
    readr::write_tsv(full_audit_df, tsv_path)
    readr::write_rds(full_audit_df, rds_path)
    
    log_msg("SAVED AUDIT TSV: ", tsv_path)
    log_msg("SAVED AUDIT RDS: ", rds_path)
  } else {
    log_msg("WARNING: No audit data generated. Check inputs.")
    return(invisible(NULL))
  }
  
  log_msg("SUCCESS: CS Optimization completed.")
  return(invisible(full_audit_df))
}
