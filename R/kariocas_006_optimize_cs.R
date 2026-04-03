#' Optimize Confidence Score using Stability Index (SI) (Step 006)
#'
#' Applies an adaptation of the Kneedle Algorithm to objectively identify the optimal
#' Confidence Score (CS) threshold. It finds the "elbow" point (Primary SI) that
#' maximizes the removal of statistical noise while minimizing the loss of true taxa.
#' It also identifies subsequent local maxima (Secondary SIs) which may represent
#' biological noise separation.
#'
#' @param project_dir Path to the project root.
#' @param tax_level Taxonomic rank to analyze and optimize (default: "Species").
#'
#' @return A data.frame containing the full SI audit trail (invisibly).
#' @export
#' @importFrom dplyr filter select group_by summarise arrange mutate case_when bind_rows pull lead lag n_distinct desc
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_vline geom_text scale_x_continuous scale_y_continuous labs annotate theme element_text coord_cartesian
#' @importFrom patchwork plot_layout plot_annotation
#' @importFrom readr write_tsv write_rds
#' @examples
#' # Get the path to the included toy dataset
#' toy_project <- system.file("extdata", "your_project_name", package = "karioCaS")
#'
#' # Example usage:
#' # Run the Stability Index (SI) optimization for the Species level
#' optimal_thresholds <- optimize_CS(
#'   project_dir = toy_project,
#'   tax_level = "Species"
#' )
#' 
#' # You can also run it for other taxonomic levels, like Genus
#' optimal_thresholds_genus <- optimize_CS(
#'   project_dir = toy_project,
#'   tax_level = "Genus"
#' )

optimize_CS <- function(project_dir, tax_level = "Species") {
  
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
      
      # 3.1 RAW COUNTS AGGREGATION
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
      # 3.2 THE MATH: SI & KNEEDLE ALGORITHM
      # ==========================================================================
      
      calc_df <- stats_df %>%
        dplyr::mutate(
          # A: Percent Retained (For the Plot)
          Pct_Retained = (Taxa_Count / max_taxa) * 100,
          
          # B: Normalization (0 to 1) for X and Y
          Norm_CS = (CS - min_cs) / (max_cs - min_cs),
          Norm_Taxa = (Taxa_Count - min_taxa) / (max_taxa - min_taxa),
          
          # C: Distance 'd' from Secant Line (x + y - 1 = 0)
          Distance = abs(Norm_CS + Norm_Taxa - 1) / sqrt(2)
        )
      
      # D: Find Primary SI (Global Max Distance)
      idx_primary <- which.max(calc_df$Distance)
      primary_cs  <- calc_df$CS[idx_primary]
      
      # E: Find Secondary SIs (Local Maxima strictly AFTER the Primary)
      calc_df <- calc_df %>%
        dplyr::mutate(
          Is_Local_Max = Distance > dplyr::lag(Distance, default = -1) & 
            Distance > dplyr::lead(Distance, default = -1),
          Is_Candidate = CS > primary_cs & Is_Local_Max
        )
      
      secondaries <- calc_df %>% 
        dplyr::filter(Is_Candidate) %>% 
        dplyr::arrange(dplyr::desc(Distance)) %>% 
        head(2)
      
      sec1_cs <- if(nrow(secondaries) >= 1) secondaries$CS[1] else NA
      sec2_cs <- if(nrow(secondaries) >= 2) secondaries$CS[2] else NA
      
      # F: Tag the Audit DataFrame
      calc_df <- calc_df %>%
        dplyr::mutate(
          Domain = dom,
          Sample = samp,
          SI_Type = dplyr::case_when(
            CS == primary_cs ~ "Primary_SI",
            !is.na(sec1_cs) & CS == sec1_cs ~ "Secondary_SI_1",
            !is.na(sec2_cs) & CS == sec2_cs ~ "Secondary_SI_2",
            TRUE ~ NA_character_
          )
        ) %>%
        dplyr::select(Sample, Domain, CS, Taxa_Count, Pct_Retained, Norm_CS, Norm_Taxa, Distance, SI_Type)
      
      audit_list[[length(audit_list) + 1]] <- calc_df
      
      # Log Results
      log_msg(sprintf("    %s -> Primary SI: %02d (Dist: %.3f)", dom, primary_cs, calc_df$Distance[idx_primary]))
      if(!is.na(sec1_cs)) log_msg(sprintf("       -> Sec 1 SI : %02d", sec1_cs))
      if(!is.na(sec2_cs)) log_msg(sprintf("       -> Sec 2 SI : %02d", sec2_cs))
      
      # ==========================================================================
      # 3.3 PLOTTING THE ELBOW
      # ==========================================================================
      
      spec_colors <- get_kariocas_colors("special")
      col_primary <- spec_colors[["Level Taxa"]]
      col_sec     <- spec_colors[["Parent"]]
      labels_vec  <- get_kariocas_labels()
      
      sub_txt <- sprintf("Primary SI: CS %02d", primary_cs)
      if(!is.na(sec1_cs)) sub_txt <- paste0(sub_txt, " | Sec: ", sec1_cs)
      
      p <- ggplot2::ggplot(calc_df, ggplot2::aes(x = CS, y = Pct_Retained)) +
        
        # Secant baseline concept (dotted line connecting first and last point)
        ggplot2::geom_segment(
          x = min_cs, y = 100, xend = max_cs, yend = (min_taxa/max_taxa)*100,
          color = "grey80", linetype = "dotted", linewidth = 0.8
        ) +
        
        # Secondary SI Lines
        {if(!is.na(sec1_cs)) ggplot2::geom_vline(xintercept = sec1_cs, color = col_sec, linetype = "dashed", linewidth = 0.8)} +
        {if(!is.na(sec2_cs)) ggplot2::geom_vline(xintercept = sec2_cs, color = col_sec, linetype = "dashed", linewidth = 0.8)} +
        
        # Primary SI Line
        ggplot2::geom_vline(xintercept = primary_cs, color = col_primary, linetype = "solid", linewidth = 1) +
        
        # Actual Data Curve
        ggplot2::geom_line(color = "black", linewidth = 1) +
        ggplot2::geom_point(color = "black", size = 2.5, shape = get_kariocas_shapes("ranks")[["Level Taxa"]]) +        
        # Primary SI Point Highlight
        ggplot2::geom_point(
          data = calc_df %>% dplyr::filter(CS == primary_cs),
          ggplot2::aes(x = CS, y = Pct_Retained),
          color = col_primary, size = 4, shape = 16
        ) +
        
        # Scales & Labels
        ggplot2::scale_x_continuous(breaks = seq(0, 100, 20), limits = c(0, 100)) +
        ggplot2::scale_y_continuous(limits = c(0, 105)) +
        ggplot2::labs(
          title = dom,
          subtitle = sub_txt,
          x = labels_vec$y_confidence,
          y = "**% Retained**"
        ) +
        theme_kariocas()
      
      plot_list[[dom]] <- p
    }
    
    # 3.4 ASSEMBLE PANEL
    if (length(plot_list) > 0) {
      final_layout <- (plot_list[["Bacteria"]] | plot_list[["Archaea"]]) /
        (plot_list[["Eukaryota"]] | plot_list[["Viruses"]]) +
        patchwork::plot_annotation(
          title = paste(samp, "- CS Optimization (SI)"),
          subtitle = paste("Level:", tax_level, "| Solid Red: Primary SI | Dashed: Secondary SI"),
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
