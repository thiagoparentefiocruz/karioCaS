#' Generate Cutoff Saturation Analysis (Step 003)
#'
#' Generates publication-quality saturation curves showing the retention of reads
#' and taxa across increasing read count cutoffs.
#'
#' @param project_dir Path to the project root (e.g., "/Users/name/project_x").
#' @param analysis_level Taxonomic level to analyze (e.g., "Species", "Genus"). Default: "Species".
#' @param cutoff_mode Analysis mode: "Extinction", "Stabilization", or "Both" (Default).
#' @param import_script_path Deprecated. The package uses internal functions for import.
#'
#' @export
#' @importFrom readr read_rds
#' @importFrom dplyr rename_with filter mutate case_when select bind_rows
#' @importFrom tidyr pivot_longer
#' @importFrom stringr str_detect str_remove str_extract
#' @importFrom ggplot2 ggplot aes geom_line geom_point scale_color_manual scale_shape_manual scale_linetype_manual scale_y_continuous labs theme_classic theme element_text element_line element_blank margin coord_cartesian scale_x_log10 scale_x_continuous
#' @importFrom scales percent_format breaks_log label_number
#' @importFrom patchwork plot_layout plot_annotation
#' @importFrom tidyselect any_of
#' @importFrom stats setNames

reads_per_taxa <- function(project_dir,
                           analysis_level = "Species",
                           cutoff_mode = "Both",
                           import_script_path = NULL) {

  # ==============================================================================
  # 1. SETUP & DIRECTORIES
  # ==============================================================================

  input_dir   <- file.path(project_dir, "000_mpa_original")
  output_dir  <- file.path(project_dir, "003_cutoffs")
  log_dir     <- file.path(project_dir, "log_archives")
  log_file    <- file.path(log_dir, "log_003_cutoff_analysis.txt")

  if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  if(!dir.exists(log_dir))    dir.create(log_dir, recursive = TRUE)

  # ==============================================================================
  # 2. START LOGGING
  # ==============================================================================

  log_con <- file(log_file, open = "wt")
  sink(log_con, type = "output")
  sink(log_con, type = "message")

  on.exit({
    sink(type = "output")
    sink(type = "message")
    close(log_con)
  }, add = TRUE)

  cat("====================================================\n")
  cat("LOG: 003_CUTOFF_ANALYSIS (reads_per_taxa)\n")
  cat("DATE:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("PROJECT DIR:", project_dir, "\n")
  cat("ANALYSIS LEVEL:", analysis_level, "\n")
  cat("CUTOFF MODE:", cutoff_mode, "\n")
  cat("====================================================\n\n")

  tryCatch({

    # ==============================================================================
    # 3. INPUT INTELIGENTE
    # ==============================================================================
    cat("STEP 1: Checking Input Data...\n")

    file_rds <- file.path(input_dir, "combined_mpas_matrix_taxon_level.rds")
    INPUT_MATRIX <- NULL

    if (file.exists(file_rds)) {
      cat("  -> Input file found:", file_rds, "\n")
      cat("  -> Loading RDS object...\n")
      INPUT_MATRIX <- readr::read_rds(file_rds)

    } else {
      cat("  -> [WARNING] Input RDS NOT found.\n")
      cat("  -> Attempting to run Import Function (Step 000)...\n")

      # PACKAGE ADAPTATION: Call internal function directly instead of sourcing
      tryCatch({
        INPUT_MATRIX <- import_mpa_data(project_dir = project_dir)
      }, error = function(e) {
        stop("CRITICAL: Input RDS missing and could not run internal 'import_mpa_data'.\nError: ", e$message)
      })
    }

    cat("\n  -> MATRIX VALIDATION:\n")
    cat("     Observações (Taxa):", nrow(INPUT_MATRIX), "\n")
    cat("     Variáveis (Cols):", ncol(INPUT_MATRIX), "\n")
    cat("     As variáveis são:\n")
    col_names <- colnames(INPUT_MATRIX)
    for (i in 1:length(col_names)) {
      cat(sprintf("     [%d] %s\n", i, col_names[i]))
    }
    cat("\n")

    # ==============================================================================
    # 4. DATA PROCESSING
    # ==============================================================================
    cat("STEP 2: Processing Data...\n")

    df_long <- INPUT_MATRIX %>%
      dplyr::rename_with(~ "Taxonomy",  .cols = tidyselect::any_of(c("Taxonomia", "Taxonomy"))) %>%
      dplyr::rename_with(~ "Tax_Level", .cols = tidyselect::any_of(c("tax_level", "Tax_Level", "Nivel_Final"))) %>%
      tidyr::pivot_longer(cols = -c(Taxonomy, Tax_Level), names_to = "Raw_Sample_Col", values_to = "Counts") %>%
      dplyr::filter(Counts > 0) %>%
      dplyr::mutate(
        Domain = dplyr::case_when(
          stringr::str_detect(Taxonomy, "d__Viruses")   ~ "Viruses",
          stringr::str_detect(Taxonomy, "d__Bacteria")  ~ "Bacteria",
          stringr::str_detect(Taxonomy, "d__Archaea")   ~ "Archaea",
          stringr::str_detect(Taxonomy, "d__Eukaryota") ~ "Eukaryota",
          TRUE ~ "Unknown"
        ),
        Sample_Name = stringr::str_remove(Raw_Sample_Col, "_CS[0-9.]+$"),
        CS_Label    = stringr::str_extract(Raw_Sample_Col, "CS[0-9.]+$")
      ) %>%
      dplyr::filter(Tax_Level == analysis_level)

    unique_samples <- unique(df_long$Sample_Name)
    cat("  -> Samples to process:", length(unique_samples), "\n")

    # ==============================================================================
    # 5. HELPER FUNCTIONS
    # ==============================================================================

    generate_master_cutoffs <- function() {
      c(1:15, seq(16, 30, 2), seq(35, 100, 5), seq(110, 500, 10), seq(550, 2000, 50), seq(2100, 10000, 100))
    }
    MASTER_CUTOFFS <- generate_master_cutoffs()

    calc_retention_curve <- function(sub_df, cutoffs) {
      total_reads_abs <- sum(sub_df$Counts)
      total_taxa_abs  <- nrow(sub_df)
      if(total_reads_abs == 0) return(NULL)

      results <- list()
      max_reads_in_sample <- max(sub_df$Counts)

      for (cut in cutoffs) {
        if (cut > max_reads_in_sample) {
          results[[length(results)+1]] <- data.frame(Cutoff = cut, Pct_Reads = 0, Pct_Taxa = 0, Abs_Taxa = 0)
          break
        }
        passed <- sub_df %>% dplyr::filter(Counts >= cut)
        rem_reads <- sum(passed$Counts)
        rem_taxa  <- nrow(passed)
        results[[length(results)+1]] <- data.frame(
          Cutoff = cut,
          Pct_Reads = (rem_reads / total_reads_abs) * 100,
          Pct_Taxa  = (rem_taxa / total_taxa_abs) * 100,
          Abs_Taxa  = rem_taxa
        )
      }
      res_df <- dplyr::bind_rows(results)
      res_df$Abs_Reads_Total <- total_reads_abs
      res_df$Abs_Taxa_Total  <- total_taxa_abs
      return(res_df)
    }

    create_dummy_plot <- function(dom_name, level_name) {
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 50, label = "No Data", color = "grey80", size = 5) +
        ggplot2::labs(
          title = dom_name,
          subtitle = paste("0 Reads | 0", level_name),
          x = "Reads", y = "% Retained"
        ) +
        ggplot2::scale_y_continuous(limits = c(0, 100)) +
        ggplot2::scale_x_continuous(limits = c(0, 1)) +
        ggplot2::theme_classic(base_size = 10) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 11, color = "black"),
          plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 9, color = "grey30"),
          axis.line = ggplot2::element_line(color = "grey80"),
          axis.text = ggplot2::element_text(color = "grey80"),
          axis.ticks = ggplot2::element_line(color = "grey80"),
          axis.title = ggplot2::element_text(color = "grey80"),
          panel.grid = ggplot2::element_blank()
        )
    }

    # ==============================================================================
    # 6. GENERATION LOOP
    # ==============================================================================
    cat("STEP 3: Generating Curves...\n")

    DOMAINS_ORDER <- c("Bacteria", "Archaea", "Eukaryota", "Viruses")
    COLOR_MAP    <- stats::setNames(c("#E64B35", "#4DBBD5"), c("Reads", analysis_level))
    SHAPE_MAP    <- stats::setNames(c(16, 17), c("Reads", analysis_level))
    LINETYPE_MAP <- stats::setNames(c("solid", "dashed"), c("Reads", analysis_level))

    modes_to_run <- if(cutoff_mode == "Both") c("Extinction", "Stabilization") else cutoff_mode

    for (samp in unique_samples) {
      cat("  Processing Sample:", samp, "...\n")

      samp_out_dir <- file.path(output_dir, samp)
      if(!dir.exists(samp_out_dir)) dir.create(samp_out_dir)

      df_samp <- df_long %>% dplyr::filter(Sample_Name == samp)
      unique_cs <- sort(unique(df_samp$CS_Label))

      for (cs in unique_cs) {

        df_cs <- df_samp %>% dplyr::filter(CS_Label == cs)

        for (current_mode in modes_to_run) {

          plot_list <- list()

          for (dom in DOMAINS_ORDER) {

            df_dom <- df_cs %>% dplyr::filter(Domain == dom)

            if (nrow(df_dom) == 0) {
              plot_list[[dom]] <- create_dummy_plot(dom, analysis_level)
              next
            }

            # --- CALCULATION ---
            curve_data <- calc_retention_curve(df_dom, MASTER_CUTOFFS)

            if (is.null(curve_data) || nrow(curve_data) == 0) {
              plot_list[[dom]] <- create_dummy_plot(dom, analysis_level)
              next
            }

            # --- TRIMMING LOGIC ---
            limit_idx <- nrow(curve_data)

            if (current_mode == "Extinction") {
              zero_idx <- which(curve_data$Abs_Taxa == 0)[1]
              if (!is.na(zero_idx)) limit_idx <- min(nrow(curve_data), zero_idx + 2)

            } else if (current_mode == "Stabilization") {
              taxa_vals <- curve_data$Abs_Taxa
              if(length(taxa_vals) > 1) {
                ratios <- taxa_vals[2:length(taxa_vals)] / taxa_vals[1:(length(taxa_vals)-1)]
                ratios <- c(0, ratios)
                stable_indices <- which(ratios >= 0.95 & seq_along(ratios) > 2)

                if (length(stable_indices) > 0) {
                  limit_idx <- min(nrow(curve_data), stable_indices[1] + 2)
                } else {
                  zero_idx <- which(curve_data$Abs_Taxa == 0)[1]
                  if (!is.na(zero_idx)) limit_idx <- min(nrow(curve_data), zero_idx + 2)
                }
              }
            }

            plot_data_final <- curve_data[1:limit_idx, ]

            # --- PLOT PREP ---
            abs_reads <- unique(plot_data_final$Abs_Reads_Total)
            abs_taxa  <- unique(plot_data_final$Abs_Taxa_Total)

            subtitle_str <- paste0(format(abs_reads, big.mark = ",", scientific = FALSE), " Reads | ",
                                   format(abs_taxa, big.mark = ",", scientific = FALSE), " ", analysis_level)

            df_gg <- plot_data_final %>%
              dplyr::select(Cutoff, Pct_Reads, Pct_Taxa) %>%
              tidyr::pivot_longer(cols = c(Pct_Reads, Pct_Taxa), names_to = "Metric_Type", values_to = "Percentage") %>%
              dplyr::mutate(Metric_Label = ifelse(Metric_Type == "Pct_Reads", "Reads", analysis_level))

            # --- CONDITIONAL SCALE (USER FIX APPLIED HERE) ---
            max_x_val <- max(plot_data_final$Cutoff, na.rm=TRUE)
            use_log_scale <- (max_x_val >= 1000)

            p <- ggplot2::ggplot(df_gg, ggplot2::aes(x = Cutoff, y = Percentage, color = Metric_Label, group = Metric_Label)) +
              ggplot2::geom_line(ggplot2::aes(linetype = Metric_Label), linewidth = 0.8) +
              ggplot2::geom_point(ggplot2::aes(shape = Metric_Label), size = 2) +
              ggplot2::scale_color_manual(values = COLOR_MAP) +
              ggplot2::scale_shape_manual(values = SHAPE_MAP) +
              ggplot2::scale_linetype_manual(values = LINETYPE_MAP) +
              ggplot2::scale_y_continuous(labels = scales::percent_format(scale = 1), limits = c(0, 105), expand = c(0,0)) +
              ggplot2::labs(title = dom, subtitle = subtitle_str, y = "% Retained", color = NULL, shape = NULL, linetype = NULL) +
              ggplot2::theme_classic(base_size = 10) +
              ggplot2::theme(
                plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 11),
                plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 9, color = "grey30"),
                axis.text = ggplot2::element_text(color = "black"),
                axis.line = ggplot2::element_line(linewidth = 0.4, color = "black"),
                legend.position = "bottom",
                panel.grid.major.y = ggplot2::element_line(color = "grey90", linetype = "dotted"),
                plot.margin = ggplot2::margin(t = 10, r = 15, b = 10, l = 10)
              ) +
              ggplot2::coord_cartesian(clip = "off")

            # --------------------------------------------------------------------------
            # APPLIED USER FIX FOR X-AXIS
            # --------------------------------------------------------------------------
            if (use_log_scale) {
              p <- p +
                ggplot2::scale_x_log10(
                  breaks = scales::breaks_log(n = 10, base = 10),   # breaks em potências de 10
                  labels = scales::label_number(accuracy = 1, big.mark = "", trim = TRUE)
                ) +
                ggplot2::labs(x = "Reads (log10 scale)")
            } else {
              x_vals <- plot_data_final$Cutoff
              if (length(x_vals) > 15) {
                x_breaks <- x_vals[seq(1, length(x_vals), length.out = 10)]
                x_breaks <- unique(round(x_breaks))
              } else {
                x_breaks <- x_vals
              }
              p <- p +
                ggplot2::scale_x_continuous(breaks = x_breaks, labels = x_breaks) +
                ggplot2::labs(x = "Reads")
            }
            # --------------------------------------------------------------------------

            plot_list[[dom]] <- p
          }

          # --- SAVE ---
          final_layout <- (plot_list[["Bacteria"]] | plot_list[["Archaea"]]) /
            (plot_list[["Eukaryota"]] | plot_list[["Viruses"]]) +
            patchwork::plot_layout(guides = "collect") &
            ggplot2::theme(legend.position = "bottom", legend.text = ggplot2::element_text(size = 10, face = "bold")) &
            patchwork::plot_annotation(
              title = paste(samp, "-", cs, "| Saturation:", current_mode),
              subtitle = paste("Retention of Reads vs", analysis_level, "with", current_mode, "focus"),
              theme = ggplot2::theme(plot.title = ggplot2::element_text(face="bold", size=14, hjust=0.5))
            )

          file_name <- paste0(samp, "_", cs, "_Cutoff_", analysis_level, "_", current_mode, ".pdf")
          out_path <- file.path(samp_out_dir, file_name)
          ggplot2::ggsave(out_path, final_layout, width = 11.69, height = 8.27, units = "in")
          cat("    -> Saved [", current_mode, "]:", file_name, "\n")
        }
      }
    }

    cat("\n====================================================\n")
    cat("SUCCESS: All Saturation plots generated.\n")
    cat("====================================================\n")

    return(TRUE)

  }, error = function(e) {
    cat("\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
    cat("CRITICAL ERROR DETECTED:\n")
    cat(e$message, "\n")
    cat("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
    stop(e$message)
  })
}
