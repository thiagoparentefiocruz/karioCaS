#' Run Confidence Score Retention Analysis (Step 001)
#'
#' Executes read retention analysis based on Confidence Score (Kraken/Bracken).
#' This function expects the standard folder structure created by step 000.
#' It generates 7 PDFs per sample (1 "All Levels" + 6 Taxonomic Levels).
#'
#' @param project_dir Root path of the project (e.g., "/Users/name/project_x").
#' @param import_script_path Deprecated/Ignored in package version. The import function is now internal.
#'
#' @return Generates PDF files in "001_cs_retention" and returns invisible NULL.
#' @export
#' @importFrom readr read_rds
#' @importFrom dplyr rename mutate filter select group_by summarise n left_join bind_rows pull case_when distinct any_of %>%
#' @importFrom tidyr pivot_longer replace_na
#' @importFrom stringr str_detect str_remove str_extract
#' @importFrom ggplot2 ggplot aes geom_line geom_point scale_color_manual scale_shape_manual scale_linetype_manual scale_y_log10 scale_x_continuous labs theme_classic guides guide_legend theme element_text element_line element_blank margin coord_cartesian annotate ggsave
#' @importFrom patchwork plot_layout plot_annotation wrap_plots
#' @import patchwork

taxa_retention <- function(project_dir, import_script_path = NULL) {

  # ==============================================================================
  # 1. SETUP & DIRECTORIES
  # ==============================================================================

  input_dir   <- file.path(project_dir, "000_mpa_original")
  output_dir  <- file.path(project_dir, "001_cs_retention")
  log_dir     <- file.path(project_dir, "log_archives")
  log_file    <- file.path(log_dir, "log_001_cs_retention.txt")

  if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  if(!dir.exists(log_dir))    dir.create(log_dir, recursive = TRUE)

  # ==============================================================================
  # 2. LOGGING SETUP
  # ==============================================================================

  # Initialize Log
  sink(log_file, append = FALSE)

  # Helper for log/console
  log_msg <- function(...) {
    cat(paste0("[", format(Sys.time(), "%H:%M:%S"), "] "), ..., "\n")
  }

  log_msg("INIT: Starting CS Retention Analysis (Step 001)")
  log_msg("Project Directory:", project_dir)

  # TryCatch block ensures sink() closes on error
  tryCatch({

    # ==============================================================================
    # 3. INTELLIGENT INPUT (LOAD OR RUN STEP 000)
    # ==============================================================================
    log_msg("STEP 1: Checking Input Data...")
    file_rds <- file.path(input_dir, "combined_mpas_matrix_taxon_level.rds")

    INPUT_MATRIX <- NULL

    if (file.exists(file_rds)) {
      log_msg("  -> Loading RDS:", file_rds)
      INPUT_MATRIX <- readr::read_rds(file_rds)
    } else {
      log_msg("  -> RDS not found. Calling internal Import Function (Step 000)...")

      # PACKAGE ADAPTATION:
      # Instead of sourcing a file, we call the package function directly.
      # We assume import_mpa_data is exported or available in the namespace.
      tryCatch({
        INPUT_MATRIX <- import_mpa_data(project_dir = project_dir)
      }, error = function(e) {
        stop("CRITICAL: Could not run 'import_mpa_data'. Please ensure the previous step logic is correct.\nError: ", e$message)
      })
    }

    log_msg("  -> Data Loaded. Rows:", nrow(INPUT_MATRIX))

    # ==============================================================================
    # 4. DATA PROCESSING
    # ==============================================================================
    log_msg("STEP 2: Processing Data Structure...")

    # Column mapping
    cols_map <- c("Taxonomy" = "Taxonomia", "Taxonomy" = "Taxonomy",
                  "Tax_Level" = "tax_level", "Tax_Level" = "Tax_Level", "Tax_Level" = "Nivel_Final")

    data_renamed <- INPUT_MATRIX %>%
      dplyr::rename(dplyr::any_of(cols_map))

    if (!"Taxonomy" %in% names(data_renamed) || !"Tax_Level" %in% names(data_renamed)) {
      stop("ERROR: Input data missing required columns (Taxonomy/Tax_Level).")
    }

    df_long <- data_renamed %>%
      tidyr::pivot_longer(cols = -c(Taxonomy, Tax_Level), names_to = "Raw_Sample_Col", values_to = "Counts") %>%
      dplyr::mutate(
        Domain = dplyr::case_when(
          stringr::str_detect(Taxonomy, "d__Viruses")   ~ "Viruses",
          stringr::str_detect(Taxonomy, "d__Bacteria")  ~ "Bacteria",
          stringr::str_detect(Taxonomy, "d__Archaea")   ~ "Archaea",
          stringr::str_detect(Taxonomy, "d__Eukaryota") ~ "Eukaryota",
          TRUE ~ "Unknown"
        ),
        Sample_Name = stringr::str_remove(Raw_Sample_Col, "_CS[0-9.]+$"),
        CS_Label    = stringr::str_extract(Raw_Sample_Col, "CS[0-9.]+$"),
        CS_Num      = as.numeric(stringr::str_extract(CS_Label, "\\d+")) * 10
      )

    unique_samples <- unique(df_long$Sample_Name)
    min_cs_val <- min(df_long$CS_Num, na.rm=TRUE)

    df_domains_only <- df_long %>% dplyr::filter(Tax_Level == "Domain")

    # Global Filter (Remove Kingdom)
    df_taxa_only <- df_long %>%
      dplyr::filter(
        Tax_Level != "Domain",
        Tax_Level != "Kingdom",
        !is.na(Tax_Level),
        Tax_Level != "unclassified",
        Counts > 0
      )

    # ==============================================================================
    # 5. AESTHETICS & PLOTTING
    # ==============================================================================
    log_msg("STEP 3: Generating Plots...")

    PLURAL_ORDER <- c("Phyla", "Classes", "Orders", "Families", "Genera", "Species")
    PLURAL_MAP <- c(
      "Phylum"  = "Phyla", "Class"   = "Classes", "Order"   = "Orders",
      "Family"  = "Families","Genus"   = "Genera", "Species" = "Species"
    )
    TAX_COLORS <- c(
      "Phyla"    = "#4DBBD5", "Classes"  = "#00A087", "Orders"   = "#3C5488",
      "Families" = "#F39B7F", "Genera"   = "#8491B4", "Species"  = "#91D1C2",
      "Total Reads" = "#7E6148", "Assigned Reads" = "#DC0000"
    )
    TAX_SHAPES <- c(
      "Phyla"    = 16, "Classes"  = 17, "Orders"   = 18,
      "Families" = 8,  "Genera"   = 3, "Species"  = 4,
      "Total Reads" = 15, "Assigned Reads" = 17
    )

    LEVELS_TO_PROCESS <- list(NULL, "Species", "Genus", "Family", "Order", "Class", "Phylum")
    DOMAINS_ORDER <- c("Bacteria", "Archaea", "Eukaryota", "Viruses")
    y_axis_label <- "% Retained\n(axis scaled to log10)"

    # Internal helper for empty plots
    create_dummy_plot <- function(dom_name, legend_levels = NULL, l_colors = NULL, l_shapes = NULL) {
      p <- ggplot2::ggplot() +
        ggplot2::annotate("text", x = 50, y = 1, label = "No Data", color = "grey80", size = 5) +
        ggplot2::labs(title = dom_name, subtitle = "Reads: 0", x = "Confidence Score (%)", y = y_axis_label) +
        ggplot2::scale_y_log10(limits = c(0.01, 100), breaks = c(0.01, 0.1, 1, 10, 100), labels = c("0.01", "0.1", "1", "10", "100")) +
        ggplot2::scale_x_continuous(limits = c(0, 100)) +
        ggplot2::theme_classic(base_size = 10) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(face="bold", hjust=0.5, size=11),
          plot.subtitle = ggplot2::element_text(hjust=0.5, size=9, color="grey30"),
          axis.line = ggplot2::element_line(color="grey80"),
          axis.text = ggplot2::element_text(color="grey80"),
          axis.title = ggplot2::element_text(color="grey80"),
          panel.grid = ggplot2::element_blank()
        )

      if(!is.null(legend_levels)) {
        dummy_df <- data.frame(
          CS_Num = 0, Percentage = 0,
          Legend_Label = factor(legend_levels, levels = legend_levels)
        )
        p <- p +
          ggplot2::geom_point(data = dummy_df, ggplot2::aes(x=CS_Num, y=Percentage, color=Legend_Label, shape=Legend_Label), alpha=0) +
          ggplot2::scale_color_manual(values = l_colors, drop = FALSE) +
          ggplot2::scale_shape_manual(values = l_shapes, drop = FALSE) +
          ggplot2::guides(color = ggplot2::guide_legend(nrow=1), shape = ggplot2::guide_legend(nrow=1))
      }
      return(p)
    }

    # --- LOOP Samples ---
    for (samp in unique_samples) {
      log_msg(paste("  Processing Sample:", samp))

      samp_out_dir <- file.path(output_dir, samp)
      if(!dir.exists(samp_out_dir)) dir.create(samp_out_dir)

      # --- LOOP Modes ---
      for (iter_level in LEVELS_TO_PROCESS) {

        ANALYSIS_LEVEL <- iter_level
        mode_tag <- ifelse(is.null(ANALYSIS_LEVEL), "All_Levels", PLURAL_MAP[ANALYSIS_LEVEL])
        if(is.na(mode_tag)) mode_tag <- ANALYSIS_LEVEL

        plot_list <- list()
        specific_levels_order <- NULL

        if (!is.null(ANALYSIS_LEVEL)) {
          lbl_plural <- PLURAL_MAP[ANALYSIS_LEVEL]
          lbl_total  <- "Total Reads"
          lbl_assign <- paste("Reads assigned to", ANALYSIS_LEVEL)
          specific_levels_order <- c(lbl_plural, lbl_total, lbl_assign)
        }

        for (dom in DOMAINS_ORDER) {
          taxa_dom  <- df_taxa_only %>% dplyr::filter(Sample_Name == samp, Domain == dom)
          dom_stats <- df_domains_only %>% dplyr::filter(Sample_Name == samp, Domain == dom)
          baseline_reads_row <- dom_stats %>% dplyr::filter(CS_Num == min_cs_val)

          # Case: No Data
          if(nrow(baseline_reads_row) == 0) {
            if (!is.null(ANALYSIS_LEVEL)) {
              plot_list[[dom]] <- create_dummy_plot(dom, specific_levels_order, TAX_COLORS, TAX_SHAPES)
            } else {
              plot_list[[dom]] <- create_dummy_plot(dom, PLURAL_ORDER, TAX_COLORS, TAX_SHAPES)
            }
            next
          }

          TOTAL_READS_BASELINE <- baseline_reads_row$Counts[1]

          # MODE 1: Specific Level
          if (!is.null(ANALYSIS_LEVEL)) {
            MODE1_COLORS <- TAX_COLORS; MODE1_COLORS[lbl_assign] <- TAX_COLORS[["Assigned Reads"]]
            MODE1_SHAPES <- TAX_SHAPES; MODE1_SHAPES[lbl_assign] <- TAX_SHAPES[["Assigned Reads"]]

            curve_total_reads <- dom_stats %>%
              dplyr::mutate(Percentage = (Counts / TOTAL_READS_BASELINE) * 100, Legend_Label = "Total Reads") %>%
              dplyr::select(CS_Num, Percentage, Legend_Label)

            spec_level_data <- taxa_dom %>% dplyr::filter(Tax_Level == ANALYSIS_LEVEL)

            if(nrow(spec_level_data) > 0) {
              spec_reads_sums <- spec_level_data %>% dplyr::group_by(CS_Num) %>% dplyr::summarise(Abs_Reads = sum(Counts), .groups="drop")
              base_spec_reads <- spec_reads_sums %>% dplyr::filter(CS_Num == min_cs_val) %>% dplyr::pull(Abs_Reads)
              if(length(base_spec_reads)==0) base_spec_reads <- 1

              curve_spec_reads <- spec_reads_sums %>%
                dplyr::mutate(Percentage = (Abs_Reads / base_spec_reads) * 100, Legend_Label = lbl_assign) %>%
                dplyr::select(CS_Num, Percentage, Legend_Label)

              spec_taxa_counts <- spec_level_data %>% dplyr::group_by(CS_Num) %>% dplyr::summarise(Abs_Taxa = dplyr::n(), .groups="drop")
              base_spec_taxa <- spec_taxa_counts %>% dplyr::filter(CS_Num == min_cs_val) %>% dplyr::pull(Abs_Taxa)

              curve_spec_taxa <- spec_taxa_counts %>%
                dplyr::mutate(Percentage = (Abs_Taxa / base_spec_taxa) * 100, Legend_Label = PLURAL_MAP[ANALYSIS_LEVEL]) %>%
                dplyr::select(CS_Num, Percentage, Legend_Label)

              df_plot <- dplyr::bind_rows(curve_total_reads, curve_spec_reads, curve_spec_taxa)
              fmt_ass <- format(base_spec_reads, big.mark = ".", scientific = FALSE)
              fmt_tax <- format(base_spec_taxa, big.mark = ".", scientific = FALSE)
            } else {
              df_plot <- curve_total_reads; fmt_ass <- "0"; fmt_tax <- "0"
            }

            df_plot$Legend_Label <- factor(df_plot$Legend_Label, levels = specific_levels_order)
            fmt_tot <- format(TOTAL_READS_BASELINE, big.mark = ".", scientific = FALSE)
            subtitle_str <- paste0(fmt_tot, " Total Reads; ", fmt_ass, " ", lbl_assign, " | ", fmt_tax, " ", lbl_plural)

            p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = CS_Num, y = Percentage, color = Legend_Label, group = Legend_Label, shape = Legend_Label, linetype = Legend_Label)) +
              ggplot2::geom_line(linewidth = 0.8) +
              ggplot2::geom_point(size = 2.5) +
              ggplot2::scale_color_manual(values = MODE1_COLORS, drop = FALSE) +
              ggplot2::scale_shape_manual(values = MODE1_SHAPES, drop = FALSE) +
              ggplot2::scale_linetype_manual(values = setNames(c("solid", "dashed", "dashed"), specific_levels_order), drop = FALSE) +
              ggplot2::scale_y_log10(limits = c(0.01, 100), breaks = c(0.01, 0.1, 1, 10, 100), labels = c("0.01", "0.1", "1", "10", "100")) +
              ggplot2::scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
              ggplot2::labs(title = dom, subtitle = subtitle_str, x = "Confidence Score (%)", y = y_axis_label, color = NULL, shape = NULL, linetype = NULL) +
              ggplot2::theme_classic(base_size = 10) +
              ggplot2::guides(color = ggplot2::guide_legend(nrow=1), shape = ggplot2::guide_legend(nrow=1), linetype = ggplot2::guide_legend(nrow=1)) +
              ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 11), plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 8, color = "grey30"), legend.position = "bottom", panel.grid.major.y = ggplot2::element_line(color = "grey90", linetype = "dotted"), plot.margin = ggplot2::margin(t = 10, r = 15, b = 10, l = 10)) +
              ggplot2::coord_cartesian(clip = "off")

            plot_list[[dom]] <- p

          } else {
            # MODE 2: All Levels
            stats_taxa <- taxa_dom %>% dplyr::group_by(CS_Num, Tax_Level) %>% dplyr::summarise(Abs_Taxa = dplyr::n(), .groups="drop")

            if(nrow(stats_taxa) == 0) {
              df_plot <- data.frame(CS_Num = 0, Percentage = 0, Legend_Label = factor("Species", levels = PLURAL_ORDER)) %>% dplyr::filter(CS_Num == -1)
            } else {
              base_taxa <- stats_taxa %>% dplyr::filter(CS_Num == min_cs_val) %>% dplyr::rename(Base_Taxa = Abs_Taxa) %>% dplyr::select(-CS_Num)
              df_plot <- stats_taxa %>% dplyr::left_join(base_taxa, by = "Tax_Level") %>% dplyr::mutate(Percentage = (Abs_Taxa / Base_Taxa) * 100, Legend_Label = PLURAL_MAP[Tax_Level]) %>% tidyr::replace_na(list(Percentage = 0))
            }

            df_plot$Legend_Label <- factor(df_plot$Legend_Label, levels = PLURAL_ORDER)
            fmt_reads <- format(TOTAL_READS_BASELINE, big.mark = ".", scientific = FALSE)
            sub_parts <- c(paste("Reads:", fmt_reads))
            tax_codes <- c("Phyla"="P", "Classes"="C", "Orders"="O", "Families"="F", "Genera"="G", "Species"="S")
            existing_levels_data <- unique(df_plot$Legend_Label)

            if(nrow(stats_taxa) > 0) {
              for(lvl_name in PLURAL_ORDER) {
                if (lvl_name %in% existing_levels_data) {
                  orig_lvl <- names(PLURAL_MAP)[which(PLURAL_MAP == lvl_name)]
                  val <- base_taxa %>% dplyr::filter(Tax_Level == orig_lvl) %>% dplyr::pull(Base_Taxa)
                  if(length(val) > 0) sub_parts <- c(sub_parts, paste0(tax_codes[lvl_name], ": ", val))
                }
              }
            }
            subtitle_str <- paste(sub_parts, collapse = " | ")

            p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = CS_Num, y = Percentage, color = Legend_Label, group = Legend_Label, shape = Legend_Label)) +
              ggplot2::geom_line(linetype = "solid", linewidth = 0.8) + ggplot2::geom_point(size = 2.5) +
              ggplot2::scale_color_manual(values = TAX_COLORS, drop = FALSE) + ggplot2::scale_shape_manual(values = TAX_SHAPES, drop = FALSE) +
              ggplot2::scale_y_log10(limits = c(0.01, 100), breaks = c(0.01, 0.1, 1, 10, 100), labels = c("0.01", "0.1", "1", "10", "100")) +
              ggplot2::scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
              ggplot2::labs(title = dom, subtitle = subtitle_str, x = "Confidence Score (%)", y = y_axis_label, color = NULL, shape = NULL) +
              ggplot2::theme_classic(base_size = 10) + ggplot2::guides(color = ggplot2::guide_legend(nrow=1), shape = ggplot2::guide_legend(nrow=1)) +
              ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 11), plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 7, color = "grey30"), legend.position = "bottom", panel.grid.major.y = ggplot2::element_line(color = "grey90", linetype = "dotted"), plot.margin = ggplot2::margin(t = 10, r = 15, b = 10, l = 10)) +
              ggplot2::coord_cartesian(clip = "off")
            plot_list[[dom]] <- p
          }
        } # End Domain Loop

        # Patchwork Layout
        final_layout <- (plot_list[["Bacteria"]] | plot_list[["Archaea"]]) / (plot_list[["Eukaryota"]] | plot_list[["Viruses"]]) +
          patchwork::plot_layout(guides = "collect") &
          ggplot2::theme(legend.position = "bottom", legend.text = ggplot2::element_text(size = 9, face = "bold")) &
          patchwork::plot_annotation(
            title = paste(samp, "- CS Retention Analysis"),
            subtitle = paste("Retention relative to Baseline (CS00) | Mode:", ifelse(is.null(ANALYSIS_LEVEL), "All Levels", ANALYSIS_LEVEL)),
            theme = ggplot2::theme(plot.title = ggplot2::element_text(face="bold", size=14, hjust=0.5))
          )

        file_name <- paste0(samp, "_CS_Retention_", mode_tag, ".pdf")
        out_path <- file.path(samp_out_dir, file_name)
        ggplot2::ggsave(out_path, final_layout, width = 11.69, height = 8.27, units = "in")

      } # End Mode Loop
    } # End Sample Loop

    log_msg("SUCCESS: All CS Retention plots generated.")

  }, error = function(e) {
    log_msg("CRITICAL ERROR:", e$message)
    stop(e)
  }, finally = {
    sink() # Close log explicitly
  })
}
