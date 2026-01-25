#' Run Confidence Score Retention Analysis (Step 001)
#'
#' Executes read retention analysis based on Confidence Score (Kraken/Bracken).
#'
#' @param project_dir Root path of the project.
#' @param import_script_path Deprecated.
#'
#' @export
#' @importFrom readr read_rds
#' @importFrom dplyr rename mutate filter select group_by summarise n distinct left_join bind_rows pull case_when any_of
#' @importFrom tidyr pivot_longer replace_na
#' @importFrom stringr str_detect str_remove str_extract
#' @importFrom ggplot2 ggplot aes geom_line geom_point scale_color_manual scale_shape_manual scale_linetype_manual scale_y_log10 scale_x_continuous labs guides guide_legend element_line coord_cartesian annotate ggsave theme_classic theme element_text element_blank margin
#' @importFrom patchwork plot_layout plot_annotation wrap_plots
#' @import patchwork

taxa_retention <- function(project_dir, import_script_path = NULL) {

  # ==============================================================================
  # 1. SETUP & DIRECTORIES
  # ==============================================================================
  input_dir   <- file.path(project_dir, "000_mpa_original")
  output_dir  <- file.path(project_dir, "001_cs_retention")

  if (!dir.exists(input_dir)) stop("Input directory not found: ", input_dir)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # ==============================================================================
  # 2. LOAD DATA
  # ==============================================================================
  rds_file <- file.path(input_dir, "combined_mpas_matrix_taxon_level.rds")

  if (!file.exists(rds_file)) {
    rds_old <- file.path(input_dir, "000_unified_MPA_matrix.rds")
    if(file.exists(rds_old)) { rds_file <- rds_old }
    else { stop("Unified RDS file not found.") }
  }

  df_long <- load_and_process_mpa_archive(rds_file)

  # ==============================================================================
  # 3. ANALYSIS LOOP
  # ==============================================================================
  SAMPLES <- unique(df_long$Sample)
  DOMAINS <- c("Bacteria", "Archaea", "Viruses", "Eukaryota")

  file_suffix_map <- c("Phylum"="Phyla", "Class"="Classes", "Order"="Orders",
                       "Family"="Families", "Genus"="Genera", "Species"="Species")
  TAX_LEVELS <- c(NA, "Phylum", "Class", "Order", "Family", "Genus", "Species")

  # --- DEFINIÇÃO DE ESTILOS ---

  # 1. Cores e Shapes para ALL LEVELS
  orig_colors <- c(
    "Phylum"  = "#000000", # Preto
    "Class"   = "#E69F00", # Laranja
    "Order"   = "#56B4E9", # Azul Céu
    "Family"  = "#009E73", # Verde Azulado
    "Genus"   = "#F0E442", # Amarelo
    "Species" = "#D55E00"  # Vermelho
  )
  orig_shapes <- c(
    "Phylum"  = 15, "Class" = 16, "Order" = 17,
    "Family"  = 18, "Genus" = 25, "Species" = 8
  )

  cat("\n=== Starting Retention Analysis (Step 001) ===\n")

  for (samp in SAMPLES) {
    cat("  Processing Sample:", samp, "\n")
    samp_out_dir <- file.path(output_dir, samp)
    if (!dir.exists(samp_out_dir)) dir.create(samp_out_dir)

    # Preparar Dataset
    df_samp <- df_long %>%
      dplyr::filter(Sample == samp) %>%
      dplyr::mutate(
        CS_Num_Raw = as.numeric(stringr::str_remove(CS, "CS")),
        CS_Pct = CS_Num_Raw * 100
      )

    # Baselines (CS0.0)
    df_baseline_domain <- df_samp %>%
      dplyr::filter(CS == "CS0.0", Level == "Domain") %>%
      dplyr::select(Domain, Base_Domain_Reads = Counts)

    df_baseline_levels <- df_samp %>%
      dplyr::filter(CS == "CS0.0") %>%
      dplyr::group_by(Domain, Level) %>%
      dplyr::summarise(
        Base_Lvl_Reads = sum(Counts),
        Base_Lvl_Taxa = dplyr::n_distinct(Taxon[Counts > 0]),
        .groups = "drop"
      )

    # Calcular Métricas Globais
    df_stats <- df_samp %>%
      dplyr::group_by(Domain, CS, CS_Pct, Level) %>%
      dplyr::summarise(
        Curr_Reads = sum(Counts),
        Curr_Taxa = dplyr::n_distinct(Taxon[Counts > 0]),
        .groups = "drop"
      ) %>%
      dplyr::left_join(df_baseline_levels, by = c("Domain", "Level")) %>%
      dplyr::left_join(df_baseline_domain, by = "Domain") %>%
      dplyr::mutate(
        Ret_Lvl_Reads = (Curr_Reads / Base_Lvl_Reads) * 100,
        Ret_Lvl_Taxa  = (Curr_Taxa / Base_Lvl_Taxa) * 100
      ) %>%
      # Ponto 1: Substituir 0 por NA para interromper a linha no gráfico
      dplyr::mutate(
        Ret_Lvl_Reads = ifelse(Ret_Lvl_Reads == 0, NA, Ret_Lvl_Reads),
        Ret_Lvl_Taxa  = ifelse(Ret_Lvl_Taxa == 0, NA, Ret_Lvl_Taxa)
      )

    df_domain_trend <- df_stats %>%
      dplyr::filter(Level == "Domain") %>%
      dplyr::select(Domain, CS, Ret_Total_Reads_Trend = Ret_Lvl_Reads)

    # LOOP DE PLOTAGEM
    for (lvl in TAX_LEVELS) {

      plot_list <- list()

      if (is.na(lvl)) {
        ANALYSIS_LEVEL <- NULL
        mode_tag <- "All_Levels"
      } else {
        ANALYSIS_LEVEL <- lvl
        mode_tag <- file_suffix_map[[lvl]]
      }

      for (dom in DOMAINS) {

        # Subtítulo
        base_dom_val <- df_baseline_domain$Base_Domain_Reads[df_baseline_domain$Domain == dom]
        if(length(base_dom_val)==0) base_dom_val <- 0
        reads_str <- format(base_dom_val, big.mark=".", decimal.mark=",")

        get_cnt <- function(l) {
          val <- df_baseline_levels$Base_Lvl_Taxa[df_baseline_levels$Domain == dom & df_baseline_levels$Level == l]
          if(length(val)==0) 0 else val
        }

        # --- MODE: ALL LEVELS ---
        if (is.null(ANALYSIS_LEVEL)) {

          subtitle_stats <- paste0("Reads: ", reads_str, " | P: ", get_cnt("Phylum"), " | C: ", get_cnt("Class"),
                                   " | O: ", get_cnt("Order"), " | F: ", get_cnt("Family"),
                                   " | G: ", get_cnt("Genus"), " | S: ", get_cnt("Species"))

          dat_plot <- df_stats %>%
            dplyr::filter(Domain == dom, Level %in% names(orig_colors))

          # Ordenação Legenda
          dat_plot$Level <- factor(dat_plot$Level, levels = c("Phylum", "Class", "Order", "Family", "Genus", "Species"))

          # O ggplot remove automaticamente NAs da linha e do ponto
          p <- ggplot2::ggplot(dat_plot, ggplot2::aes(x = CS_Pct, y = Ret_Lvl_Taxa, color = Level, shape = Level)) +
            ggplot2::geom_line(linewidth = 0.8) +
            ggplot2::geom_point(size = 2.5) +
            ggplot2::scale_color_manual(values = orig_colors) +
            ggplot2::scale_shape_manual(values = orig_shapes) +
            ggplot2::guides(color = ggplot2::guide_legend(nrow = 1), shape = ggplot2::guide_legend(nrow = 1)) +
            ggplot2::labs(subtitle = subtitle_stats)

        } else {
          # --- MODE: SPECIFIC LEVEL ---

          curr_reads_base <- df_baseline_levels$Base_Lvl_Reads[df_baseline_levels$Domain == dom & df_baseline_levels$Level == lvl]
          curr_reads_str <- format(curr_reads_base, big.mark=".", decimal.mark=",")

          subtitle_stats <- paste0(reads_str, " Total Reads; ",
                                   curr_reads_str, " assigned to ", lvl, " | ",
                                   get_cnt(lvl), " ", mode_tag)

          dat_lvl <- df_stats %>%
            dplyr::filter(Domain == dom, Level == lvl) %>%
            dplyr::left_join(df_domain_trend, by = c("Domain", "CS")) %>%
            dplyr::select(CS_Pct, Ret_Total_Reads_Trend, Ret_Lvl_Reads, Ret_Lvl_Taxa) %>%
            tidyr::pivot_longer(cols = -CS_Pct, names_to = "Metric", values_to = "Pct")

          # Pontos 2 e 3: Ajuste de Ordem e Labels da Legenda
          # Nova Ordem: 1. Taxa (Nome do Nível), 2. Total Reads, 3. Reads do Nível

          new_levels <- c("Ret_Lvl_Taxa", "Ret_Total_Reads_Trend", "Ret_Lvl_Reads")
          new_labels <- c(lvl, "Total Reads", paste(lvl, "Reads")) # ex: "Species", "Total Reads", "Species Reads"

          dat_lvl$Metric <- factor(dat_lvl$Metric, levels = new_levels, labels = new_labels)

          # Definição dos Estilos baseada na ORDEM dos levels (1, 2, 3)
          # 1. Taxa (Vermelho, Solid, Bola)
          # 2. Total (Preto, Longdash, Quadrado)
          # 3. Reads (Azul, Longdash, Triangulo)

          final_colors <- c("#E64B35", "#000000", "#3C5488")
          final_lines  <- c("solid", "longdash", "longdash")
          final_shapes <- c(16, 15, 17)

          p <- ggplot2::ggplot(dat_lvl, ggplot2::aes(x = CS_Pct, y = Pct, group = Metric)) +
            ggplot2::geom_line(ggplot2::aes(color = Metric, linetype = Metric), linewidth = 1) +
            ggplot2::geom_point(ggplot2::aes(color = Metric, shape = Metric), size = 3) +

            # Aplicação Manual na Ordem do Fator
            ggplot2::scale_color_manual(values = final_colors) +
            ggplot2::scale_linetype_manual(values = final_lines) +
            ggplot2::scale_shape_manual(values = final_shapes) +

            ggplot2::guides(
              color = ggplot2::guide_legend(nrow = 1),
              shape = ggplot2::guide_legend(nrow = 1),
              linetype = ggplot2::guide_legend(nrow = 1)
            ) +
            ggplot2::labs(subtitle = subtitle_stats)
        }

        # --- ESTILO COMUM ---
        p <- p +
          ggplot2::labs(
            title = dom,
            x = "Confidence Score (%)",
            y = "% Retained\n(axis scaled to log10)"
          ) +
          ggplot2::scale_x_continuous(breaks = seq(0, 100, 20), limits = c(0, 100)) +
          ggplot2::scale_y_log10(
            limits = c(0.01, 105),
            breaks = c(0.01, 0.1, 1, 10, 100),
            labels = c("0.01", "0.1", "1", "10", "100")
          ) +
          ggplot2::theme_classic(base_size = 12, base_family = "sans") +
          ggplot2::theme(
            plot.title = ggplot2::element_text(face = "bold", size = 14, hjust = 0.5, color = "black"),
            plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10, color = "grey30"),
            axis.text = ggplot2::element_text(color = "black", size = 10),
            axis.title = ggplot2::element_text(color = "black", size = 11, face = "bold"),
            axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 10)),
            panel.grid.major.y = ggplot2::element_line(color = "grey90", linetype = "dotted"),
            legend.position = "bottom",
            legend.title = ggplot2::element_blank(),
            legend.text = ggplot2::element_text(size = 10),
            legend.background = ggplot2::element_blank(),
            legend.box.background = ggplot2::element_blank(),
            plot.margin = ggplot2::margin(10, 15, 10, 10)
          ) +
          ggplot2::coord_cartesian(clip = "off")

        plot_list[[dom]] <- p
      }

      # Layout Final
      final_layout <- (plot_list[["Bacteria"]] | plot_list[["Archaea"]]) /
        (plot_list[["Eukaryota"]] | plot_list[["Viruses"]]) +
        patchwork::plot_annotation(
          title = paste(samp, "- CS Retention Analysis"),
          subtitle = paste("Retention relative to Baseline (CS00) | Mode:", mode_tag),
          theme = ggplot2::theme(plot.title = ggplot2::element_text(face="bold", size=16, hjust=0.5))
        ) +
        patchwork::plot_layout(guides = "collect") &
        ggplot2::theme(legend.position = "bottom")

      file_name <- paste0(samp, "_CS_Retention_", mode_tag, ".pdf")
      out_path <- file.path(samp_out_dir, file_name)
      ggplot2::ggsave(out_path, final_layout, width = 11.69, height = 8.27, units = "in")
    }
  }

  cat("\nSUCCESS: Retention analysis completed.\n")
  return(invisible(NULL))
}
