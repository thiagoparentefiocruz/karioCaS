#' Generate Taxa Resolution Plots (Step 004)
#'
#' Generates stacked bar plots showing child vs parent level reads.
#'
#' @param project_dir Path to the project root.
#' @param analysis_level "Species" (default) or "Genus".
#' @param import_script_path Deprecated.
#' @param max_bars Max taxa per plot.
#'
#' @export
#' @importFrom readr read_rds
#' @importFrom dplyr rename_with filter mutate group_by summarise arrange left_join select bind_rows desc case_when
#' @importFrom tidyr pivot_longer replace_na
#' @importFrom ggplot2 ggplot aes geom_col scale_fill_manual scale_x_continuous scale_y_discrete labs theme expansion ggsave
#' @importFrom scales label_number cut_short_scale
#' @importFrom patchwork plot_layout plot_annotation
#' @importFrom utils head

taxa_resolution <- function(project_dir,
                            analysis_level = "Species",
                            import_script_path = NULL,
                            max_bars = 10) {

  # ==============================================================================
  # 1. SETUP
  # ==============================================================================
  input_dir   <- file.path(project_dir, "000_mpa_original")
  output_dir  <- file.path(project_dir, "004_taxa_resolution_per_ConfidenceScore")

  if (!dir.exists(input_dir)) stop("Input directory not found.")
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  rds_file <- file.path(input_dir, "combined_mpas_matrix_taxon_level.rds")
  df_long <- load_and_process_mpa(rds_file)

  # Determine Levels
  if (analysis_level == "Species") {
    LVL_CHILD <- "Species"; LVL_PARENT <- "Genus"; LVL_CHILD_NAME <- "Species"; LVL_PARENT_NAME <- "Genus"
  } else if (analysis_level == "Genus") {
    LVL_CHILD <- "Genus"; LVL_PARENT <- "Family"; LVL_CHILD_NAME <- "Genus"; LVL_PARENT_NAME <- "Family"
  } else {
    stop("Invalid analysis_level. Use 'Species' or 'Genus'.")
  }

  SAMPLES <- unique(df_long$Sample)
  CS_LIST <- unique(df_long$CS)
  DOMAINS <- c("Bacteria", "Archaea", "Viruses", "Eukaryota")

  cat("\n=== Starting Taxa Resolution (Step 004) ===\n")

  for (samp in SAMPLES) {
    cat("  Processing:", samp, "\n")
    samp_out_dir <- file.path(output_dir, samp)
    if (!dir.exists(samp_out_dir)) dir.create(samp_out_dir)

    df_samp <- df_long %>% dplyr::filter(Sample == samp)

    for (cs in CS_LIST) {
      plot_list <- list()

      for (dom in DOMAINS) {
        # Data Prep (Condensed for brevity, logic identical to original)
        df_child <- df_samp %>% dplyr::filter(CS==cs, Domain==dom, Level==LVL_CHILD)
        df_parent <- df_samp %>% dplyr::filter(CS==cs, Domain==dom, Level==LVL_PARENT)

        if(nrow(df_parent) == 0) {
          plot_list[[dom]] <- patchwork::wrap_plots(ggplot2::ggplot() + ggplot2::theme_void())
          next
        }

        # Parent Selection (Top N)
        top_parents <- df_parent %>%
          dplyr::arrange(dplyr::desc(Counts)) %>%
          utils::head(max_bars) %>%
          dplyr::pull(Taxon)

        # Matched Data construction
        plot_data <- df_parent %>%
          dplyr::filter(Taxon %in% top_parents) %>%
          dplyr::rename(Parent_Taxon = Taxon, Parent_Total = Counts) %>%
          dplyr::left_join(
            df_child %>%
              dplyr::mutate(Parent_Taxon = stringr::str_extract(Taxon, "^\\w+")) %>% # Simple regex for Genus extraction
              dplyr::group_by(Parent_Taxon) %>%
              dplyr::summarise(Child_Sum = sum(Counts), .groups="drop"),
            by = "Parent_Taxon"
          ) %>%
          tidyr::replace_na(list(Child_Sum = 0)) %>%
          dplyr::mutate(Parent_Exclusive = Parent_Total - Child_Sum) %>%
          dplyr::mutate(Parent_Exclusive = ifelse(Parent_Exclusive < 0, 0, Parent_Exclusive)) %>%
          dplyr::select(Parent_Taxon, Child_Sum, Parent_Exclusive) %>%
          tidyr::pivot_longer(cols = c("Child_Sum", "Parent_Exclusive"), names_to = "Category", values_to = "Reads") %>%
          dplyr::mutate(Label = ifelse(Category == "Child_Sum", LVL_CHILD_NAME, paste(LVL_PARENT_NAME, "(Exclusive)")))

        # Factor reordering
        plot_data$Category <- factor(plot_data$Category, levels = c("Parent_Exclusive", "Child_Sum"))
        plot_data$Parent_Taxon <- factor(plot_data$Parent_Taxon, levels = rev(top_parents))

        # --- PLOTTING (Refactored) ---
        p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Reads, y = Parent_Taxon, fill = Label)) +
          ggplot2::geom_col(width = 0.7) +
          # Use Centralized Resolution Colors
          ggplot2::scale_fill_manual(values = setNames(karioCaS_cols$resolution, c(LVL_CHILD_NAME, paste(LVL_PARENT_NAME, "(Exclusive)")))) +
          ggplot2::scale_x_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
          ggplot2::labs(title = dom, x = "Reads", y = NULL) +
          theme_karioCaS() # Central Theme

        plot_list[[dom]] <- p
      }

      final_layout <- (plot_list[["Bacteria"]] | plot_list[["Archaea"]]) /
        (plot_list[["Eukaryota"]] | plot_list[["Viruses"]]) +
        patchwork::plot_annotation(
          title = paste(samp, "-", cs),
          subtitle = paste0("Taxa resolution: ", LVL_PARENT_NAME, " vs ", LVL_CHILD_NAME),
          theme = ggplot2::theme(plot.title = ggplot2::element_text(face="bold", size=14, hjust=0.5))
        )

      file_name <- paste0(samp, "_", cs, "_resolution_", LVL_PARENT_NAME, "_vs_", LVL_CHILD_NAME, ".pdf")
      out_file <- file.path(samp_out_dir, file_name)
      ggplot2::ggsave(out_file, final_layout, width = 11.69, height = 8.27, units = "in")
    }
  }
  cat("\nSUCCESS: Resolution figures generated.\n")
  return(TRUE)
}
