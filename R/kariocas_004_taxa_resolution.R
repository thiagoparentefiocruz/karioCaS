#' Generate Taxa Resolution Analysis (Step 004)
#'
#' Creates stacked bar plots showing the resolution efficiency between a Parent Rank
#' (e.g., Genus) and a Child Rank (e.g., Species).
#' Performs strict subtraction: (Parent Row Value) - (Sum of Child Row Values).
#' Uses max() for Parent aggregation to prevent double-counting of children in cumulative data.
#'
#' @param project_dir Path to the project root.
#' @param parent_level Name of the parent rank (default: "Genus").
#' @param child_level Name of the child rank (default: "Species").
#' @param top_n Number of top taxa to display per domain (default: 10).
#'
#' @export
#' @importFrom dplyr filter mutate select group_by summarise arrange slice_head left_join bind_rows distinct case_when pull rename
#' @importFrom tidyr pivot_longer replace_na
#' @importFrom ggplot2 ggplot aes geom_col scale_fill_manual scale_x_continuous scale_y_discrete labs coord_flip theme element_text guides guide_legend
#' @importFrom patchwork plot_layout plot_annotation
#' @importFrom scales label_number
#' @importFrom ggtext element_markdown
#' @importFrom SummarizedExperiment rowData

taxa_resolution <- function(project_dir,
                            parent_level = "Genus",
                            child_level = "Species",
                            top_n = 10) {

  # ==============================================================================
  # 1. SETUP & LOGGING
  # ==============================================================================
  output_dir <- file.path(project_dir, "004_taxa_resolution")
  log_dir    <- file.path(project_dir, "logs")

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)

  log_file <- file.path(log_dir, "log_004_taxa_resolution.txt")

  log_msg <- function(...) {
    msg <- paste0(...)
    message(msg)
    cat(paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", msg, "\n"),
        file = log_file, append = TRUE)
  }

  cat("====================================================\n", file = log_file)
  cat("LOG: 004_TAXA_RESOLUTION (Max-Based + Doc Fix)\n", file = log_file, append = TRUE)
  cat("PROJECT DIR: ", project_dir, "\n", file = log_file, append = TRUE)
  cat("ANALYSIS: ", parent_level, " vs ", child_level, "\n", file = log_file, append = TRUE)
  cat("====================================================\n", file = log_file, append = TRUE)

  # ==============================================================================
  # 2. DATA LOADING & ENRICHMENT
  # ==============================================================================
  log_msg(">>> Loading Data...")
  df_long <- .get_tidy_data(project_dir)

  # 2.1 Enrich Data if needed
  required_cols <- c(parent_level, child_level)
  if (!all(required_cols %in% colnames(df_long))) {
    log_msg(">>> Enriching data with full taxonomy from TSE...")
    tse_path <- file.path(project_dir, "000_karioCaS_input_matrix", "karioCaS_TSE.rds")
    if (!file.exists(tse_path)) stop("TSE file missing for enrichment.")
    tse <- readRDS(tse_path)
    tax_df <- SummarizedExperiment::rowData(tse) %>% as.data.frame()

    key_col <- "Taxonomy_Full"
    if (!key_col %in% colnames(df_long)) {
      if ("Taxonomy" %in% colnames(df_long)) df_long <- df_long %>% dplyr::rename(Taxonomy_Full = Taxonomy)
      else stop("Could not find Taxonomy key column.")
    }
    cols_to_add <- setdiff(colnames(tax_df), colnames(df_long))
    tax_subset <- tax_df %>% dplyr::select(Taxonomy_Full, dplyr::all_of(cols_to_add))
    df_long <- df_long %>% dplyr::left_join(tax_subset, by = "Taxonomy_Full")
  }

  # 2.3 Create Generic Columns
  df_proc <- df_long %>%
    dplyr::mutate(
      Parent_Name = .data[[parent_level]],
      Child_Name  = .data[[child_level]],
      # Sanitize Rank to avoid matching errors
      Rank = trimws(Rank)
    )

  SAMPLES <- unique(df_proc$sample)
  DOMAINS <- names(kariocas_colors$domains)
  CS_LIST <- unique(df_proc$CS)

  log_msg(">>> Starting Analysis Loop...")

  # ==============================================================================
  # 3. ANALYSIS LOOP
  # ==============================================================================
  for (samp in SAMPLES) {
    log_msg("------------------------------------------------")
    log_msg("  Processing Sample: ", samp)

    for (cs in CS_LIST) {

      plot_list <- list()
      df_curr <- df_proc %>% dplyr::filter(sample == samp, CS == cs)

      if (nrow(df_curr) == 0) next

      for (dom in DOMAINS) {

        df_dom <- df_curr %>% dplyr::filter(Domain == dom)

        # 3.1 Identify Top Parents (Use max to avoid double counting children)
        parent_totals <- df_dom %>%
          dplyr::filter(Rank == parent_level) %>%
          dplyr::group_by(Parent_Name) %>%
          dplyr::summarise(Total_Clade_Reads = max(Counts, na.rm = TRUE), .groups = "drop") %>%
          dplyr::filter(!is.na(Parent_Name)) %>%
          dplyr::arrange(dplyr::desc(Total_Clade_Reads)) %>%
          dplyr::slice_head(n = top_n)

        top_parents <- parent_totals$Parent_Name

        if (length(top_parents) == 0) {
          plot_list[[dom]] <- plot_kariocas_empty(dom, "No Data")
          next
        }

        # ======================================================================
        # 3.2 CALCULATION - THE FIX
        # ======================================================================

        # A: Get Parent Total using MAX (Corrects the 1.19M error)
        df_parent_stats <- df_dom %>%
          dplyr::filter(Parent_Name %in% top_parents, Rank == parent_level) %>%
          dplyr::group_by(Parent_Name) %>%
          dplyr::summarise(Parent_Cumulative_Total = max(Counts, na.rm = TRUE), .groups = "drop")

        # B: Get Children Sum (Standard Sum)
        df_child_stats <- df_dom %>%
          dplyr::filter(Parent_Name %in% top_parents, Rank == child_level) %>%
          dplyr::group_by(Parent_Name) %>%
          dplyr::summarise(Child_Sum_Resolved = sum(Counts, na.rm = TRUE), .groups = "drop")

        # C: Join and Subtract
        calc_df <- df_parent_stats %>%
          dplyr::left_join(df_child_stats, by = "Parent_Name") %>%
          dplyr::mutate(
            Child_Sum_Resolved = tidyr::replace_na(Child_Sum_Resolved, 0),

            # Exclusive = Max(Genus) - Sum(Species)
            Parent_Exclusive = Parent_Cumulative_Total - Child_Sum_Resolved,
            Parent_Exclusive = pmax(0, Parent_Exclusive)
          ) %>%
          dplyr::arrange(dplyr::desc(Parent_Cumulative_Total))

        # --- AUDIT LOG ---
        audit_row <- calc_df %>% dplyr::filter(Parent_Name == "Limnohabitans")
        if (nrow(audit_row) > 0) {
          log_msg(sprintf("    AUDIT [Limnohabitans]: Genus(Max)=%d | Species(Sum)=%d | Exclusive(Diff)=%d",
                          audit_row$Parent_Cumulative_Total, audit_row$Child_Sum_Resolved, audit_row$Parent_Exclusive))
        } else {
          top_1 <- calc_df[1, ]
          log_msg(sprintf("    AUDIT [Top1 %s]: Genus(Max)=%d | Species(Sum)=%d | Exclusive(Diff)=%d",
                          top_1$Parent_Name, top_1$Parent_Cumulative_Total, top_1$Child_Sum_Resolved, top_1$Parent_Exclusive))
        }

        # 3.3 RESHAPE
        label_exclusive <- paste0(parent_level, "-exclusive")
        label_child     <- child_level

        stack_data <- calc_df %>%
          dplyr::select(Parent_Name, Child_Sum_Resolved, Parent_Exclusive) %>%
          tidyr::pivot_longer(
            cols = c(Child_Sum_Resolved, Parent_Exclusive),
            names_to = "Category_Raw",
            values_to = "Reads"
          ) %>%
          dplyr::mutate(Category = dplyr::case_when(
            Category_Raw == "Child_Sum_Resolved" ~ label_child,
            Category_Raw == "Parent_Exclusive" ~ label_exclusive
          ))

        # 3.4 GHOST BARS
        missing_slots <- top_n - length(top_parents)
        final_levels  <- top_parents

        if (missing_slots > 0) {
          ghost_names <- paste0("Ghost_", seq_len(missing_slots))
          ghost_df <- data.frame(Parent_Name = ghost_names, Category = label_exclusive, Reads = 0)
          stack_data <- dplyr::bind_rows(stack_data, ghost_df)
          final_levels <- c(final_levels, ghost_names)
        }

        stack_data$Parent_Name <- factor(stack_data$Parent_Name, levels = rev(final_levels))
        stack_data$Category    <- factor(stack_data$Category, levels = c(label_exclusive, label_child))

        # 3.5 PLOT
        spec_colors <- get_kariocas_colors("special")
        fill_values <- setNames(
          c(spec_colors[["Parent"]], spec_colors[["Child"]]),
          c(label_exclusive, label_child)
        )

        p <- ggplot2::ggplot(stack_data, ggplot2::aes(x = Reads, y = Parent_Name, fill = Category)) +
          ggplot2::geom_col(width = 0.7) +
          ggplot2::scale_fill_manual(values = fill_values) +
          ggplot2::scale_x_continuous(labels = label_k_number) +
          ggplot2::labs(title = dom, x = "Reads", y = NULL) +
          theme_kariocas() +
          ggplot2::theme(
            axis.text.y = ggplot2::element_text(face = "italic"),
            panel.grid.major.y = ggplot2::element_blank(),
            legend.position = "bottom",
            legend.title = ggplot2::element_blank()
          )

        clean_y_labels <- function(x) { ifelse(grepl("^Ghost_", x), "", x) }
        p <- p + ggplot2::scale_y_discrete(labels = clean_y_labels)

        plot_list[[dom]] <- p
      }

      # 3.6 ASSEMBLE
      final_layout <- (plot_list[["Bacteria"]] | plot_list[["Archaea"]]) /
        (plot_list[["Eukaryota"]] | plot_list[["Viruses"]]) +
        patchwork::plot_annotation(
          title = paste(samp, "- CS", sprintf("%02d", cs)),
          subtitle = paste("Taxa Resolution:", parent_level, "vs", child_level),
          theme = ggplot2::theme(
            plot.title = ggplot2::element_text(face = "bold", size = 16, hjust = 0.5),
            plot.subtitle = ggplot2::element_text(size = 12, hjust = 0.5, color = "grey30")
          )
        ) +
        patchwork::plot_layout(guides = "collect") &
        ggplot2::theme(legend.position = "bottom")

      file_name <- paste0(samp, "_CS", sprintf("%02d", cs), "_Resolution_", parent_level, "_vs_", child_level, ".pdf")
      save_path <- file.path(output_dir, file_name)

      ggplot2::ggsave(save_path, final_layout, width = kariocas_dims$width, height = kariocas_dims$height)
      log_msg("    -> Generated: ", file_name)
    }
  }

  log_msg("SUCCESS: Resolution analysis completed.")
  return(invisible(NULL))
}
