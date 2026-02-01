#' Generate Taxa Resolution Analysis (Step 004)
#'
#' Creates stacked bar plots showing the resolution efficiency between a Parent Rank
#' (e.g., Genus) and a Child Rank (e.g., Species).
#' Displays how many reads were fully resolved to the Child level vs. those that
#' remained ambiguous at the Parent level (Exclusive).
#'
#' @param project_dir Path to the project root.
#' @param parent_level Name of the parent rank (default: "Genus").
#' @param child_level Name of the child rank (default: "Species").
#' @param top_n Number of top taxa to display per domain (default: 10).
#'
#' @export
#' @importFrom readr read_rds
#' @importFrom dplyr filter mutate select group_by summarise arrange slice_head left_join bind_rows distinct case_when pull
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_col scale_fill_manual scale_x_continuous scale_y_discrete labs coord_flip theme element_text guides guide_legend
#' @importFrom patchwork plot_layout plot_annotation
#' @importFrom scales label_number
#' @importFrom ggtext element_markdown

taxa_resolution <- function(project_dir,
                            parent_level = "Genus",
                            child_level = "Species",
                            top_n = 10) {

  # ==============================================================================
  # 1. SETUP & DIRECTORIES
  # ==============================================================================
  input_dir  <- file.path(project_dir, "000_karioCaS_input_matrix")
  output_dir <- file.path(project_dir, "004_taxa_resolution")

  if (!dir.exists(input_dir)) stop("CRITICAL ERROR: Input directory not found: ", input_dir)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  rds_file <- file.path(input_dir, "karioCaS_input_matrix.rds")
  if (!file.exists(rds_file)) stop("Input RDS file not found: ", rds_file)

  # ==============================================================================
  # 2. DATA LOADING
  # ==============================================================================
  message(">>> Loading Tidy Data from: ", rds_file)
  df_long <- readr::read_rds(rds_file)

  # Dynamically select columns for Parent and Child names
  df_proc <- df_long %>%
    dplyr::mutate(
      Parent_Name = .data[[parent_level]],
      Child_Name  = .data[[child_level]]
    )

  SAMPLES <- unique(df_proc$sample)
  DOMAINS <- names(kariocas_colors$domains)
  CS_LIST <- unique(df_proc$CS)

  message(">>> Starting Resolution Analysis (", parent_level, " vs ", child_level, ")")

  # ==============================================================================
  # 3. ANALYSIS LOOP
  # ==============================================================================
  for (samp in SAMPLES) {
    message("  Processing Sample: ", samp)

    for (cs in CS_LIST) {

      plot_list <- list()
      df_curr <- df_proc %>% dplyr::filter(sample == samp, CS == cs)

      if (nrow(df_curr) == 0) next

      for (dom in DOMAINS) {

        df_dom <- df_curr %>% dplyr::filter(Domain == dom)

        # 3.1 IDENTIFY TOP PARENTS (Based on Total/Cumulative Reads) -----------
        # Since input is cumulative, the row where Rank == Parent represents the Total.
        parent_totals <- df_dom %>%
          dplyr::filter(Rank == parent_level) %>%
          dplyr::group_by(Parent_Name) %>%
          dplyr::summarise(Total_Clade_Reads = sum(Counts), .groups = "drop") %>%
          dplyr::filter(!is.na(Parent_Name)) %>%
          dplyr::arrange(dplyr::desc(Total_Clade_Reads)) %>%
          dplyr::slice_head(n = top_n)

        top_parents <- parent_totals$Parent_Name

        if (length(top_parents) == 0) {
          plot_list[[dom]] <- plot_kariocas_empty(
            title_text = dom,
            subtitle_text = "No data for this rank",
            x_label = kariocas_labels$x_log10_reads
          )
          next
        }

        # 3.2 CALCULATE EXCLUSIVE VS RESOLVED ----------------------------------
        # Logic:
        # Parent_Exclusive = (Reads at Parent Rank) - (Sum of Reads at Child Rank)

        calc_df <- df_dom %>%
          dplyr::filter(Parent_Name %in% top_parents) %>%
          dplyr::group_by(Parent_Name) %>%
          dplyr::summarise(
            # Reads specifically at the Parent row (Cumulative Total)
            Parent_Total = sum(Counts[Rank == parent_level]),
            # Reads sum of all children rows
            Child_Sum    = sum(Counts[Rank == child_level]),
            .groups = "drop"
          ) %>%
          dplyr::mutate(
            Parent_Exclusive = Parent_Total - Child_Sum,
            # Safety: Ensure no negative numbers if data structure varies
            Parent_Exclusive = pmax(0, Parent_Exclusive)
          )

        # 3.3 RESHAPE FOR PLOTTING ---------------------------------------------
        stack_data <- calc_df %>%
          dplyr::select(Parent_Name, Child_Sum, Parent_Exclusive) %>%
          tidyr::pivot_longer(
            cols = c(Child_Sum, Parent_Exclusive),
            names_to = "Category_Raw",
            values_to = "Reads"
          ) %>%
          dplyr::mutate(Category = dplyr::case_when(
            Category_Raw == "Child_Sum" ~ "Species",
            Category_Raw == "Parent_Exclusive" ~ "Genus-exclusive"
          ))

        # 3.4 GHOST BARS LOGIC (Consistent Thickness) --------------------------
        missing_slots <- top_n - length(top_parents)
        final_levels  <- top_parents

        if (missing_slots > 0) {
          ghost_names <- paste0("Ghost_", seq_len(missing_slots))

          # Add ghosts with 0 reads, assigned to "Genus-exclusive" to avoid new legend key
          ghost_df <- data.frame(
            Parent_Name = ghost_names,
            Category = "Genus-exclusive",
            Reads = 0
          )
          stack_data <- dplyr::bind_rows(stack_data, ghost_df)
          final_levels <- c(final_levels, ghost_names)
        }

        # Enforce Factors
        stack_data$Parent_Name <- factor(stack_data$Parent_Name, levels = rev(final_levels))
        stack_data$Category    <- factor(stack_data$Category, levels = c("Genus-exclusive", "Species"))

        # Map to Style Colors
        # Using Style keys: "Parent" (for exclusive) and "Child" (for species)
        fill_values <- c(
          "Genus-exclusive" = kariocas_colors$special[["Parent"]],
          "Species"         = kariocas_colors$special[["Child"]]
        )

        # 3.5 PLOT CONSTRUCTION ------------------------------------------------
        p <- ggplot2::ggplot(stack_data, ggplot2::aes(x = Reads, y = Parent_Name, fill = Category)) +
          ggplot2::geom_col(width = 0.7) +

          ggplot2::scale_fill_manual(values = fill_values) +

          ggplot2::scale_x_continuous(labels = label_k_number) +

          ggplot2::labs(
            title = dom,
            x = "Reads",
            y = NULL
          ) +

          theme_kariocas() +

          ggplot2::theme(
            axis.text.y = ggplot2::element_text(face = "italic"),
            panel.grid.major.y = ggplot2::element_blank(),
            legend.position = "bottom",
            # Ensure legend title is gone
            legend.title = ggplot2::element_blank()
          )

        # Remove "Ghost" labels from Y axis
        clean_y_labels <- function(x) { ifelse(grepl("^Ghost_", x), "", x) }
        p <- p + ggplot2::scale_y_discrete(labels = clean_y_labels)

        plot_list[[dom]] <- p

      } # End Domain Loop

      # 3.6 ASSEMBLE PANEL -----------------------------------------------------
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

    } # End CS Loop
  } # End Sample Loop

  message("\nSUCCESS: Resolution analysis completed. Files saved in: ", output_dir)
  return(invisible(NULL))
}
