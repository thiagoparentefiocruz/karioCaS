#' Generate Taxa Resolution Plots (Step 004)
#'
#' Generates publication-quality (Nature-style) stacked bar plots showing the
#' proportion of reads assigned to a child level (e.g., Species) vs. the
#' exclusive reads of its parent level (e.g., Genus).
#'
#' @param project_dir Path to the project root (e.g., "/Users/name/project_x").
#' @param analysis_level The taxonomic level to analyze as the "Child".
#'        Options: "Species" (default, compares to Genus) or "Genus" (compares to Family).
#' @param import_script_path Deprecated. The package uses internal functions.
#' @param max_bars Integer. Maximum number of taxa to show per plot (Default: 10).
#'
#' @export
#' @importFrom readr read_rds
#' @importFrom dplyr rename_with filter mutate group_by summarise arrange left_join select bind_rows desc case_when
#' @importFrom tidyr pivot_longer replace_na
#' @importFrom stringr str_detect str_remove str_extract
#' @importFrom ggplot2 ggplot aes geom_col scale_fill_manual scale_x_continuous scale_y_discrete labs theme_classic theme element_text element_line element_blank expansion ggsave
#' @importFrom scales label_number cut_short_scale
#' @importFrom patchwork plot_layout plot_annotation
#' @importFrom tidyselect any_of all_of
#' @importFrom stats setNames
#' @importFrom utils head

taxa_resolution <- function(project_dir,
                            analysis_level = "Species",
                            import_script_path = NULL,
                            max_bars = 10) {

  # ==============================================================================
  # 1. SETUP & DIRECTORIES
  # ==============================================================================

  input_dir   <- file.path(project_dir, "000_mpa_original")
  output_dir  <- file.path(project_dir, "004_taxa_resolution_per_ConfidenceScore")
  log_dir     <- file.path(project_dir, "log_archives")
  log_file    <- file.path(log_dir, "log_004_taxa_resolution.txt")

  if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  if(!dir.exists(log_dir))    dir.create(log_dir, recursive = TRUE)

  # ==============================================================================
  # 2. TAXONOMIC MAPPING LOGIC
  # ==============================================================================

  tax_map <- list(
    "Species" = list(code = "s", parent = "Genus",  parent_code = "g"),
    "Genus"   = list(code = "g", parent = "Family", parent_code = "f"),
    "Family"  = list(code = "f", parent = "Order",  parent_code = "o"),
    "Order"   = list(code = "o", parent = "Class",  parent_code = "c"),
    "Class"   = list(code = "c", parent = "Phylum", parent_code = "p"),
    "Phylum"  = list(code = "p", parent = "Kingdom",parent_code = "k")
  )

  if (is.null(analysis_level) || analysis_level == "") analysis_level <- "Species"
  if (is.null(tax_map[[analysis_level]])) stop("Taxonomic level not supported: ", analysis_level)

  LVL_CHILD_NAME  <- analysis_level
  LVL_PARENT_NAME <- tax_map[[analysis_level]]$parent
  LVL_PARENT_CODE <- tax_map[[analysis_level]]$parent_code

  # Labels for Legend
  LBL_EXCLUSIVE <- paste0(LVL_PARENT_NAME, "-exclusive")
  LBL_CHILD     <- LVL_CHILD_NAME

  # Dynamic Color Palette (Nature Style: Blue vs Grey)
  COLOR_PALETTE <- stats::setNames(c("#3C5488", "#B0B0B0"), c(LBL_CHILD, LBL_EXCLUSIVE))

  # ==============================================================================
  # 3. START LOGGING
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
  cat("LOG: 004_TAXA_RESOLUTION (taxa_resolution)\n")
  cat("DATE:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("PROJECT DIR:", project_dir, "\n")
  cat("ANALYSIS:", LVL_CHILD_NAME, "vs", LVL_PARENT_NAME, "\n")
  cat("====================================================\n\n")

  tryCatch({

    # ==============================================================================
    # 4. INPUT INTELIGENTE (STANDARDIZED)
    # ==============================================================================
    cat("STEP 1: Checking Input Data...\n")

    file_rds <- file.path(input_dir, "combined_mpas_matrix_taxon_level.rds")
    INPUT_MATRIX <- NULL

    # 4.1. Check RDS
    if (file.exists(file_rds)) {
      cat("  -> Input file found:", file_rds, "\n")
      cat("  -> Loading RDS object...\n")
      INPUT_MATRIX <- readr::read_rds(file_rds)

    } else {
      cat("  -> [WARNING] Input RDS NOT found.\n")
      cat("  -> Attempting to run Import Function (Step 000)...\n")

      # PACKAGE ADAPTATION: Call internal function directly
      tryCatch({
        INPUT_MATRIX <- import_mpa_data(project_dir = project_dir)
      }, error = function(e) {
        stop("CRITICAL: Input RDS missing and could not run internal 'import_mpa_data'.\nError: ", e$message)
      })
    }

    # 4.3. Validate Matrix
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
    # 5. DATA PROCESSING (OPTIMIZED)
    # ==============================================================================
    cat("STEP 2: Processing Data...\n")

    regex_parent <- paste0(LVL_PARENT_CODE, "__[^|]+")

    # Initial Pivot & Metadata
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
      )

    # --- CALCULATION LOGIC ---

    # 1. Process CHILDREN (e.g., Species)
    # Extract Parent name from Child string
    df_child_sum <- df_long %>%
      dplyr::filter(Tax_Level == LVL_CHILD_NAME) %>%
      dplyr::mutate(
        Parent_Taxon_Raw = stringr::str_extract(Taxonomy, regex_parent),
        Parent_Name      = stringr::str_remove(Parent_Taxon_Raw, paste0("^", LVL_PARENT_CODE, "__"))
      ) %>%
      dplyr::group_by(Sample_Name, Domain, CS_Label, Parent_Name) %>%
      dplyr::summarise(Sum_Child_Counts = sum(Counts), .groups="drop")

    # 2. Process PARENTS (e.g., Genus)
    # Extract Parent name directly from Parent string
    df_parent_rows <- df_long %>%
      dplyr::filter(Tax_Level == LVL_PARENT_NAME) %>%
      dplyr::mutate(
        Parent_Name = stringr::str_remove(stringr::str_extract(Taxonomy, "[^|]+$"), "^[kpcofgs]__")
      ) %>%
      dplyr::rename(Total_Parent_Count = Counts)

    # 3. Join & Subtraction
    df_calc <- df_parent_rows %>%
      dplyr::left_join(df_child_sum, by = c("Sample_Name", "Domain", "CS_Label", "Parent_Name")) %>%
      dplyr::mutate(Sum_Child_Counts = tidyr::replace_na(Sum_Child_Counts, 0)) %>%
      dplyr::mutate(
        Parent_Exclusive_Reads = ifelse(Total_Parent_Count >= Sum_Child_Counts,
                                        Total_Parent_Count - Sum_Child_Counts, 0)
      )

    # 4. Format for Plotting
    df_plot_ready <- df_calc %>%
      dplyr::select(Sample_Name, Domain, CS_Label, Parent_Name,
                    !!LBL_CHILD     := Sum_Child_Counts,
                    !!LBL_EXCLUSIVE := Parent_Exclusive_Reads) %>%
      tidyr::pivot_longer(cols = c(tidyselect::all_of(LBL_CHILD), tidyselect::all_of(LBL_EXCLUSIVE)),
                          names_to = "Assignment_Type",
                          values_to = "Final_Counts")

    unique_samples <- unique(df_plot_ready$Sample_Name)
    cat("  -> Processed Samples:", length(unique_samples), "\n\n")

    # ==============================================================================
    # 6. PLOT GENERATION (CENTERED & PUBLICATION READY)
    # ==============================================================================
    cat("STEP 3: Generating Figures...\n")

    DOMAINS_ORDER <- c("Bacteria", "Archaea", "Eukaryota", "Viruses")

    for (samp in unique_samples) {
      cat("  Processing Sample:", samp, "...\n")

      samp_out_dir <- file.path(output_dir, samp)
      if(!dir.exists(samp_out_dir)) dir.create(samp_out_dir)

      df_samp <- df_plot_ready %>% dplyr::filter(Sample_Name == samp)
      unique_cs_list <- sort(unique(df_samp$CS_Label))

      for (cs in unique_cs_list) {

        df_cs <- df_samp %>% dplyr::filter(CS_Label == cs)
        plot_list <- list()

        for (dom in DOMAINS_ORDER) {

          df_dom <- df_cs %>% dplyr::filter(Domain == dom)

          # Empty Plot Placeholder
          if (nrow(df_dom) == 0 || sum(df_dom$Final_Counts) == 0) {
            plot_list[[dom]] <- ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::labs(title = dom)
            next
          }

          # --- RANKING & SELECTION ---
          rank_data <- df_dom %>%
            dplyr::group_by(Parent_Name) %>%
            dplyr::summarise(Total = sum(Final_Counts), .groups="drop") %>%
            dplyr::arrange(dplyr::desc(Total)) %>%
            utils::head(max_bars)

          top_taxa <- rank_data$Parent_Name
          df_final <- df_dom %>% dplyr::filter(Parent_Name %in% top_taxa)

          # --- CENTERING LOGIC (DUMMY BARS) ---
          n_current <- length(top_taxa)

          if (n_current < max_bars) {
            total_missing <- max_bars - n_current
            pad_top    <- floor(total_missing / 2)
            pad_bottom <- ceiling(total_missing / 2)

            dummy_top    <- if(pad_top > 0) paste0("DUMMY_TOP_", 1:pad_top) else c()
            dummy_bottom <- if(pad_bottom > 0) paste0("DUMMY_BOT_", 1:pad_bottom) else c()

            all_dummies <- c(dummy_top, dummy_bottom)

            dummy_df <- expand.grid(
              Sample_Name = samp, Domain = dom, CS_Label = cs,
              Parent_Name = all_dummies,
              Assignment_Type = c(LBL_CHILD, LBL_EXCLUSIVE),
              Final_Counts = 0
            )

            df_final <- dplyr::bind_rows(df_final, dummy_df)
            final_levels <- c(dummy_bottom, rev(top_taxa), dummy_top)

          } else {
            final_levels <- rev(top_taxa)
          }

          df_final$Parent_Name <- factor(df_final$Parent_Name, levels = final_levels)
          df_final$Assignment_Type <- factor(df_final$Assignment_Type,
                                             levels = c(LBL_EXCLUSIVE, LBL_CHILD))

          # --- PLOTTING ---
          # Pct Labels (Logic kept from source but not plotting text as it was commented out in input)
          df_labels <- df_final %>%
            dplyr::group_by(Parent_Name) %>%
            dplyr::mutate(Total_Taxon = sum(Final_Counts)) %>%
            dplyr::filter(Assignment_Type == LBL_CHILD, Total_Taxon > 0) %>%
            dplyr::mutate(Pct = Final_Counts / Total_Taxon) %>%
            dplyr::filter(!stringr::str_detect(Parent_Name, "DUMMY"))

          p <- ggplot2::ggplot(df_final, ggplot2::aes(x = Final_Counts, y = Parent_Name, fill = Assignment_Type)) +
            ggplot2::geom_col(position = "stack", width = 0.75) +

            # Smart Labels (Inside Bar) - Commented out as in source
            # ggplot2::geom_text(
            #   data = df_labels %>% dplyr::filter(Pct > 0.10),
            #   aes(label = scales::percent(Pct, accuracy = 1)),
            #   position = position_stack(vjust = 0.5),
            #   color = "white", fontface = "bold", size = 2.5
            # ) +

            ggplot2::scale_fill_manual(values = COLOR_PALETTE) +
            ggplot2::scale_x_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()),
                                        expand = ggplot2::expansion(mult = c(0, 0.1))) +

            # Hide Dummies on Y-axis
            ggplot2::scale_y_discrete(breaks = top_taxa) +

            ggplot2::labs(title = dom, x = "Reads", y = NULL, fill = NULL) +

            # THEME: Nature Style
            ggplot2::theme_classic(base_size = 10) +
            ggplot2::theme(
              plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 11),
              axis.text.y = ggplot2::element_text(color = "black", face = "italic", size = 9),
              axis.text.x = ggplot2::element_text(color = "black", size = 8),
              axis.line   = ggplot2::element_line(linewidth = 0.3, color = "black"),
              legend.position = "none",
              panel.grid = ggplot2::element_blank()
            )

          plot_list[[dom]] <- p
        }

        # --- FINAL LAYOUT ---
        final_layout <- (plot_list[["Bacteria"]] | plot_list[["Archaea"]]) /
          (plot_list[["Eukaryota"]] | plot_list[["Viruses"]]) +
          patchwork::plot_layout(guides = "collect") &
          ggplot2::theme(legend.position = "bottom", legend.text = ggplot2::element_text(size = 9)) &
          patchwork::plot_annotation(
            title = paste(samp, "-", cs),
            subtitle = paste0("Taxa resolution: Top ", max_bars, " ", LVL_PARENT_NAME, " per Domain"),
            theme = ggplot2::theme(plot.title = ggplot2::element_text(face="bold", size=14, hjust=0.5))
          )

        # Save
        file_name <- paste0(samp, "_", cs, "_resolution_", LVL_PARENT_NAME, "_vs_", LVL_CHILD_NAME, ".pdf")
        out_file <- file.path(samp_out_dir, file_name)

        ggplot2::ggsave(out_file, final_layout, width = 11.69, height = 8.27, units = "in")
        cat("    -> Saved:", file_name, "\n")
      }
    }

    cat("\n====================================================\n")
    cat("SUCCESS: All Taxa resolution figures generated.\n")
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
