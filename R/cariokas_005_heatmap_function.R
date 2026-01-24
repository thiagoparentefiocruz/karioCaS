#' Generate Taxonomic Heatmaps (Step 005)
#'
#' This function processes the unified matrix, handles missing inputs by triggering
#' the internal import step (Step 000) if necessary, and generates detailed heatmaps.
#'
#' @param project_dir Path to the project root (e.g., "/Users/name/project_x").
#' @param import_script_path Deprecated. The package uses internal functions for import.
#' @param confidence_score Numeric threshold (0 to 1). If NULL, auto-detected from the last column.
#' @param top_n Integer. Number of top taxa to display. If NULL, defaults to 20.
#'
#' @return Returns TRUE invisibly if successful. Generates PDFs in "005_taxonomic_heatmaps_CS...".
#' @export
#' @importFrom readr read_rds
#' @importFrom dplyr rename_with filter mutate group_by summarise arrange left_join bind_rows pull desc n
#' @importFrom tidyr pivot_longer replace_na
#' @importFrom stringr str_detect str_replace_all str_extract str_remove
#' @importFrom ggplot2 ggplot aes geom_tile facet_grid scale_fill_gradientn scale_y_discrete labs theme_minimal theme element_text element_rect element_blank ggsave
#' @importFrom scales rescale label_number
#' @importFrom grid unit
#' @importFrom tidyselect any_of
#' @importFrom stats setNames

heatmaps_karioCaS <- function(project_dir,
                                        import_script_path = NULL,
                                        confidence_score = NULL,
                                        top_n = NULL) {

  # ==============================================================================
  # 1. SETUP & DIRECTORIES
  # ==============================================================================

  input_dir   <- file.path(project_dir, "000_mpa_original")
  log_dir     <- file.path(project_dir, "log_archives")
  log_file    <- file.path(log_dir, "log_005_heatmaps.txt")

  if(!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)

  # ==============================================================================
  # 2. INTERNAL HELPER FUNCTIONS
  # ==============================================================================

  get_formatted_labels <- function(labels, level_name) {
    sapply(labels, function(l) {
      is_aggregated <- stringr::str_detect(l, "^[0-9]+")
      if (level_name == "Family" || is_aggregated) {
        return(parse(text = paste0("plain(\"", l, "\")")))
      } else {
        l_clean <- stringr::str_replace_all(l, "\"", "")
        return(parse(text = paste0("italic(\"", l_clean, "\")")))
      }
    })
  }

  plot_detailed_heatmap <- function(df_long, level_name, sample_name, threshold_label) {
    domains <- unique(df_long$Domain)
    ordered_levels <- c()

    for(dom in domains) {
      sub_df <- df_long %>% dplyr::filter(Domain == dom)
      unique_taxa <- unique(sub_df$Taxon)

      recovered_labels <- grep("Recovered only", unique_taxa, value = TRUE)
      recovered_labels <- sort(recovered_labels, decreasing = FALSE)

      lowest_labels <- grep("Lowest abundance", unique_taxa, value = TRUE)

      regular_taxa <- unique_taxa[!unique_taxa %in% c(recovered_labels, lowest_labels)]

      sorted_regular <- c()
      if(length(regular_taxa) > 0) {
        sums <- sub_df %>%
          dplyr::filter(Taxon %in% regular_taxa) %>%
          dplyr::group_by(Taxon) %>%
          dplyr::summarise(T = sum(Rel_Abund)) %>%
          dplyr::arrange(T)
        sorted_regular <- sums$Taxon
      }
      ordered_levels <- c(ordered_levels, recovered_labels, lowest_labels, sorted_regular)
    }

    df_long$Taxon <- factor(df_long$Taxon, levels = unique(ordered_levels))

    n_taxa <- length(unique(df_long$Taxon))
    calc_height <- 3.5 + (n_taxa * 0.18)
    if(calc_height > 100) calc_height <- 100

    p <- ggplot2::ggplot(df_long, ggplot2::aes(x = factor(CS), y = Taxon, fill = Rel_Abund)) +
      ggplot2::geom_tile(color = "white", linewidth = 0.05) +
      ggplot2::facet_grid(Domain ~ ., scales = "free_y", space = "free_y", switch = "y") +
      ggplot2::scale_fill_gradientn(
        colors = c("#ffffff", "#fff7bc", "#fec44f", "#d95f0e", "#993404"),
        values = scales::rescale(c(0, 1, 10, 50, 100)),
        name = "Rel. Abund. (%)",
        labels = scales::label_number(suffix = "%")
      ) +
      ggplot2::scale_y_discrete(labels = function(x) get_formatted_labels(x, level_name)) +
      ggplot2::labs(
        title = paste(level_name, "Relative Abundance -", sample_name),
        subtitle = paste0("Detailed aggregation (Threshold: CS ", threshold_label, ")"),
        x = "Confidence Score",
        y = NULL
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(face = "bold"),
        axis.text.y = ggplot2::element_text(size = 9, color = "black"),
        strip.text.y.left = ggplot2::element_text(angle = 0, face = "bold", size = 11, color="black"),
        strip.placement = "outside",
        strip.background = ggplot2::element_rect(fill = "grey95", color = NA),
        panel.grid = ggplot2::element_blank(),
        panel.spacing = grid::unit(0.1, "cm"),
        legend.position = "right",
        plot.title = ggplot2::element_text(face="bold", hjust = 0.5)
      )
    return(list(plot = p, height = calc_height))
  }

  # ==============================================================================
  # 3. START LOGGING
  # ==============================================================================

  # Initialize Log
  log_con <- file(log_file, open = "wt")
  sink(log_con, type = "output")
  sink(log_con, type = "message")

  # Ensure log closes on exit
  on.exit({
    sink(type = "output")
    sink(type = "message")
    close(log_con)
  }, add = TRUE)

  cat("====================================================\n")
  cat("LOG: 005_HEATMAP_GENERATION\n")
  cat("DATE/TIME:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("PROJECT DIR:", project_dir, "\n")
  cat("====================================================\n\n")

  tryCatch({

    # ==============================================================================
    # 4. INPUT LOGIC
    # ==============================================================================
    cat("STEP 1: Checking Input Data...\n")

    file_rds <- file.path(input_dir, "combined_mpas_matrix_taxon_level.rds")

    INPUT_MATRIX <- NULL

    # 4.1. Check if RDS exists
    if (file.exists(file_rds)) {
      cat("  -> Input file found:", file_rds, "\n")
      cat("  -> Loading RDS object...\n")
      INPUT_MATRIX <- readr::read_rds(file_rds)

    } else {
      cat("  -> [WARNING] Input file NOT found at:", file_rds, "\n")
      cat("  -> Attempting to run Import Function (Step 000)...\n")

      # 4.2. Fallback to Internal Function (STRICT CHANGE FOR PACKAGE)
      tryCatch({
        INPUT_MATRIX <- import_mpa_data(project_dir = project_dir)
      }, error = function(e) {
        stop("CRITICAL: Could not run 'import_mpa_data'. Please ensure the previous step logic is correct.\nError: ", e$message)
      })
    }

    # 4.3. Validate and Log Matrix
    cat("\n  -> MATRIX VALIDATION:\n")
    cat("     Observations:", nrow(INPUT_MATRIX), "\n")
    cat("     Variáveis:", ncol(INPUT_MATRIX), "\n")
    cat("     As variáveis são:\n")

    col_names <- colnames(INPUT_MATRIX)
    for (i in 1:length(col_names)) {
      cat(sprintf("     [%d] %s\n", i, col_names[i]))
    }
    cat("\n")

    # ==============================================================================
    # 5. PARAMETER CONFIGURATION
    # ==============================================================================
    cat("STEP 2: Configuring Analysis Parameters...\n")

    # 5.1. Confidence Score Logic
    if (is.null(confidence_score) || is.na(confidence_score)) {
      last_var_name <- col_names[length(col_names)]
      raw_cs <- stringr::str_extract(last_var_name, "(?<=CS)[0-9.]+")

      if (!is.na(raw_cs)) {
        num_cs <- as.numeric(raw_cs)
        final_cs <- ifelse(num_cs > 1, num_cs/10, num_cs)
      } else {
        final_cs <- 0
        cat("     [WARNING] Could not extract CS from column name. Defaulting to 0.\n")
      }

      confidence_score <- final_cs
      cat("  -> CONFIDENCE_SCORE: Not provided by user.\n")
      cat("     Action: Used last variable reference (", last_var_name, ").\n")
      cat("     Value set to:", confidence_score, "\n")

    } else {
      cat("  -> CONFIDENCE_SCORE: Provided by user.\n")
      cat("     Value set to:", confidence_score, "\n")
    }

    # 5.2. Hardcut (Top N) Logic
    if (is.null(top_n) || is.na(top_n)) {
      top_n <- 20
      cat("  -> TOP_N (Hardcut): Not provided by user.\n")
      cat("     Value set to Default:", top_n, "\n")
    } else {
      cat("  -> TOP_N (Hardcut): Provided by user.\n")
      cat("     Value set to:", top_n, "\n")
    }

    # Define Output Directory
    output_dir <- file.path(project_dir, paste0("005_taxonomic_heatmaps_CS", confidence_score))
    if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    cat("  -> Output Directory:", output_dir, "\n\n")

    # ==============================================================================
    # 6. DATA PROCESSING & PLOTTING
    # ==============================================================================
    cat("STEP 3: Processing Data...\n")

    # 6.1. Header Sanitization
    df_std <- INPUT_MATRIX %>%
      dplyr::rename_with(~ "Taxonomy",  .cols = tidyselect::any_of(c("Taxonomia", "Taxonomy"))) %>%
      dplyr::rename_with(~ "Tax_Level", .cols = tidyselect::any_of(c("tax_level", "Tax_Level", "Nivel_Final")))

    # 6.2. Pivot
    cat("  -> Pivoting to long format...\n")
    df_master <- df_std %>%
      tidyr::pivot_longer(
        cols = -c(Taxonomy, Tax_Level),
        names_to = "Sample_Column_Name",
        values_to = "Counts"
      ) %>%
      dplyr::filter(Counts > 0)

    # 6.3. Metadata Extraction
    df_master <- df_master %>%
      dplyr::mutate(
        Domain = dplyr::case_when(
          stringr::str_detect(Taxonomy, "d__Viruses")   ~ "Viruses",
          stringr::str_detect(Taxonomy, "d__Bacteria")  ~ "Bacteria",
          stringr::str_detect(Taxonomy, "d__Archaea")   ~ "Archaea",
          stringr::str_detect(Taxonomy, "d__Eukaryota") ~ "Eukaryota",
          TRUE ~ "Unknown"
        ),
        Sample_Name = stringr::str_remove(Sample_Column_Name, "_CS[0-9.]+$"),
        CS_Label_Raw = stringr::str_extract(Sample_Column_Name, "CS[0-9.]+$"),
        CS_Val_Raw   = as.numeric(stringr::str_remove(CS_Label_Raw, "CS"))
      ) %>%
      dplyr::mutate(
        CS_Num = ifelse(!is.na(CS_Val_Raw) & CS_Val_Raw > 1, CS_Val_Raw/10, CS_Val_Raw),
        CS_Final_Label = sprintf("%.2f", CS_Num)
      ) %>%
      dplyr::mutate(
        Taxon_Clean_Name = stringr::str_extract(Taxonomy, "[^|]+$"),
        Taxon_Clean_Name = stringr::str_remove(Taxon_Clean_Name, "^[kpcofgs]__")
      )

    unique_samples <- unique(df_master$Sample_Name)
    cat("  -> Unique Samples found:", length(unique_samples), "\n\n")

    # ==============================================================================
    # 7. GENERATION LOOP
    # ==============================================================================
    cat("STEP 4: Generating Heatmaps...\n")

    TARGET_LEVELS <- c("Species", "Genus", "Family")

    for (samp in unique_samples) {

      cat("  Processing Sample:", samp, "...\n")
      samp_out_path <- file.path(output_dir, samp)
      if(!dir.exists(samp_out_path)) dir.create(samp_out_path)

      df_samp <- df_master %>% dplyr::filter(Sample_Name == samp)

      max_cs_in_sample <- max(df_samp$CS_Num, na.rm = TRUE)
      if (max_cs_in_sample < confidence_score) {
        cat("    [SKIP] Sample never reached CS threshold", confidence_score, "\n")
        next
      }

      levels_in_sample <- unique(df_samp$Tax_Level)
      levels_to_plot   <- intersect(levels_in_sample, TARGET_LEVELS)

      for (lvl_name in levels_to_plot) {

        df_lvl <- df_samp %>% dplyr::filter(Tax_Level == lvl_name)
        if(nrow(df_lvl) == 0) next

        df_clean <- df_lvl %>%
          dplyr::group_by(Domain, Taxon_Clean_Name, CS_Final_Label, CS_Num) %>%
          dplyr::summarise(Counts = sum(Counts), .groups="drop") %>%
          dplyr::rename(Taxon_Name = Taxon_Clean_Name, CS = CS_Final_Label)

        max_cs_info <- df_clean %>%
          dplyr::filter(Counts > 0) %>%
          dplyr::group_by(Domain, Taxon_Name) %>%
          dplyr::summarise(Max_CS = max(CS_Num), .groups="drop")

        survivors_list <- max_cs_info %>% dplyr::filter(Max_CS >= confidence_score) %>% dplyr::pull(Taxon_Name)
        non_survivors  <- max_cs_info %>% dplyr::filter(Max_CS < confidence_score)

        plot_data_list <- list()

        # A. Non-Survivors
        if(nrow(non_survivors) > 0) {
          ns_groups <- non_survivors %>%
            dplyr::group_by(Domain, Max_CS) %>%
            dplyr::summarise(Count = dplyr::n(), Taxa = list(Taxon_Name), .groups="drop")

          for(r in 1:nrow(ns_groups)) {
            dom <- ns_groups$Domain[r]
            mcs <- ns_groups$Max_CS[r]
            cnt <- ns_groups$Count[r]
            taxa_in_group <- ns_groups$Taxa[[r]]

            group_label <- paste0(cnt, " ", lvl_name, " Recovered only in CS", sprintf("%.2f", mcs))

            subset_data <- df_clean %>%
              dplyr::filter(Domain == dom, Taxon_Name %in% taxa_in_group) %>%
              dplyr::group_by(Domain, CS, CS_Num) %>%
              dplyr::summarise(Counts = sum(Counts), .groups="drop") %>%
              dplyr::mutate(Taxon = group_label)
            plot_data_list[[length(plot_data_list)+1]] <- subset_data
          }
        }

        # B. Survivors
        unique_domains_here <- unique(df_clean$Domain)

        for(dom in unique_domains_here) {
          dom_survivors <- survivors_list[survivors_list %in% df_clean$Taxon_Name[df_clean$Domain == dom]]
          if(length(dom_survivors) == 0) next

          surv_data_raw <- df_clean %>% dplyr::filter(Domain == dom, Taxon_Name %in% dom_survivors)

          limit <- top_n

          ranking <- surv_data_raw %>%
            dplyr::group_by(Taxon_Name) %>% dplyr::summarise(Total = sum(Counts)) %>% dplyr::arrange(desc(Total))

          if(nrow(ranking) <= limit) {
            plot_data_list[[length(plot_data_list)+1]] <- surv_data_raw %>% dplyr::rename(Taxon = Taxon_Name)
          } else {
            top_names <- ranking$Taxon_Name[1:limit]
            lowest_names <- ranking$Taxon_Name[(limit+1):nrow(ranking)]
            count_lowest <- length(lowest_names)

            min_surv_cs <- min(df_clean$CS_Num[df_clean$CS_Num >= confidence_score])
            low_label <- paste0(count_lowest, " Lowest abundance ", lvl_name, " in CS", sprintf("%.2f", min_surv_cs))

            plot_data_list[[length(plot_data_list)+1]] <- surv_data_raw %>%
              dplyr::filter(Taxon_Name %in% top_names) %>% dplyr::rename(Taxon = Taxon_Name)
            plot_data_list[[length(plot_data_list)+1]] <- surv_data_raw %>%
              dplyr::filter(Taxon_Name %in% lowest_names) %>%
              dplyr::group_by(Domain, CS, CS_Num) %>%
              dplyr::summarise(Counts = sum(Counts), .groups="drop") %>%
              dplyr::mutate(Taxon = low_label)
          }
        }

        if(length(plot_data_list) == 0) next
        df_all_counts <- dplyr::bind_rows(plot_data_list)

        total_reads_per_cs <- df_clean %>%
          dplyr::group_by(Domain, CS) %>%
          dplyr::summarise(Total_Domain = sum(Counts), .groups="drop")

        df_final_perc <- df_all_counts %>%
          dplyr::left_join(total_reads_per_cs, by=c("Domain", "CS")) %>%
          dplyr::mutate(Rel_Abund = (Counts / Total_Domain) * 100) %>%
          tidyr::replace_na(list(Rel_Abund = 0))

        cat("    -> Generating Heatmap:", lvl_name, "\n")
        res <- plot_detailed_heatmap(df_final_perc, lvl_name, samp, confidence_score)
        out_file <- file.path(samp_out_path, paste0(samp, "_", lvl_name, "_Heatmap.pdf"))
        ggplot2::ggsave(out_file, res$plot, width = 8.27, height = 11.69, units = "in")
      }
    }

    cat("\n====================================================\n")
    cat("SUCCESS: All heatmaps generated successfully.\n")
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
