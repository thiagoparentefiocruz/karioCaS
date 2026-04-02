#' Retrieve Selected Taxa with Domain-Specific Thresholds (Final Mosaic Step)
#'
#' This function creates a "biological mosaic" for each sample using 'karioCaS_TSE.rds'.
#' It seamlessly integrates with the Stability Index (SI) audit file from Step 006.
#' Users can specify "auto" (Primary SI), "secondary" (Secondary SI), or manual numeric
#' thresholds for each domain.
#' ENFORCES strict > 0 filter to prevent zero-count taxa from appearing in output.
#'
#' @param project_dir Path to the project root.
#' @param tax_level Specific taxonomic level to filter (e.g., "Species", "Genus"). If NULL, keeps all levels.
#' @param CS_A Character/Numeric. CS for Archaea: "auto", "secondary", or numeric value. Default: "auto".
#' @param reads_min_A Integer. Min reads for Archaea. Default: 0.
#' @param CS_B Character/Numeric. CS for Bacteria: "auto", "secondary", or numeric value. Default: "auto".
#' @param reads_min_B Integer. Min reads for Bacteria. Default: 0.
#' @param CS_E Character/Numeric. CS for Eukaryota: "auto", "secondary", or numeric value. Default: "auto".
#' @param reads_min_E Integer. Min reads for Eukaryota. Default: 0.
#' @param CS_V Character/Numeric. CS for Viruses: "auto", "secondary", or numeric value. Default: "auto".
#' @param reads_min_V Integer. Min reads for Viruses. Default: 0.
#'
#' @export
#' @importFrom readr read_rds write_tsv write_delim
#' @importFrom dplyr filter mutate select group_by summarise bind_rows rename left_join case_when
#' @importFrom stringr str_detect str_remove str_extract
#' @importFrom SummarizedExperiment assay rowData

retrieve_selected_taxa <- function(project_dir,
                                   tax_level = NULL,
                                   CS_A = "auto", reads_min_A = 0,
                                   CS_B = "auto", reads_min_B = 0,
                                   CS_E = "auto", reads_min_E = 0,
                                   CS_V = "auto", reads_min_V = 0) {
  
  # 1. SETUP & LOGGING
  output_dir  <- file.path(project_dir, "1000_final_selection")
  log_dir     <- file.path(project_dir, "logs")
  
  if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  if(!dir.exists(log_dir))    dir.create(log_dir, recursive = TRUE)
  
  log_file <- file.path(log_dir, "log_1000_mosaic_retrieval.txt")
  
  log_msg <- function(...) {
    msg <- paste0(...)
    message(msg)
    cat(paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", msg, "\n"), file = log_file, append = TRUE)
  }
  
  files_generated_count <- 0
  
  tryCatch({
    
    # 2. LOAD TSE DATA
    log_msg("====================================================")
    log_msg("STEP 1: Loading karioCaS TSE Object...")
    
    file_tse <- file.path(project_dir, "000_karioCaS_input_matrix", "karioCaS_TSE.rds")
    
    if (!file.exists(file_tse)) {
      log_msg("CRITICAL ERROR: TSE file not found at: ", file_tse)
      stop("TSE file missing.")
    }
    
    tse <- readRDS(file_tse)
    
    # Extract matrix and metadata
    count_matrix <- SummarizedExperiment::assay(tse, 1)
    row_meta <- SummarizedExperiment::rowData(tse) %>% as.data.frame()
    
    # Metadata Validation
    if (!"Taxonomy" %in% names(row_meta)) {
      if ("Taxonomy_Full" %in% names(row_meta)) {
        row_meta$Taxonomy <- row_meta$Taxonomy_Full
      } else {
        stop("TSE RowData missing 'Taxonomy' column.")
      }
    }
    
    # Normalize Rank
    rank_col <- if("Nivel_Final" %in% names(row_meta)) "Nivel_Final" else "Rank"
    if (rank_col %in% names(row_meta)) {
      row_meta$Rank <- row_meta[[rank_col]]
    }
    
    # Domain Classification
    row_meta$Domain_Code <- dplyr::case_when(
      stringr::str_detect(row_meta$Taxonomy, "d__Bacteria") ~ "Bacteria",
      stringr::str_detect(row_meta$Taxonomy, "d__Archaea") ~ "Archaea",
      stringr::str_detect(row_meta$Taxonomy, "d__Eukaryota") ~ "Eukaryota",
      stringr::str_detect(row_meta$Taxonomy, "d__Viruses") ~ "Viruses",
      TRUE ~ "Other"
    )
    
    # 3. IDENTIFY SAMPLES
    all_cols <- colnames(count_matrix)
    sample_cs_pattern <- "_CS[0-9.]+$"
    cols_with_cs <- grep(sample_cs_pattern, all_cols, value = TRUE)
    
    if(length(cols_with_cs) == 0) stop("No '_CSxx' columns found in TSE.")
    SAMPLES <- unique(stringr::str_remove(cols_with_cs, sample_cs_pattern))
    log_msg("  -> Data loaded. Found ", length(SAMPLES), " samples.")
    
    # 4. LOAD AUDIT DATA (If needed)
    audit_df <- NULL
    needs_audit <- any(c(CS_A, CS_B, CS_E, CS_V) %in% c("auto", "secondary"))
    
    if (needs_audit) {
      log_msg("STEP 1.5: Loading SI Audit Data...")
      # Assume audit matches the requested tax_level, or defaults to "Species" if NULL
      audit_tax <- if(is.null(tax_level)) "Species" else tax_level
      audit_file <- file.path(project_dir, "006_optimize_CS", paste0("SI_Audit_", audit_tax, ".rds"))
      
      if (file.exists(audit_file)) {
        audit_df <- readRDS(audit_file)
        log_msg("  -> SI Audit loaded successfully for level: ", audit_tax)
      } else {
        log_msg("  [WARNING] Audit file not found: ", audit_file)
        log_msg("  Ensure you ran optimize_CS() (Step 006). Fallback to CS = 0 for auto/secondary requests.")
      }
    }
    
    # 5. CONFIG
    configs <- list(
      "Archaea"   = list(val = CS_A, min = as.numeric(reads_min_A)),
      "Bacteria"  = list(val = CS_B, min = as.numeric(reads_min_B)),
      "Eukaryota" = list(val = CS_E, min = as.numeric(reads_min_E)),
      "Viruses"   = list(val = CS_V, min = as.numeric(reads_min_V))
    )
    
    # 6. MOSAIC LOOP
    log_msg("STEP 2: Processing Mosaics...")
    
    for (samp in SAMPLES) {
      log_msg("----------------------------------------------------")
      log_msg("  Sample: ", samp)
      
      samp_cols <- grep(paste0("^", samp, "_CS"), all_cols, value = TRUE)
      available_suffixes <- stringr::str_remove(samp_cols, paste0("^", samp, "_CS"))
      
      sample_taxa_list <- list()
      
      for (dom in names(configs)) {
        cfg <- configs[[dom]]
        user_val <- tolower(as.character(cfg$val))
        
        requested_val <- NULL
        info_tag <- ""
        
        # --- DYNAMIC THRESHOLD RESOLUTION ---
        if (user_val %in% c("auto", "secondary")) {
          
          if (!is.null(audit_df)) {
            # Filter audit for current Sample and Domain
            audit_sub <- audit_df %>% dplyr::filter(Sample == samp, Domain == dom)
            
            if (nrow(audit_sub) > 0) {
              if (user_val == "auto") {
                requested_val <- audit_sub$CS[audit_sub$SI_Type == "Primary_SI"]
                info_tag <- "[SI: Primary]"
              } else if (user_val == "secondary") {
                sec_val <- audit_sub$CS[audit_sub$SI_Type == "Secondary_SI_1"]
                if (length(sec_val) > 0 && !is.na(sec_val[1])) {
                  requested_val <- sec_val
                  info_tag <- "[SI: Secondary]"
                } else {
                  # Fallback if no secondary exists
                  requested_val <- audit_sub$CS[audit_sub$SI_Type == "Primary_SI"]
                  info_tag <- "[SI: Fallback to Primary]"
                  log_msg(sprintf("    [INFO] No Secondary SI for %s in %s. Using Primary.", dom, samp))
                }
              }
              requested_val <- requested_val[1] # Ensure single value
            }
          }
          
          # Ultimate fallback if audit lookup fails completely
          if (is.null(requested_val) || is.na(requested_val)) {
            requested_val <- 0
            info_tag <- "[SI: Audit Fail -> CS0]"
          }
          
        } else {
          # Manual Numeric Input
          requested_val <- as.numeric(user_val)
          info_tag <- "[Manual]"
        }
        
        # --- FLEXIBLE MATCHING LOGIC ---
        target_col <- NULL
        match_suffix <- NULL
        
        for (suf in available_suffixes) {
          suf_num <- as.numeric(suf)
          
          if (!is.na(suf_num) && suf_num == requested_val) {
            match_suffix <- suf
            break
          }
          if (!is.na(suf_num) && (requested_val * 10) == suf_num) {
            match_suffix <- suf
            break
          }
          if (!is.na(suf_num) && requested_val == (suf_num * 10)) {
            match_suffix <- suf
            break
          }
        }
        
        if (!is.null(match_suffix)) {
          target_col <- paste0(samp, "_CS", match_suffix)
        } else {
          log_msg(sprintf("    [ERROR] CS input '%s' (Resolved: %s) not matched in available: %s", 
                          user_val, requested_val, paste(available_suffixes, collapse=", ")))
          next
        }
        
        # --- FILTERING (ZERO READS CORRECTION) ---
        counts_vec <- count_matrix[, target_col]
        pass_reads <- which(counts_vec >= cfg$min & counts_vec > 0)
        
        if (length(pass_reads) == 0) next
        
        subset_meta <- row_meta[pass_reads, ]
        subset_counts <- counts_vec[pass_reads]
        
        is_domain <- subset_meta$Domain_Code == dom
        
        is_rank <- TRUE
        if (!is.null(tax_level)) {
          is_rank <- subset_meta$Rank == tax_level
        }
        
        final_mask <- is_domain & is_rank
        final_mask[is.na(final_mask)] <- FALSE
        
        if (sum(final_mask) > 0) {
          part_df <- data.frame(
            Taxonomy = subset_meta$Taxonomy[final_mask],
            Counts   = subset_counts[final_mask],
            stringsAsFactors = FALSE
          )
          sample_taxa_list[[dom]] <- part_df
          
          # Custom Log
          cs_display <- as.numeric(match_suffix) / 100
          log_msg(sprintf("    -> Added %d %s taxa (CS = %.2f %s | min_reads = %d)",
                          nrow(part_df), dom, cs_display, info_tag, cfg$min))
        }
      }
      
      # 7. SAVE OUTPUT
      if (length(sample_taxa_list) > 0) {
        final_df <- dplyr::bind_rows(sample_taxa_list) %>%
          dplyr::group_by(Taxonomy) %>%
          dplyr::summarise(Counts = sum(Counts), .groups = "drop") %>%
          dplyr::rename(!!samp := Counts)
        
        base_name <- paste0(samp, "_karioCaS_Mosaic")
        file_path_mpa <- file.path(output_dir, paste0(base_name, ".mpa"))
        file_path_tsv <- file.path(output_dir, paste0(base_name, ".tsv"))
        
        readr::write_delim(final_df, file_path_mpa, delim = "\t")
        readr::write_tsv(final_df, file_path_tsv)
        
        log_msg("    -> GENERATED: ", base_name, ".mpa")
        files_generated_count <- files_generated_count + 1
        
      } else {
        log_msg("    -> FAILED: No output generated for ", samp)
      }
    }
    
    if (files_generated_count > 0) {
      log_msg("\nSUCCESS: Process completed. Generated ", files_generated_count, " mosaic files.")
    } else {
      log_msg("\nFAILURE: No files were generated.")
    }
    
    return(invisible(TRUE))
    
  }, error = function(e) {
    cat("\nCRITICAL ERROR: ", e$message, "\n", file = log_file, append = TRUE)
    stop(e$message)
  })
}