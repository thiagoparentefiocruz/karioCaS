# ==============================================================================
# PRIVATE HELPERS — retrieve_selected_taxa()
# ==============================================================================

#' @noRd
.rst_setup <- function(project_dir) {
    output_dir <- file.path(project_dir, "1000_final_selection")
    log_dir <- file.path(project_dir, "logs")
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)
    log_file <- file.path(log_dir, "log_1000_mosaic_retrieval.txt")
    log_msg <- function(...) {
        msg <- paste0(...)
        message(msg)
        write(
            paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", msg),
            file = log_file, append = TRUE
        )
    }
    list(output_dir = output_dir, log_msg = log_msg)
}

#' @noRd
.rst_load_tse <- function(project_dir, log_msg) {
    log_msg("====================================================")
    log_msg("STEP 1: Loading karioCaS TSE Object...")
    tse_path <- file.path(
        project_dir, "000_karioCaS_input_matrix", "karioCaS_TSE.rds"
    )
    if (!file.exists(tse_path)) {
        log_msg("CRITICAL ERROR: TSE file not found at: ", tse_path)
        stop("TSE file missing.")
    }
    tse <- readRDS(tse_path)
    count_matrix <- SummarizedExperiment::assay(tse, 1)
    row_meta <- SummarizedExperiment::rowData(tse) |> as.data.frame()
    if (!"Taxonomy" %in% names(row_meta)) {
        if ("Taxonomy_Full" %in% names(row_meta)) {
            row_meta$Taxonomy <- row_meta$Taxonomy_Full
        } else {
            stop("TSE RowData missing 'Taxonomy' column.")
        }
    }
    rank_col <- if ("Nivel_Final" %in% names(row_meta)) "Nivel_Final" else "Rank"
    if (rank_col %in% names(row_meta)) row_meta$Rank <- row_meta[[rank_col]]
    row_meta$Domain_Code <- dplyr::case_when(
        stringr::str_detect(row_meta$Taxonomy, "d__Bacteria") ~ "Bacteria",
        stringr::str_detect(row_meta$Taxonomy, "d__Archaea") ~ "Archaea",
        stringr::str_detect(row_meta$Taxonomy, "d__Eukaryota") ~ "Eukaryota",
        stringr::str_detect(row_meta$Taxonomy, "d__Viruses") ~ "Viruses",
        TRUE ~ "Other"
    )
    sample_cs_pattern <- "_CS[0-9.]+$"
    all_cols <- colnames(count_matrix)
    cols_with_cs <- grep(sample_cs_pattern, all_cols, value = TRUE)
    if (length(cols_with_cs) == 0) stop("No '_CSxx' columns found in TSE.")
    SAMPLES <- unique(stringr::str_remove(cols_with_cs, sample_cs_pattern))
    log_msg("  -> Data loaded. Found ", length(SAMPLES), " samples.")
    list(
        count_matrix = count_matrix, row_meta = row_meta,
        SAMPLES = SAMPLES, all_cols = all_cols
    )
}

#' @noRd
.rst_load_audit <- function(project_dir, tax_level, configs, log_msg) {
    needs_audit <- any(
        vapply(configs, function(x) {
            tolower(as.character(x$val)) %in%
                c("auto", "secondary")
        }, logical(1))
    )
    if (!needs_audit) {
        return(NULL)
    }
    log_msg("STEP 1.5: Loading SI Audit Data...")
    audit_tax <- if (is.null(tax_level)) "Species" else tax_level
    audit_file <- file.path(
        project_dir, "006_optimize_CS",
        paste0("SI_Audit_", audit_tax, ".rds")
    )
    if (!file.exists(audit_file)) {
        log_msg("  [WARNING] Audit file not found: ", audit_file)
        log_msg("  Ensure you ran optimize_CS() (Step 006). Fallback to CS = 0.")
        return(NULL)
    }
    audit_df <- readRDS(audit_file)
    log_msg("  -> SI Audit loaded successfully for level: ", audit_tax)
    audit_df
}

#' @noRd
.rst_resolve_cs <- function(user_val, dom, samp, audit_df, log_msg) {
    uv <- tolower(as.character(user_val))
    if (!uv %in% c("auto", "secondary")) {
        return(list(val = as.numeric(uv), tag = "[Manual]"))
    }
    if (!is.null(audit_df)) {
        sub <- dplyr::filter(audit_df, .data$Sample == samp, .data$Domain == dom)
        if (nrow(sub) > 0) {
            if (uv == "auto") {
                v <- sub$CS[which(sub$SI_Type == "Primary_SI")]
                if (length(v) > 0) {
                    return(list(val = v[1], tag = "[SI: Primary]"))
                }
            } else {
                v <- sub$CS[which(sub$SI_Type == "Secondary_SI_1")]
                if (length(v) > 0 && !is.na(v[1])) {
                    return(list(val = v[1], tag = "[SI: Secondary]"))
                }
                v <- sub$CS[which(sub$SI_Type == "Primary_SI")]
                log_msg(sprintf(
                    "    [INFO] No Secondary SI for %s in %s. Using Primary.", dom, samp
                ))
                if (length(v) > 0) {
                    return(list(val = v[1], tag = "[SI: Fallback to Primary]"))
                }
            }
        }
    }
    list(val = 0, tag = "[SI: Audit Fail -> CS0]")
}

#' @noRd
.rst_match_column <- function(samp, requested_val, available_suffixes, log_msg) {
    for (suf in available_suffixes) {
        suf_num <- tryCatch(as.numeric(suf), warning = function(w) NA_real_)
        if (is.na(suf_num)) next
        if (suf_num == requested_val) {
            return(suf)
        }
        if ((requested_val * 10) == suf_num) {
            return(suf)
        }
        if (requested_val == (suf_num * 10)) {
            return(suf)
        }
    }
    log_msg(sprintf(
        "    [ERROR] CS input (Resolved: %s) not matched in available: %s",
        requested_val, paste(available_suffixes, collapse = ", ")
    ))
    NULL
}

#' @noRd
.rst_filter_domain <- function(count_matrix, row_meta, target_col,
                               dom, tax_level, min_reads) {
    counts_vec <- count_matrix[, target_col]
    pass <- which(counts_vec >= min_reads & counts_vec > 0)
    if (length(pass) == 0) {
        return(NULL)
    }
    sub_meta <- row_meta[pass, ]
    sub_counts <- counts_vec[pass]
    mask <- sub_meta$Domain_Code == dom
    if (!is.null(tax_level)) mask <- mask & (sub_meta$Rank == tax_level)
    mask[is.na(mask)] <- FALSE
    if (sum(mask) == 0) {
        return(NULL)
    }
    data.frame(
        Taxonomy = sub_meta$Taxonomy[mask],
        Counts = sub_counts[mask],
        stringsAsFactors = FALSE
    )
}

#' @noRd
.rst_process_sample <- function(samp, count_matrix, row_meta, all_cols,
                                configs, audit_df, tax_level,
                                output_dir, log_msg) {
    log_msg("----------------------------------------------------")
    log_msg("  Sample: ", samp)
    samp_cols <- grep(paste0("^", samp, "_CS"), all_cols, value = TRUE)
    available_suffixes <- stringr::str_remove(samp_cols, paste0("^", samp, "_CS"))
    taxa_list <- list()
    for (dom in names(configs)) {
        cfg <- configs[[dom]]
        resolved <- .rst_resolve_cs(cfg$val, dom, samp, audit_df, log_msg)
        match_suf <- .rst_match_column(samp, resolved$val, available_suffixes, log_msg)
        if (is.null(match_suf)) next
        target_col <- paste0(samp, "_CS", match_suf)
        part_df <- .rst_filter_domain(
            count_matrix, row_meta, target_col, dom, tax_level, cfg$min
        )
        if (is.null(part_df)) next
        cs_num <- tryCatch(as.numeric(match_suf), warning = function(w) NA_real_)
        cs_display <- if (!is.na(cs_num)) cs_num / 100 else 0
        log_msg(sprintf(
            "    -> Added %d %s taxa (CS = %.2f %s | min_reads = %d)",
            nrow(part_df), dom, cs_display, resolved$tag, cfg$min
        ))
        taxa_list[[dom]] <- part_df
    }
    if (length(taxa_list) == 0) {
        log_msg("    -> FAILED: No output generated for ", samp)
        return(FALSE)
    }
    final_df <- dplyr::bind_rows(taxa_list) |>
        dplyr::group_by(.data$Taxonomy) |>
        dplyr::summarise(Counts = sum(.data$Counts), .groups = "drop") |>
        dplyr::rename(!!samp := "Counts")
    base_name <- paste0(samp, "_karioCaS_Mosaic")
    readr::write_delim(final_df, file.path(output_dir, paste0(base_name, ".mpa")),
        delim = "\t"
    )
    readr::write_tsv(final_df, file.path(output_dir, paste0(base_name, ".tsv")))
    log_msg("    -> GENERATED: ", base_name, ".mpa")
    TRUE
}

# ==============================================================================
# EXPORTED FUNCTION
# ==============================================================================

#' Retrieve Selected Taxa with Domain-Specific Thresholds (Final Mosaic Step)
#'
#' Creates a "biological mosaic" for each sample using the
#' \code{karioCaS_TSE.rds} object from Step 000 and, optionally, the
#' Stability Index audit from Step 006. Supports \code{"auto"},
#' \code{"secondary"}, or manual numeric thresholds per domain.
#' Enforces a strict \code{> 0} read filter to prevent zero-count taxa.
#'
#' @param project_dir Path to the project root.
#' @param tax_level Taxonomic level to filter (e.g., \code{"Species"}).
#'   If \code{NULL}, all levels are retained.
#' @param CS_A Character or numeric. CS for Archaea:
#'   \code{"auto"}, \code{"secondary"}, or a numeric value. Default: \code{"auto"}.
#' @param reads_min_A Integer. Minimum reads for Archaea. Default: 0.
#' @param CS_B Character or numeric. CS for Bacteria. Default: \code{"auto"}.
#' @param reads_min_B Integer. Minimum reads for Bacteria. Default: 0.
#' @param CS_E Character or numeric. CS for Eukaryota. Default: \code{"auto"}.
#' @param reads_min_E Integer. Minimum reads for Eukaryota. Default: 0.
#' @param CS_V Character or numeric. CS for Viruses. Default: \code{"auto"}.
#' @param reads_min_V Integer. Minimum reads for Viruses. Default: 0.
#'
#' @return Invisibly returns \code{TRUE}. Mosaic \code{.mpa} and \code{.tsv}
#'   files are saved to \code{<project_dir>/1000_final_selection/}.
#' @export
#' @importFrom readr read_rds write_tsv write_delim
#' @importFrom dplyr filter mutate select group_by summarise bind_rows
#'   rename left_join case_when
#' @importFrom stringr str_detect str_remove str_extract
#' @importFrom SummarizedExperiment assay rowData
#' @examples
#' toy_project <- system.file("extdata", "your_project_name", package = "karioCaS")
#'
#' # retrieve_selected_taxa(
#' #   project_dir = toy_project,
#' #   tax_level   = "Species",
#' #   CS_B        = "auto",
#' #   CS_A        = 20,
#' #   CS_E        = 40,
#' #   CS_V        = 0
#' # )
retrieve_selected_taxa <- function(project_dir,
                                   tax_level = NULL,
                                   CS_A = "auto", reads_min_A = 0,
                                   CS_B = "auto", reads_min_B = 0,
                                   CS_E = "auto", reads_min_E = 0,
                                   CS_V = "auto", reads_min_V = 0) {
    setup <- .rst_setup(project_dir)
    log_msg <- setup$log_msg
    configs <- list(
        Archaea   = list(val = CS_A, min = as.numeric(reads_min_A)),
        Bacteria  = list(val = CS_B, min = as.numeric(reads_min_B)),
        Eukaryota = list(val = CS_E, min = as.numeric(reads_min_E)),
        Viruses   = list(val = CS_V, min = as.numeric(reads_min_V))
    )
    tryCatch(
        {
            tse_data <- .rst_load_tse(project_dir, log_msg)
            audit_df <- .rst_load_audit(project_dir, tax_level, configs, log_msg)
            log_msg("STEP 2: Processing Mosaics...")
            n_generated <- sum(vapply(tse_data$SAMPLES, function(samp) {
                .rst_process_sample(
                    samp, tse_data$count_matrix, tse_data$row_meta,
                    tse_data$all_cols, configs, audit_df, tax_level,
                    setup$output_dir, log_msg
                )
            }, logical(1)))
            if (n_generated > 0) {
                log_msg(
                    "\nSUCCESS: Process completed. Generated ",
                    n_generated, " mosaic files."
                )
            } else {
                log_msg("\nFAILURE: No files were generated.")
            }
            invisible(TRUE)
        },
        error = function(e) {
            write(paste0("\nCRITICAL ERROR: ", e$message),
                file = file.path(
                    project_dir, "logs",
                    "log_1000_mosaic_retrieval.txt"
                ),
                append = TRUE
            )
            stop(e$message)
        }
    )
}
