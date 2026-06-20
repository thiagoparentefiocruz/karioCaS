# ==============================================================================
# PRIVATE HELPERS - retrieve_selected_taxa()
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
    audit_name <- paste0("SI_Audit_", audit_tax, ".rds")
    audit_file <- file.path(project_dir, "001_taxa_retention", audit_name)
    if (!file.exists(audit_file)) {
        # Backward compatibility: audit used to live in 006_optimize_CS/.
        legacy_file <- file.path(project_dir, "006_optimize_CS", audit_name)
        if (file.exists(legacy_file)) {
            audit_file <- legacy_file
        } else {
            log_msg("  [WARNING] Audit file not found: ", audit_file)
            log_msg(
                "  Ensure you ran taxa_retention() (Step 001).",
                " Fallback to CS = 0."
            )
            return(NULL)
        }
    }
    audit_df <- readRDS(audit_file)
    log_msg("  -> SI Audit loaded successfully for level: ", audit_tax)
    audit_df
}

#' @noRd
.rst_resolve_cs <- function(user_val, dom, samp, audit_df, log_msg) {
    uv <- tolower(as.character(user_val))
    if (!uv %in% c("auto", "secondary")) {
        return(list(val = .cs_arg_to_percent(uv), tag = "[Manual]"))
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
    target_pct <- .cs_arg_to_percent(requested_val)
    if (is.na(target_pct)) {
        log_msg(sprintf(
            "    [ERROR] Invalid CS input (%s) for sample %s.",
            requested_val, samp
        ))
        return(NULL)
    }
    for (suf in available_suffixes) {
        # Column suffixes are already stored as canonical integer percent
        # (e.g. "00", "90", "100"); parse them directly, no re-encoding.
        suf_pct <- suppressWarnings(as.numeric(suf))
        if (!is.na(suf_pct) && suf_pct == target_pct) {
            return(suf)
        }
    }
    log_msg(sprintf(
        "    [ERROR] CS input (Resolved: %d%%) not matched in available: %s",
        target_pct, paste(available_suffixes, collapse = ", ")
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
.rst_load_reads_audit <- function(project_dir, tax_level, configs, log_msg) {
    needs <- any(vapply(configs, function(x) {
        tolower(as.character(x$min_val)) %in% c("auto", "secondary")
    }, logical(1)))
    if (!needs) {
        return(NULL)
    }
    log_msg("STEP 1.6: Loading Reads Audit Data...")
    audit_tax <- if (is.null(tax_level)) "Species" else tax_level
    f <- file.path(
        project_dir, "003_cutoffs", paste0("Reads_Audit_", audit_tax, ".rds")
    )
    if (!file.exists(f)) {
        log_msg("  [WARNING] Reads audit not found: ", f)
        log_msg("  Ensure you ran reads_per_taxa() (Step 003). Fallback to 0.")
        return(NULL)
    }
    log_msg("  -> Reads audit loaded for level: ", audit_tax)
    readRDS(f)
}

#' @noRd
.rst_resolve_reads <- function(user_val, dom, samp, resolved_cs,
                               reads_audit, log_msg) {
    uv <- tolower(as.character(user_val))
    if (!uv %in% c("auto", "secondary")) {
        n <- suppressWarnings(as.numeric(uv))
        return(list(val = if (is.na(n)) 0 else n, tag = "[Manual]"))
    }
    if (!is.null(reads_audit)) {
        sub <- dplyr::filter(
            reads_audit, .data$Sample == samp,
            .data$Domain == dom, .data$CS == resolved_cs
        )
        if (nrow(sub) > 0) {
            want <- if (uv == "auto") "Primary_SI" else "Secondary_SI_1"
            v <- sub$Cutoff[which(sub$SI_Type == want)]
            if (length(v) > 0 && !is.na(v[1])) {
                return(list(val = v[1], tag = paste0("[Reads:", uv, "]")))
            }
            if (uv == "secondary") {
                v <- sub$Cutoff[which(sub$SI_Type == "Primary_SI")]
                if (length(v) > 0) {
                    return(list(val = v[1], tag = "[Reads:Sec->Primary]"))
                }
            }
        }
    }
    log_msg(sprintf(
        "    [INFO] No optimal reads for %s @ CS%02d in %s; using 0.",
        dom, resolved_cs, samp
    ))
    list(val = 0, tag = "[Reads:none->0]")
}

#' @noRd
.rst_process_sample <- function(samp, count_matrix, row_meta, all_cols,
                                configs, audit_df, reads_audit, tax_level,
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
        reads_res <- .rst_resolve_reads(
            cfg$min_val, dom, samp, resolved$val, reads_audit, log_msg
        )
        part_df <- .rst_filter_domain(
            count_matrix, row_meta, target_col, dom, tax_level, reads_res$val
        )
        if (is.null(part_df)) next
        cs_num <- tryCatch(as.numeric(match_suf), warning = function(w) NA_real_)
        cs_display <- if (!is.na(cs_num)) cs_num / 100 else 0
        log_msg(sprintf(
            "    -> Added %d %s taxa (CS = %.2f %s | min_reads = %g %s)",
            nrow(part_df), dom, cs_display, resolved$tag,
            reads_res$val, reads_res$tag
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
#' \code{karioCaS_TSE.rds} object from Step 000 and, optionally, the optimal
#' thresholds computed earlier: the optimal Confidence Score (Stability Index
#' audit from \code{taxa_retention()}, Step 001) and the optimal minimum reads
#' (\code{Reads_Audit} from \code{reads_per_taxa()}, Step 003). Both the
#' \code{CS_*} and \code{reads_min_*} arguments accept \code{"auto"},
#' \code{"secondary"}, or a manual numeric value per domain. The optimal reads is
#' looked up at each domain's resolved CS, so the mosaic combines both data-driven
#' thresholds automatically. Enforces a strict \code{> 0} read filter.
#'
#' @param project_dir Path to the project root.
#' @param tax_level Taxonomic level to filter (e.g., \code{"Species"}).
#'   If \code{NULL}, all levels are retained. Also selects which
#'   \code{SI_Audit}/\code{Reads_Audit} (\code{"Species"} when \code{NULL}).
#' @param CS_A Character or numeric. CS for Archaea:
#'   \code{"auto"}, \code{"secondary"}, or a numeric value. Default: \code{"auto"}.
#' @param reads_min_A Minimum reads for Archaea: \code{"auto"}/\code{"secondary"}
#'   (pulled from the Reads_Audit at the resolved CS) or a numeric value.
#'   Default: 0.
#' @param CS_B Character or numeric. CS for Bacteria. Default: \code{"auto"}.
#' @param reads_min_B Minimum reads for Bacteria. Default: 0.
#' @param CS_E Character or numeric. CS for Eukaryota. Default: \code{"auto"}.
#' @param reads_min_E Minimum reads for Eukaryota. Default: 0.
#' @param CS_V Character or numeric. CS for Viruses. Default: \code{"auto"}.
#' @param reads_min_V Minimum reads for Viruses. Default: 0.
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
#' # Fully data-driven mosaic: optimal CS and optimal min-reads per domain
#' # retrieve_selected_taxa(
#' #   project_dir = toy_project,
#' #   tax_level   = "Species",
#' #   CS_B        = "auto", reads_min_B = "auto",
#' #   CS_A        = "auto", reads_min_A = "auto",
#' #   CS_E        = 40,     reads_min_E = 10,
#' #   CS_V        = 0,      reads_min_V = 0
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
        Archaea   = list(val = CS_A, min_val = reads_min_A),
        Bacteria  = list(val = CS_B, min_val = reads_min_B),
        Eukaryota = list(val = CS_E, min_val = reads_min_E),
        Viruses   = list(val = CS_V, min_val = reads_min_V)
    )
    tryCatch(
        {
            tse_data <- .rst_load_tse(project_dir, log_msg)
            audit_df <- .rst_load_audit(project_dir, tax_level, configs, log_msg)
            reads_audit <- .rst_load_reads_audit(
                project_dir, tax_level, configs, log_msg
            )
            log_msg("STEP 2: Processing Mosaics...")
            n_generated <- sum(vapply(tse_data$SAMPLES, function(samp) {
                .rst_process_sample(
                    samp, tse_data$count_matrix, tse_data$row_meta,
                    tse_data$all_cols, configs, audit_df, reads_audit,
                    tax_level, setup$output_dir, log_msg
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
