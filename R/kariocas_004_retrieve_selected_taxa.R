# ==============================================================================
# PRIVATE HELPERS - retrieve_selected_taxa()
# ==============================================================================

#' @noRd
.rst_setup <- function(project_dir, export = TRUE) {
    .kcs_setup_step(
        project_dir, "004_final_mosaic", "log_004_final_mosaic.txt",
        export = export
    )
}

#' Validate the per-domain threshold arguments up front
#' @noRd
.rst_validate_configs <- function(configs) {
    for (dom in names(configs)) {
        for (arg in c("val", "min_val")) {
            v <- configs[[dom]][[arg]]
            label <- if (arg == "val") "CS_" else "reads_min_"
            label <- paste0(label, substr(dom, 1, 1))
            if (length(v) != 1 || is.na(v)) {
                stop("'", label, "' must be a single value.", call. = FALSE)
            }
            if (tolower(as.character(v)) %in% c("auto", "secondary")) next
            n <- suppressWarnings(as.numeric(v))
            if (is.na(n) || n < 0) {
                stop(
                    "Invalid '", label, "': ", v,
                    ". Use \"auto\", \"secondary\" or a non-negative number.",
                    call. = FALSE
                )
            }
            if (arg == "val" && is.na(.cs_arg_to_percent(v))) {
                stop(
                    "Invalid '", label, "': ", v,
                    ". Use a Kraken fraction (0-1) or a percentage (0-100).",
                    call. = FALSE
                )
            }
        }
    }
    invisible(TRUE)
}

#' @noRd
.rst_load_tse <- function(project_dir, log_msg) {
    log_msg("====================================================")
    log_msg("STEP 1: Loading karioCaS TSE Object...")
    tse_path <- .kcs_path(
        project_dir,
        file.path("001_imported_matrix", "karioCaS_TSE.rds"),
        file.path("000_karioCaS_input_matrix", "karioCaS_TSE.rds")
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
    # Canonical 002_taxa_retention/, with fallbacks to earlier folder names.
    audit_file <- .kcs_path(
        project_dir,
        file.path("002_taxa_retention", audit_name),
        file.path("001_taxa_retention", audit_name),
        file.path("006_optimize_CS", audit_name)
    )
    if (!file.exists(audit_file)) {
        stop(
            "CS_* = \"auto\"/\"secondary\" needs the Stability Index audit ",
            "'", audit_name, "', which was not found in ", project_dir, ". ",
            "Run taxa_retention(project_dir, tax_level = \"", audit_tax,
            "\", export = TRUE) first, or give numeric CS values.",
            call. = FALSE
        )
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
                               dom, min_reads) {
    counts_vec <- count_matrix[, target_col]
    pass <- which(counts_vec >= min_reads & counts_vec > 0)
    if (length(pass) == 0) {
        return(NULL)
    }
    sub_meta <- row_meta[pass, ]
    sub_counts <- counts_vec[pass]
    # The mosaic always keeps ALL taxonomic ranks present in the .mpa (no rank
    # filter); tax_level only governs which optimization audit "auto" reads.
    mask <- sub_meta$Domain_Code == dom
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
    reads_name <- paste0("Reads_Audit_", audit_tax, ".rds")
    f <- .kcs_path(
        project_dir,
        file.path("003_reads_saturation", reads_name),
        file.path("003_cutoffs", reads_name)
    )
    if (!file.exists(f)) {
        stop(
            "reads_min_* = \"auto\"/\"secondary\" needs the reads audit ",
            "'", reads_name, "', which was not found in ", project_dir, ". ",
            "Run reads_per_taxa(project_dir, analysis_level = \"", audit_tax,
            "\", export = TRUE) first, or give numeric minimum reads.",
            call. = FALSE
        )
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
                                configs, audit_df, reads_audit, setup) {
    log_msg <- setup$log_msg
    log_msg("----------------------------------------------------")
    log_msg("  Sample: ", samp)
    issues <- character(0)
    samp_cols <- grep(paste0("^", samp, "_CS"), all_cols, value = TRUE)
    available_suffixes <- stringr::str_remove(samp_cols, paste0("^", samp, "_CS"))
    taxa_list <- list()
    for (dom in names(configs)) {
        cfg <- configs[[dom]]
        resolved <- .rst_resolve_cs(cfg$val, dom, samp, audit_df, log_msg)
        if (grepl("Fail", resolved$tag)) {
            issues <- c(issues, paste0(
                samp, "/", dom, ": no ", cfg$val, " CS in the SI audit, used CS 0"
            ))
        }
        match_suf <- .rst_match_column(samp, resolved$val, available_suffixes, log_msg)
        if (is.null(match_suf)) {
            issues <- c(issues, paste0(
                samp, "/", dom, ": CS ", resolved$val,
                " not available, domain skipped"
            ))
            next
        }
        target_col <- paste0(samp, "_CS", match_suf)
        reads_res <- .rst_resolve_reads(
            cfg$min_val, dom, samp, resolved$val, reads_audit, log_msg
        )
        if (grepl("none", reads_res$tag)) {
            issues <- c(issues, paste0(
                samp, "/", dom, ": no ", cfg$min_val,
                " minimum reads at CS ", resolved$val, ", used 0"
            ))
        }
        part_df <- .rst_filter_domain(
            count_matrix, row_meta, target_col, dom, reads_res$val
        )
        if (is.null(part_df)) next
        cs_num <- tryCatch(as.numeric(match_suf), warning = function(w) NA_real_)
        cs_display <- if (!is.na(cs_num)) cs_num / 100 else 0
        log_msg(sprintf(
            "    -> Added %d %s taxa (CS = %.2f %s | min_reads = %g %s)",
            nrow(part_df), dom, cs_display, resolved$tag,
            reads_res$val, reads_res$tag
        ))
        taxa_list[[dom]] <- dplyr::mutate(
            part_df,
            Domain = dom, CS = resolved$val, CS_Source = resolved$tag,
            Min_Reads = reads_res$val, Min_Reads_Source = reads_res$tag
        )
    }
    if (length(taxa_list) == 0) {
        log_msg("    -> FAILED: No output generated for ", samp)
        return(list(data = NULL, paths = character(0), issues = issues))
    }
    long_df <- dplyr::bind_rows(taxa_list)
    final_df <- long_df |>
        dplyr::group_by(.data$Taxonomy) |>
        dplyr::summarise(Counts = sum(.data$Counts), .groups = "drop") |>
        dplyr::rename(!!samp := "Counts")
    paths <- character(0)
    if (isTRUE(setup$export)) {
        base_name <- paste0(samp, "_karioCaS_Mosaic")
        mpa_dir <- file.path(setup$output_dir, "mpa")
        tsv_dir <- file.path(setup$output_dir, "tsv")
        if (!dir.exists(mpa_dir)) dir.create(mpa_dir, recursive = TRUE)
        if (!dir.exists(tsv_dir)) dir.create(tsv_dir, recursive = TRUE)
        paths <- c(
            file.path(mpa_dir, paste0(base_name, ".mpa")),
            file.path(tsv_dir, paste0(base_name, ".tsv"))
        )
        readr::write_delim(final_df, paths[1], delim = "\t")
        readr::write_tsv(final_df, paths[2])
        log_msg("    -> GENERATED: mpa/", base_name, ".mpa")
    }
    list(
        data = dplyr::mutate(long_df, sample = samp, .before = 1),
        paths = paths, issues = issues
    )
}

# ==============================================================================
# EXPORTED FUNCTION
# ==============================================================================

#' Retrieve Selected Taxa with Domain-Specific Thresholds (Step 004)
#'
#' Creates a "biological mosaic" for each sample using the
#' \code{karioCaS_TSE.rds} object from Step 001 and, optionally, the optimal
#' thresholds computed earlier: the optimal Confidence Score (Stability Index
#' audit from \code{taxa_retention()}, Step 002) and the optimal minimum reads
#' (\code{Reads_Audit} from \code{reads_per_taxa()}, Step 003). Both the
#' \code{CS_*} and \code{reads_min_*} arguments accept \code{"auto"},
#' \code{"secondary"}, or a manual numeric value per domain. The optimal reads is
#' looked up at each domain's resolved CS, so the mosaic combines both data-driven
#' thresholds automatically. The mosaic always retains \strong{all} taxonomic
#' ranks present in the input (as an MPA profile does); enforces a strict
#' \code{> 0} read filter.
#'
#' @param project_dir Path to the project root.
#' @param tax_level Which rank's optimization audit the \code{"auto"} /
#'   \code{"secondary"} thresholds are read from, i.e. \code{SI_Audit_<tax_level>}
#'   and \code{Reads_Audit_<tax_level>} (\code{NULL}, the default, uses
#'   \code{"Species"}). This does \emph{not} filter the output, which always
#'   contains all ranks.
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
#' @param export Logical (default \code{TRUE}). When \code{TRUE}, one mosaic
#'   per sample is written to \code{<project_dir>/004_final_mosaic/}
#'   (\code{.mpa} files under \code{mpa/}, \code{.tsv} files under
#'   \code{tsv/}) together with a log. When \code{FALSE}, nothing is written
#'   and the mosaic is only returned; note that \code{\link{taxa_resolution}}
#'   and \code{\link{group_upset}} read the exported mosaic by default.
#'
#' @details
#' \code{"auto"} and \code{"secondary"} read the audit files exported by
#' \code{\link{taxa_retention}} and \code{\link{reads_per_taxa}}; if the
#' required file is missing the function stops with an explanation rather
#' than silently applying no threshold. When the audit exists but has no
#' optimal value for a given sample and domain (typically a domain with too
#' few taxa to fit a curve), CS 0 or 0 minimum reads is used for that domain
#' and a single warning lists every such case.
#'
#' @return A \code{\link{kariocas_result}} object. \code{$data} is the mosaic
#'   in long format: one row per sample and taxon, with the read
#'   \code{Counts} at the selected CS, the \code{Domain}, the \code{CS} and
#'   \code{Min_Reads} applied, and where each threshold came from
#'   (\code{CS_Source}, \code{Min_Reads_Source}). \code{$plots} is empty and
#'   \code{$paths} lists the files written when \code{export = TRUE}.
#' @export
#' @importFrom readr read_rds write_tsv write_delim
#' @importFrom dplyr filter mutate select group_by summarise bind_rows
#'   rename left_join case_when
#' @importFrom stringr str_detect str_remove str_extract
#' @importFrom SummarizedExperiment assay rowData
#' @examples
#' # Copy the bundled toy project to a temporary folder and import it
#' toy_project <- file.path(tempdir(), "toy_karioCaS")
#' dir.create(toy_project, showWarnings = FALSE)
#' file.copy(
#'     system.file("extdata", "your_project_name", "000_mpa_original",
#'         package = "karioCaS"
#'     ),
#'     toy_project,
#'     recursive = TRUE
#' )
#' import_karioCaS(toy_project)
#'
#' # "auto" thresholds come from the exported audits of the previous steps
#' taxa_retention(toy_project)
#' reads_per_taxa(toy_project)
#'
#' # Data-driven CS and minimum reads for Bacteria and Archaea,
#' # manual thresholds for Eukaryota and Viruses
#' mosaic <- retrieve_selected_taxa(
#'     toy_project,
#'     CS_B = "auto", reads_min_B = "auto",
#'     CS_A = "auto", reads_min_A = "auto",
#'     CS_E = 40, reads_min_E = 10,
#'     CS_V = 0, reads_min_V = 0
#' )
#' mosaic
#' unique(mosaic$data[, c("Domain", "CS", "CS_Source", "Min_Reads")])
#'
#' unlink(toy_project, recursive = TRUE)
retrieve_selected_taxa <- function(project_dir,
                                   tax_level = NULL,
                                   CS_A = "auto", reads_min_A = 0,
                                   CS_B = "auto", reads_min_B = 0,
                                   CS_E = "auto", reads_min_E = 0,
                                   CS_V = "auto", reads_min_V = 0,
                                   export = TRUE) {
    configs <- list(
        Archaea   = list(val = CS_A, min_val = reads_min_A),
        Bacteria  = list(val = CS_B, min_val = reads_min_B),
        Eukaryota = list(val = CS_E, min_val = reads_min_E),
        Viruses   = list(val = CS_V, min_val = reads_min_V)
    )
    .rst_validate_configs(configs)
    setup <- .rst_setup(project_dir, export)
    log_msg <- setup$log_msg
    withCallingHandlers(
        {
            tse_data <- .rst_load_tse(project_dir, log_msg)
            audit_df <- .rst_load_audit(project_dir, tax_level, configs, log_msg)
            reads_audit <- .rst_load_reads_audit(
                project_dir, tax_level, configs, log_msg
            )
            log_msg("STEP 2: Processing Mosaics...")
            per_sample <- lapply(tse_data$SAMPLES, function(samp) {
                .rst_process_sample(
                    samp, tse_data$count_matrix, tse_data$row_meta,
                    tse_data$all_cols, configs, audit_df, reads_audit, setup
                )
            })
        },
        error = function(e) {
            if (!is.null(setup$log_file)) {
                write(paste0("\nCRITICAL ERROR: ", conditionMessage(e)),
                    file = setup$log_file, append = TRUE
                )
            }
        }
    )
    issues <- unlist(lapply(per_sample, `[[`, "issues"))
    if (length(issues) > 0) {
        log_msg("  [WARNING] Threshold fallbacks: ", paste(issues, collapse = "; "))
        warning(
            "Some thresholds could not be resolved as requested:\n  ",
            paste(issues, collapse = "\n  "),
            call. = FALSE
        )
    }
    .kcs_finish_step(
        "004_final_mosaic",
        dplyr::bind_rows(lapply(per_sample, `[[`, "data")),
        list(),
        unlist(lapply(per_sample, `[[`, "paths")),
        setup, log_msg,
        what = "mosaic taxa"
    )
}
