# ==============================================================================
# PRIVATE HELPERS - group_upset()
# ==============================================================================
# Cross-sample UpSet within each biological group: which taxa are CORE (present
# in every sample of the group) vs UNIQUE/rare (present in one or a few samples)
# - the expected signature of pathogens and false positives. Also writes a
# membership TSV per group/domain.
# ==============================================================================

#' @noRd
.gup_setup <- function(project_dir, export = TRUE) {
    .kcs_setup_step(
        project_dir, "008_taxa_intersections_across_samples",
        "log_008_taxa_intersections_samples.txt",
        export = export
    )
}

#' Load presence data (sample x taxon) from the mosaic or a single CS.
#' @return list(df = distinct Group/sample/Domain/Taxon_Name, label).
#' @noRd
.gup_load <- function(project_dir, tax_level, CS, log_msg) {
    rank_cols <- c(
        "Domain", "Kingdom", "Phylum", "Class",
        "Order", "Family", "Genus", "Species"
    )
    if (!tax_level %in% rank_cols) {
        stop("Invalid 'tax_level': ", tax_level)
    }
    if (is.null(CS)) {
        log_msg(">>> Source: Final Mosaic (004_final_mosaic)...")
        mdir <- .kcs_path(project_dir, "004_final_mosaic", "1000_final_selection")
        tsv_dir <- if (dir.exists(file.path(mdir, "tsv"))) {
            file.path(mdir, "tsv")
        } else {
            mdir
        }
        files <- list.files(
            tsv_dir,
            pattern = "_karioCaS_Mosaic\\.tsv$", full.names = TRUE
        )
        if (length(files) == 0) {
            stop(
                "No mosaic files found in ", mdir,
                ". Run retrieve_selected_taxa() first."
            )
        }
        mosaic <- dplyr::bind_rows(lapply(files, function(f) {
            samp <- sub("_karioCaS_Mosaic\\.tsv$", "", basename(f))
            d <- readr::read_tsv(f, show_col_types = FALSE, progress = FALSE)
            data.frame(
                sample = samp, Taxonomy_Full = d[[1]],
                Counts = as.numeric(d[[2]]), stringsAsFactors = FALSE
            )
        }))
        tax_df <- .imp_parse_taxonomy(unique(mosaic$Taxonomy_Full), log_msg)
        df <- mosaic |>
            dplyr::left_join(tax_df, by = "Taxonomy_Full") |>
            dplyr::filter(.data$Rank == tax_level) |>
            dplyr::mutate(Taxon_Name = .data[[tax_level]])
        label <- "Final_Mosaic"
    } else {
        cs_pct <- .cs_arg_to_percent(CS)
        if (is.na(cs_pct)) {
            stop("Invalid 'CS': ", CS, ". Use a fraction (0-1) or percent (0-100).")
        }
        df_long <- .get_tidy_data(project_dir)
        avail <- sort(unique(df_long$CS))
        if (!cs_pct %in% avail) {
            stop(
                "CS ", cs_pct, "% not found. Available: ",
                paste(avail, collapse = ", ")
            )
        }
        df <- dplyr::filter(
            df_long, .data$CS == cs_pct, .data$Rank == tax_level
        )
        label <- paste0("CS", sprintf("%02d", cs_pct))
    }
    df <- dplyr::filter(df, .data$Counts > 0, !is.na(.data$Taxon_Name))
    df$Group <- .grp_parse_group(df$sample)
    list(
        df = dplyr::distinct(
            df, .data$Group, .data$sample, .data$Domain, .data$Taxon_Name
        ),
        label = label
    )
}

#' Presence/absence matrix (taxa x samples) for one group + domain.
#' @return A data.frame (Taxon_Name + one 0/1 column per sample), or NULL if
#'   fewer than 2 samples.
#' @noRd
.gup_binary <- function(df_gd) {
    samples <- unique(df_gd$sample)
    if (length(samples) < 2 || nrow(df_gd) == 0) {
        return(NULL)
    }
    df_gd |>
        dplyr::distinct(.data$Taxon_Name, .data$sample) |>
        dplyr::mutate(Present = 1L) |>
        tidyr::pivot_wider(
            names_from = "sample", values_from = "Present", values_fill = 0L
        ) |>
        as.data.frame()
}

#' Membership table: per taxon, how many samples and Core/Shared/Unique.
#' @noRd
.gup_membership <- function(binary, samples, group, dom, tax_level) {
    m <- as.matrix(binary[, samples, drop = FALSE])
    n <- rowSums(m)
    total <- length(samples)
    category <- dplyr::case_when(
        n == total ~ "Core",
        n == 1 ~ "Unique",
        TRUE ~ "Shared"
    )
    unique_sample <- ifelse(
        n == 1, samples[max.col(m, ties.method = "first")], NA_character_
    )
    out <- data.frame(
        Group = group, Domain = dom, Rank = tax_level,
        Taxon = binary$Taxon_Name, N_Samples = n,
        Category = category, Unique_Sample = unique_sample,
        stringsAsFactors = FALSE
    )
    out <- cbind(out, binary[, samples, drop = FALSE])
    out[order(-out$N_Samples, out$Taxon), ]
}

#' Build the cross-sample UpSet object (does not draw it)
#' @noRd
.gup_build <- function(binary, samples, tax_level) {
    UpSetR::upset(
        binary,
        sets = samples,
        nsets = length(samples),
        nintersects = 40,
        order.by = "freq",
        mainbar.y.label = paste(tax_level, "shared across samples"),
        sets.x.label = paste("Total", tax_level, "per sample"),
        main.bar.color = get_kariocas_colors("upset")$main,
        sets.bar.color = get_kariocas_colors("upset")$sets,
        matrix.color = "grey20", shade.color = "grey90"
    )
}

#' @noRd
.gup_process_group <- function(df, group, DOMAINS, tax_level, label, setup) {
    log_msg <- setup$log_msg
    grp_dir <- if (isTRUE(setup$export)) {
        file.path(setup$output_dir, group)
    } else {
        NULL
    }
    membership <- list()
    plots <- list()
    paths <- character(0)
    failed <- character(0)
    for (dom in DOMAINS) {
        df_gd <- dplyr::filter(df, .data$Group == group, .data$Domain == dom)
        binary <- .gup_binary(df_gd)
        if (is.null(binary)) {
            log_msg("    Skipping ", dom, ": < 2 samples or no taxa.")
            next
        }
        samples <- setdiff(colnames(binary), "Taxon_Name")
        base <- paste0(group, "_", dom, "_", tax_level, "_", label)
        memb <- .gup_membership(binary, samples, group, dom, tax_level)
        membership[[dom]] <- memb
        up <- tryCatch(
            .gup_build(binary, samples, tax_level),
            error = function(e) {
                failed <<- c(failed, paste0(base, ": ", conditionMessage(e)))
                NULL
            }
        )
        if (!is.null(grp_dir)) {
            if (!dir.exists(grp_dir)) dir.create(grp_dir, recursive = TRUE)
            tsv <- file.path(grp_dir, paste0(base, "_membership.tsv"))
            readr::write_tsv(memb, tsv)
            paths <- c(paths, tsv)
        }
        if (is.null(up)) next
        if (!is.null(grp_dir)) {
            path <- .kcs_draw_upset_pdf(
                up, file.path(grp_dir, paste0(base, "_SampleUpSet.pdf")),
                paste(
                    group, "-", dom, "|", tax_level,
                    "core vs unique across samples"
                ),
                log_msg
            )
            if (!file.exists(path)) {
                failed <- c(failed, paste0(base, ": ", attr(path, "error")))
                next
            }
            paths <- c(paths, path)
        }
        plots[[base]] <- up
    }
    list(
        data = dplyr::bind_rows(membership), plots = plots,
        paths = paths, failed = failed
    )
}

# ==============================================================================
# EXPORTED FUNCTION
# ==============================================================================

#' Cross-Sample UpSet: Core vs Unique Taxa per Biological Group (Step 008)
#'
#' For each biological group (inferred from sample names by stripping trailing
#' digits, e.g. \code{SAMPLE33}, \code{SAMPLE34} -> \code{SAMPLE}), draws an
#' UpSet plot comparing which taxa (at \code{tax_level}) are present across the
#' samples of the group, per Domain. This separates the \strong{core} taxa
#' (present in every sample) from \strong{unique}/rare taxa (present in one or a
#' few samples) - the expected pattern for pathogens and false positives. A
#' membership TSV (presence matrix plus \code{N_Samples} and a
#' Core/Shared/Unique \code{Category}) is written alongside each plot.
#'
#' By default the analysis uses the \strong{final mosaic} from
#' \code{retrieve_selected_taxa()}; set \code{CS} to compare at a single
#' Confidence Score from the imported data instead.
#'
#' @param project_dir Path to the project root.
#' @param tax_level Taxonomic rank to compare (default: \code{"Species"}).
#' @param CS Confidence Score to analyse. \code{NULL} (default) uses the final
#'   mosaic; a numeric value (Kraken fraction \code{0-1} or percentage
#'   \code{0-100}) compares the imported data at that single CS.
#' @param export Logical (default \code{TRUE}). When \code{TRUE}, an UpSet PDF
#'   and a membership TSV per group and domain are written to
#'   \code{<project_dir>/008_taxa_intersections_across_samples/<group>/},
#'   together with a log. When \code{FALSE}, nothing is written and the
#'   results are only returned.
#'
#' @details
#' Groups with a single sample cannot be compared and are skipped with a log
#' message; if no group has at least two samples, the function warns that
#' nothing was generated. \code{Category} is \code{"Core"} when a taxon is
#' present in every sample of the group, \code{"Unique"} when it is present
#' in exactly one, and \code{"Shared"} otherwise.
#'
#' @return A \code{\link{kariocas_result}} object. \code{$data} is the
#'   membership table: one row per group, domain and taxon, with
#'   \code{N_Samples}, \code{Category}, \code{Unique_Sample} and one 0/1
#'   column per sample. \code{$plots} holds one \code{UpSetR} plot per group
#'   and domain (printing an element draws it) and \code{$paths} lists the
#'   files written when \code{export = TRUE}.
#' @export
#' @importFrom dplyr filter mutate distinct bind_rows left_join case_when
#'   n_distinct
#' @importFrom tidyr pivot_wider
#' @importFrom readr read_tsv write_tsv
#' @importFrom UpSetR upset
#' @importFrom grDevices pdf dev.off
#' @importFrom grid grid.text gpar
#' @examples
#' # The toy project has a single sample; add a second sample of the same
#' # group (a copy of SAMPLE01 named SAMPLE02) so they can be compared.
#' toy_project <- file.path(tempdir(), "toy_karioCaS_groups")
#' in_dir <- file.path(toy_project, "000_mpa_original")
#' dir.create(in_dir, recursive = TRUE, showWarnings = FALSE)
#' src <- list.files(
#'     system.file("extdata", "your_project_name", "000_mpa_original",
#'         package = "karioCaS"
#'     ),
#'     full.names = TRUE
#' )
#' file.copy(src, in_dir)
#' file.copy(src, file.path(in_dir, sub("SAMPLE01", "SAMPLE02", basename(src))))
#' import_karioCaS(toy_project)
#'
#' # Core vs unique genera across the samples of each group, at CS 40
#' res <- group_upset(toy_project, tax_level = "Genus", CS = 40, export = FALSE)
#' res
#' table(res$data$Domain, res$data$Category)
#'
#' unlink(toy_project, recursive = TRUE)
group_upset <- function(project_dir, tax_level = "Species", CS = NULL,
                        export = TRUE) {
    setup <- .gup_setup(project_dir, export)
    loaded <- .gup_load(project_dir, tax_level, CS, setup$log_msg)
    df <- loaded$df
    if (nrow(df) == 0) {
        stop("No data found for rank '", tax_level, "'.")
    }
    DOMAINS <- names(get_kariocas_colors("domains"))
    GROUPS <- unique(df$Group)
    setup$log_msg(
        ">>> Cross-sample UpSet (", tax_level, " | ", loaded$label,
        ") for ", length(GROUPS), " group(s)."
    )
    results <- list()
    for (grp in GROUPS) {
        n_s <- dplyr::n_distinct(df$sample[df$Group == grp])
        setup$log_msg("  Group: ", grp, " (", n_s, " samples)")
        if (n_s < 2) {
            setup$log_msg("    Skipping ", grp, ": needs >= 2 samples.")
            next
        }
        results[[grp]] <- .gup_process_group(
            df, grp, DOMAINS, tax_level, loaded$label, setup
        )
    }
    .kcs_warn_failed_plots(unlist(lapply(results, `[[`, "failed")), setup$log_msg)
    plots <- do.call(c, unname(lapply(results, `[[`, "plots")))
    .kcs_finish_step(
        "008_taxa_intersections_across_samples",
        dplyr::bind_rows(lapply(results, `[[`, "data")),
        if (is.null(plots)) list() else plots,
        unlist(lapply(results, `[[`, "paths")),
        setup, setup$log_msg,
        what = "group comparisons (each group needs at least 2 samples)"
    )
}
