# ==============================================================================
# PRIVATE HELPERS - group_upset()
# ==============================================================================
# Cross-sample UpSet within each biological group: which taxa are CORE (present
# in every sample of the group) vs UNIQUE/rare (present in one or a few samples)
# - the expected signature of pathogens and false positives. Also writes a
# membership TSV per group/domain.
# ==============================================================================

#' @noRd
.gup_setup <- function(project_dir) {
    output_dir <- file.path(project_dir, "006_group_upset")
    log_dir <- file.path(project_dir, "logs")
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)
    log_file <- file.path(log_dir, "log_006_group_upset.txt")
    writeLines(c(
        "====================================================",
        "LOG: 006_GROUP_UPSET (core vs unique taxa across samples)",
        paste0("PROJECT DIR: ", project_dir),
        "===================================================="
    ), con = log_file)
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
        log_msg(">>> Source: Final Mosaic (1000_final_selection)...")
        mdir <- file.path(project_dir, "1000_final_selection")
        files <- list.files(
            mdir, pattern = "_karioCaS_Mosaic\\.tsv$", full.names = TRUE
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

#' @noRd
.gup_plot <- function(binary, samples, group, dom, tax_level, path, log_msg) {
    grDevices::pdf(
        file = path, width = get_kariocas_dims()$width,
        height = get_kariocas_dims()$height, onefile = FALSE
    )
    tryCatch(
        {
            show(
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
            )
            grid::grid.text(
                label = paste(
                    group, "-", dom, "|", tax_level,
                    "core vs unique across samples"
                ),
                x = 0.5, y = 0.98,
                gp = grid::gpar(fontsize = 14, fontface = "bold")
            )
            log_msg("    -> Generated: ", basename(path))
        },
        error = function(e) {
            log_msg("    ERROR plotting ", basename(path), ": ", e$message)
        }
    )
    grDevices::dev.off()
}

#' @noRd
.gup_process_group <- function(df, group, DOMAINS, tax_level, label,
                               output_dir, log_msg) {
    grp_dir <- file.path(output_dir, group)
    if (!dir.exists(grp_dir)) dir.create(grp_dir, recursive = TRUE)
    membership <- list()
    for (dom in DOMAINS) {
        df_gd <- dplyr::filter(df, .data$Group == group, .data$Domain == dom)
        binary <- .gup_binary(df_gd)
        if (is.null(binary)) {
            log_msg("    Skipping ", dom, ": < 2 samples or no taxa.")
            next
        }
        samples <- setdiff(colnames(binary), "Taxon_Name")
        base <- paste0(group, "_", dom, "_", tax_level, "_", label)
        .gup_plot(
            binary, samples, group, dom, tax_level,
            file.path(grp_dir, paste0(base, "_SampleUpSet.pdf")), log_msg
        )
        memb <- .gup_membership(binary, samples, group, dom, tax_level)
        readr::write_tsv(
            memb, file.path(grp_dir, paste0(base, "_membership.tsv"))
        )
        membership[[dom]] <- memb
    }
    dplyr::bind_rows(membership)
}

# ==============================================================================
# EXPORTED FUNCTION
# ==============================================================================

#' Cross-Sample UpSet: Core vs Unique Taxa per Biological Group (Step 006)
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
#'
#' @return Invisibly returns a \code{data.frame} with the full membership table.
#'   UpSet PDFs and membership TSVs are saved per group to
#'   \code{<project_dir>/006_group_upset/}.
#' @export
#' @importFrom dplyr filter mutate distinct bind_rows left_join case_when
#'   n_distinct
#' @importFrom tidyr pivot_wider
#' @importFrom readr read_tsv write_tsv
#' @importFrom UpSetR upset
#' @importFrom grDevices pdf dev.off
#' @importFrom grid grid.text gpar
#' @examples
#' toy_project <- system.file("extdata", "your_project_name", package = "karioCaS")
#'
#' # Core vs unique species across samples, from the final mosaic
#' # group_upset(project_dir = toy_project)
#'
#' # Compare at a single Confidence Score instead
#' # group_upset(project_dir = toy_project, CS = 40)
group_upset <- function(project_dir, tax_level = "Species", CS = NULL) {
    if (!requireNamespace("UpSetR", quietly = TRUE)) {
        stop("Package 'UpSetR' is required.")
    }
    setup <- .gup_setup(project_dir)
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
    membership <- list()
    for (grp in GROUPS) {
        n_s <- dplyr::n_distinct(df$sample[df$Group == grp])
        setup$log_msg("  Group: ", grp, " (", n_s, " samples)")
        if (n_s < 2) {
            setup$log_msg("    Skipping ", grp, ": needs >= 2 samples.")
            next
        }
        membership[[grp]] <- .gup_process_group(
            df, grp, DOMAINS, tax_level, loaded$label,
            setup$output_dir, setup$log_msg
        )
    }
    setup$log_msg("SUCCESS: Group UpSet analysis completed.")
    invisible(dplyr::bind_rows(membership))
}
