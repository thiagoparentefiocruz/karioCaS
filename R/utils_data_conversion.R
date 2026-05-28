#' @noRd
.gtd_load_tse <- function(input) {
  if (is(input, "TreeSummarizedExperiment") ||
      is(input, "SummarizedExperiment")) {
    return(input)
  }
  if (is.character(input) && file.exists(input) && !dir.exists(input)) {
    obj <- readRDS(input)
    if (is(obj, "TreeSummarizedExperiment") ||
        is(obj, "SummarizedExperiment")) {
      return(obj)
    }
    if (is.data.frame(obj)) return(obj)
  }
  if (is.character(input) && dir.exists(input)) {
    tse_path <- file.path(
      input, "000_karioCaS_input_matrix", "karioCaS_TSE.rds"
    )
    rds_path <- file.path(
      input, "000_karioCaS_input_matrix", "karioCaS_input_matrix.rds"
    )
    if (file.exists(tse_path)) return(readRDS(tse_path))
    if (file.exists(rds_path)) return(readRDS(rds_path))
    stop("No valid karioCaS data found in directory: ", input)
  }
  stop("Invalid input data.")
}

#' @noRd
.gtd_extract_tables <- function(tse_obj) {
  counts_df <- as.data.frame(
    SummarizedExperiment::assay(tse_obj, "counts")
  ) |>
    tibble::rownames_to_column("Taxonomy_Full") |>
    tidyr::pivot_longer(
      cols      = -"Taxonomy_Full",
      names_to  = "Col_ID",
      values_to = "Counts"
    )
  meta_df <- as.data.frame(
    SummarizedExperiment::colData(tse_obj)
  ) |>
    tibble::rownames_to_column("Col_ID")
  tax_df <- as.data.frame(SummarizedExperiment::rowData(tse_obj))
  if (!"Taxonomy_Full" %in% colnames(tax_df)) {
    tax_df$Taxonomy_Full <- rownames(tax_df)
  }
  list(counts_df = counts_df, meta_df = meta_df, tax_df = tax_df)
}

#' @noRd
.gtd_merge_and_pivot <- function(counts_df, meta_df, tax_df) {
  rank_cols <- c(
    "Domain", "Kingdom", "Phylum", "Class",
    "Order", "Family", "Genus", "Species"
  )
  df_base <- counts_df |>
    dplyr::left_join(meta_df, by = "Col_ID") |>
    dplyr::left_join(tax_df,  by = "Taxonomy_Full") |>
    dplyr::rename(sample = "Sample_ID", CS = "Confidence_Score") |>
    dplyr::filter(.data$Counts > 0)
  if ("Rank" %in% colnames(df_base)) {
    df_base <- dplyr::rename(df_base, Lowest_Rank = "Rank")
  }
  if ("Domain"  %in% colnames(df_base)) df_base$Domain_Group  <- df_base$Domain
  if ("Kingdom" %in% colnames(df_base)) df_base$Kingdom_Group <- df_base$Kingdom
  df_long <- df_base |>
    tidyr::pivot_longer(
      cols      = dplyr::all_of(rank_cols),
      names_to  = "Rank",
      values_to = "Taxon_Name"
    ) |>
    dplyr::filter(!is.na(.data$Taxon_Name))
  if ("Domain_Group"  %in% colnames(df_long)) df_long$Domain  <- df_long$Domain_Group
  if ("Kingdom_Group" %in% colnames(df_long)) df_long$Kingdom <- df_long$Kingdom_Group
  cols_to_keep <- c(
    "sample", "CS", "Domain", "Kingdom", "Rank",
    "Taxon_Name", "Counts", "Taxonomy_Full", "Lowest_Rank"
  )
  dplyr::select(df_long, dplyr::all_of(intersect(cols_to_keep, colnames(df_long))))
}

#' Internal: Convert TSE to Tidy Dataframe
#'
#' Helper to bridge Bioconductor TSE objects with karioCaS plotting functions.
#' Handles column name collisions and preserves Taxonomy groupings.
#'
#' @param input A path to an RDS file, a project directory, or a TSE object.
#' @return A tidy \code{data.frame} compatible with karioCaS functions.
#' @noRd
#' @importFrom SummarizedExperiment assay rowData colData
#' @importFrom dplyr left_join mutate select filter rename all_of
#' @importFrom tidyr pivot_longer
#' @importFrom tibble rownames_to_column
.get_tidy_data <- function(input) {
  result <- .gtd_load_tse(input)
  if (is.data.frame(result)) return(result)
  tables <- .gtd_extract_tables(result)
  .gtd_merge_and_pivot(tables$counts_df, tables$meta_df, tables$tax_df)
}