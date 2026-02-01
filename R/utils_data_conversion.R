#' Internal: Convert TSE to Tidy Dataframe
#'
#' Helper function to bridge Bioconductor objects (TSE) with karioCaS plotting functions.
#' Handles column name collisions (Rank vs Lowest_Rank) and preserves Taxonomy.
#'
#' @param input Input which can be a path to RDS, a directory, or a TSE object.
#' @return A tidy data.frame compatible with karioCaS functions.
#' @noRd
#' @importFrom SummarizedExperiment assay rowData colData
#' @importFrom dplyr left_join mutate select filter rename
#' @importFrom tidyr pivot_longer
#' @importFrom tibble rownames_to_column

.get_tidy_data <- function(input) {

  tse_obj <- NULL

  # 1. Handle Input Types
  if (is(input, "TreeSummarizedExperiment") || is(input, "SummarizedExperiment")) {
    tse_obj <- input
  } else if (is.character(input) && file.exists(input) && !dir.exists(input)) {
    obj <- readRDS(input)
    if (is(obj, "TreeSummarizedExperiment") || is(obj, "SummarizedExperiment")) {
      tse_obj <- obj
    } else if (is.data.frame(obj)) {
      return(obj)
    }
  } else if (is.character(input) && dir.exists(input)) {
    tse_path <- file.path(input, "000_karioCaS_input_matrix", "karioCaS_TSE.rds")
    rds_path <- file.path(input, "000_karioCaS_input_matrix", "karioCaS_input_matrix.rds")

    if (file.exists(tse_path)) {
      tse_obj <- readRDS(tse_path)
    } else if (file.exists(rds_path)) {
      return(readRDS(rds_path))
    } else {
      stop("No valid karioCaS data found in directory: ", input)
    }
  }

  if (is.null(tse_obj)) stop("Invalid input data.")

  # 2. Extract Data
  counts_df <- as.data.frame(SummarizedExperiment::assay(tse_obj, "counts")) %>%
    tibble::rownames_to_column("Taxonomy_Full") %>%
    tidyr::pivot_longer(cols = -Taxonomy_Full, names_to = "Col_ID", values_to = "Counts")

  meta_df <- as.data.frame(SummarizedExperiment::colData(tse_obj)) %>%
    tibble::rownames_to_column("Col_ID")

  tax_df <- as.data.frame(SummarizedExperiment::rowData(tse_obj))
  if (!"Taxonomy_Full" %in% colnames(tax_df)) {
    tax_df$Taxonomy_Full <- rownames(tax_df)
  }

  # 3. Merge to Base Tidy DF
  df_base <- counts_df %>%
    dplyr::left_join(meta_df, by = "Col_ID") %>%
    dplyr::left_join(tax_df, by = "Taxonomy_Full") %>%
    dplyr::rename(sample = Sample_ID, CS = Confidence_Score) %>%
    dplyr::filter(Counts > 0)

  # CRITICAL FIX: Rename existing 'Rank' from rowData to 'Lowest_Rank'
  # This prevents collision with the new 'Rank' column created by pivot_longer below
  if ("Rank" %in% colnames(df_base)) {
    df_base <- df_base %>% dplyr::rename(Lowest_Rank = Rank)
  }

  # 4. Prepare for Plotting (Pivot Ranks)

  # Create static copies of high-level groups
  # (Because pivot_longer will melt the original columns)
  if ("Domain" %in% colnames(df_base)) df_base$Domain_Group  <- df_base$Domain
  if ("Kingdom" %in% colnames(df_base)) df_base$Kingdom_Group <- df_base$Kingdom

  # Melt ranks into long format
  # This creates the "Rank" column used by plotting functions (x-axis)
  df_long <- df_base %>%
    tidyr::pivot_longer(
      cols = c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
      names_to = "Rank",
      values_to = "Taxon_Name"
    ) %>%
    dplyr::filter(!is.na(Taxon_Name))

  # Restore grouping columns
  if ("Domain_Group" %in% colnames(df_long)) df_long$Domain <- df_long$Domain_Group
  if ("Kingdom_Group" %in% colnames(df_long)) df_long$Kingdom <- df_long$Kingdom_Group

  # Final Selection
  # We include Lowest_Rank in the output as it's useful info
  cols_to_keep <- c("sample", "CS", "Domain", "Kingdom", "Rank", "Taxon_Name", "Counts", "Taxonomy_Full", "Lowest_Rank")
  # Intersect with existing names to be safe (in case Lowest_Rank didn't exist)
  cols_final <- intersect(cols_to_keep, colnames(df_long))

  df_final <- df_long %>%
    dplyr::select(dplyr::all_of(cols_final))

  return(df_final)
}
