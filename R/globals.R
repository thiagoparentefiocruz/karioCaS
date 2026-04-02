# R/globals.R

#' @importFrom methods is
#' @importFrom stats setNames
#' @importFrom utils head
#' @importFrom rlang :=
NULL

utils::globalVariables(c(
  ".data", "Base_Global", "Base_Reads", "Base_Taxa", "CS", "CS_Label",
  "Category", "Child_Sum_Resolved", "Class", "Col_ID", "Confidence_Score",
  "Counts", "Cutoff", "Distance", "Domain", "Family", "Genus", "Global_Reads",
  "Is_Candidate", "Is_Local_Max", "Kingdom", "Lowest_Rank", "Metric_Key",
  "Metric_Label", "Norm_CS", "Norm_Taxa", "Order", "Order_Index",
  "Original_CS_Code", "Parent_Cumulative_Total", "Parent_Exclusive",
  "Parent_Name", "Pct", "Pct_Domain_Retained", "Pct_Rank_Reads", "Pct_Retained",
  "Pct_Taxa", "Pct_Value", "Phylum", "Present", "Rank", "Rank_Reads",
  "Rank_Taxa", "Reads", "Rel_Abund", "Ret_Reads_Pct", "Ret_Taxa_Pct",
  "SI_Type", "Sample", "Sample_ID", "Species", "Taxa_Count", "Taxon_Name",
  "Taxonomy", "Taxonomy_Full", "Total_CS_Reads", "Total_Clade_Reads",
  "Wide_Header", "kariocas_colors", "kariocas_dims", "kariocas_labels",
  "kariocas_linetypes"
))
