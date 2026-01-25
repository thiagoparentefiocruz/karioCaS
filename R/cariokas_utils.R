#' Utilities for karioCaS
#'
#' Internal helper functions for data processing.
#' Converts MPA format to a fully parsed, tidy long format with individual taxonomic columns.
#'
#' @keywords internal
#' @noRd
#' @importFrom readr read_rds
#' @importFrom tidyr pivot_longer extract
#' @importFrom dplyr rename mutate select case_when
#' @importFrom stringr str_extract str_remove str_detect str_sub

load_and_process_mpa <- function(rds_file) {

  if (!file.exists(rds_file)) {
    stop("Data file not found: ", rds_file)
  }

  df_raw <- readr::read_rds(rds_file)

  # Check if it is already in the new Tidy format (checking a specific column like 'Phylum')
  if ("Phylum" %in% names(df_raw) && "CS" %in% names(df_raw)) {
    return(df_raw)
  }

  message("Converting Wide Matrix to Enhanced Tidy Long Format...")

  # Step 1: Pivot Wide to Long
  # Expected raw columns: Taxonomy, Tax_Level, SAMPLE01_CS00, SAMPLE01_CS02...
  df_long <- df_raw %>%
    tidyr::pivot_longer(
      cols = -c("Taxonomy", "Tax_Level"),
      names_to = "Sample_CS_Raw",
      values_to = "Counts"
    ) %>%
    # Step 2: Extract Sample and Raw CS String
    tidyr::extract(
      col = "Sample_CS_Raw",
      into = c("sample", "CS_String"),
      regex = "^(.*)_(CS[0-9\\.]+)$"
    )

  # Step 3: Parsing & Cleaning
  df_tidy <- df_long %>%
    dplyr::mutate(
      # --- TAXONOMY PARSING ---
      # Extracts content after prefix (e.g., 'p__') but before the next pipe '|' or end of string
      Domain  = stringr::str_extract(Taxonomy, "d__[^|]+") %>% stringr::str_remove("d__"),
      Kingdom = stringr::str_extract(Taxonomy, "k__[^|]+") %>% stringr::str_remove("k__"),
      Phylum  = stringr::str_extract(Taxonomy, "p__[^|]+") %>% stringr::str_remove("p__"),
      Class   = stringr::str_extract(Taxonomy, "c__[^|]+") %>% stringr::str_remove("c__"),
      Order   = stringr::str_extract(Taxonomy, "o__[^|]+") %>% stringr::str_remove("o__"),
      Family  = stringr::str_extract(Taxonomy, "f__[^|]+") %>% stringr::str_remove("f__"),
      Genus   = stringr::str_extract(Taxonomy, "g__[^|]+") %>% stringr::str_remove("g__"),
      Species = stringr::str_extract(Taxonomy, "s__[^|]+") %>% stringr::str_remove("s__"),

      # --- CS CLEANING (Percentage Integer) ---
      # 1. Remove "CS" prefix
      CS_Temp = stringr::str_remove(CS_String, "CS"),
      # 2. Normalize "00", "02" to "0.0", "0.2" if dot is missing
      CS_Temp = dplyr::case_when(
        !stringr::str_detect(CS_Temp, "\\.") & stringr::str_length(CS_Temp) == 2 ~
          paste0(stringr::str_sub(CS_Temp, 1, 1), ".", stringr::str_sub(CS_Temp, 2, 2)),
        TRUE ~ CS_Temp
      ),
      # 3. Convert to Numeric Percentage (0.2 -> 20)
      CS = as.numeric(CS_Temp) * 100
    ) %>%
    # Rename Level to Rank (Standardization)
    dplyr::rename(Rank = Tax_Level) %>%
    # Organize Columns exactly as requested
    dplyr::select(
      Taxonomy,
      Domain, Kingdom, Phylum, Class, Order, Family, Genus, Species,
      Rank,
      sample, # lowercase
      CS,     # numeric percentage
      Counts
    )

  return(df_tidy)
}
