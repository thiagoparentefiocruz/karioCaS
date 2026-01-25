#' Utilities for karioCaS
#'
#' Internal helper functions for data processing.
#' @keywords internal
#' @noRd
#' @importFrom readr read_rds
#' @importFrom tidyr pivot_longer extract
#' @importFrom dplyr rename mutate filter case_when
#' @importFrom stringr str_extract str_remove str_detect str_sub

load_and_process_mpa_archive <- function(rds_file) {

  if (!file.exists(rds_file)) {
    stop("Data file not found: ", rds_file)
  }

  df_raw <- readr::read_rds(rds_file)

  # Check if it is already in Long format (has Sample column)
  if ("Sample" %in% names(df_raw)) {
    return(df_raw)
  }

  # If Wide format (Step 000 output), process it to Long
  # Expected columns: Taxonomy, Tax_Level, SAMPLE01_CS00, SAMPLE01_CS02...

  message("Converting Wide Matrix to Long Format...")

  df_long <- df_raw %>%
    # 1. Pivot Samples to Rows
    tidyr::pivot_longer(
      cols = -c("Taxonomy", "Tax_Level"),
      names_to = "Sample_CS",
      values_to = "Counts"
    ) %>%
    # 2. Extract Sample and CS from the column name
    # Regex logic: Everything before the last '_CS' is Sample, everything after is CS
    tidyr::extract(
      col = "Sample_CS",
      into = c("Sample", "CS"),
      regex = "^(.*)_(CS[0-9\\.]+)$"
    ) %>%
    # 3. Rename Level
    dplyr::rename(Level = Tax_Level) %>%
    # 4. Clean Taxonomy and Normalize CS
    dplyr::mutate(
      # Extract Domain
      Domain = stringr::str_extract(Taxonomy, "d__[A-Za-z]+"),
      Domain = stringr::str_remove(Domain, "d__"),

      # Extract Taxon Name
      Taxon = stringr::str_extract(Taxonomy, "[^|]+$"),
      Taxon = stringr::str_remove(Taxon, "^[kpcofgs]__"),

      # --- CORREÇÃO DE FORMATO DO CS ---
      # Se o CS for "CS" seguido de 2 digitos (ex: CS00, CS02), insere o ponto (CS0.0, CS0.2)
      # Isso garante compatibilidade com os scripts de análise
      CS = dplyr::case_when(
        stringr::str_detect(CS, "^CS[0-9]{2}$") ~ paste0("CS", stringr::str_sub(CS, 3, 3), ".", stringr::str_sub(CS, 4, 4)),
        TRUE ~ CS # Se já tiver ponto (CS0.0), mantém igual
      )
    )

  return(df_long)
}
