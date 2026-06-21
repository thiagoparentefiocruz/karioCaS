# tests/testthat/test-001-si.R
# taxa_retention() computes the Stability Index audit (formerly optimize_CS)
# and writes it into the 002_taxa_retention/ folder.

test_that("taxa_retention computes the SI audit and writes outputs", {
    temp_proj_dir <- tempfile(pattern = "kariocas_test_si_")
    dir.create(temp_proj_dir)

    mock_data_src <- system.file(
        "extdata/your_project_name/000_mpa_original",
        package = "karioCaS"
    )
    if (mock_data_src == "") {
        mock_data_src <- file.path(
            "..", "..", "inst", "extdata",
            "your_project_name", "000_mpa_original"
        )
    }
    dir.create(file.path(temp_proj_dir, "000_mpa_original"), recursive = TRUE)
    file.copy(
        list.files(mock_data_src, full.names = TRUE),
        file.path(temp_proj_dir, "000_mpa_original")
    )
    suppressMessages(import_karioCaS(project_dir = temp_proj_dir))

    expect_message(
        audit_df <- taxa_retention(
            project_dir = temp_proj_dir, tax_level = "Species"
        ),
        "SUCCESS: Retention analysis completed."
    )

    # Returns the SI audit data frame with a Primary SI
    expect_s3_class(audit_df, "data.frame")
    expect_true("Primary_SI" %in% audit_df$SI_Type)

    # Audit files + group plot land in the 001 folder
    out_dir <- file.path(temp_proj_dir, "002_taxa_retention")
    expect_true(file.exists(file.path(out_dir, "SI_Audit_Species.rds")))
    expect_true(file.exists(file.path(out_dir, "SI_Audit_Species.tsv")))
    expect_true(length(list.files(out_dir, pattern = "\\.pdf$")) > 0)

    unlink(temp_proj_dir, recursive = TRUE)
})
