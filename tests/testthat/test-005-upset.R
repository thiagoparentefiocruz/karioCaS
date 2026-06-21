# tests/testthat/test-002-upset.R
# upset_kariocas(): single tax_level flag (default "Species").

test_that("upset_kariocas draws one rank per sample/domain and validates rank", {
    temp_proj_dir <- tempfile(pattern = "kariocas_test_ups_")
    dir.create(file.path(temp_proj_dir, "000_mpa_original"), recursive = TRUE)
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
    file.copy(
        list.files(mock_data_src, full.names = TRUE),
        file.path(temp_proj_dir, "000_mpa_original")
    )
    suppressMessages(import_karioCaS(project_dir = temp_proj_dir))

    expect_message(
        upset_kariocas(project_dir = temp_proj_dir),
        "SUCCESS: UpSet analysis completed."
    )
    samp_dir <- file.path(
        temp_proj_dir, "005_taxa_intersections_across_CS", "SAMPLE01"
    )
    pdfs <- list.files(samp_dir, pattern = "\\.pdf$")
    expect_true(length(pdfs) > 0)
    # Default rank only -> Species files, no Genus/Family
    expect_true(all(grepl("_Species_UpSet\\.pdf$", pdfs)))
    expect_false(any(grepl("_Genus_UpSet\\.pdf$", pdfs)))

    # Invalid rank rejected with a helpful message
    expect_error(
        suppressMessages(upset_kariocas(temp_proj_dir, tax_level = "Nope")),
        "not found"
    )

    unlink(temp_proj_dir, recursive = TRUE)
})
