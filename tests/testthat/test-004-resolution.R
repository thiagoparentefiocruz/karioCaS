# tests/testthat/test-004-resolution.R
# taxa_resolution(): default = final mosaic; CS = single Confidence Score.

.kcs_setup_resolution_proj <- function() {
    temp_proj_dir <- tempfile(pattern = "kariocas_test_res_")
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
    temp_proj_dir
}

test_that("taxa_resolution(CS = ...) analyses a single Confidence Score", {
    proj <- .kcs_setup_resolution_proj()
    out_dir <- file.path(proj, "004_taxa_resolution")

    expect_message(
        taxa_resolution(project_dir = proj, CS = 40),
        "SUCCESS: Resolution analysis completed."
    )
    pdfs <- list.files(out_dir, pattern = "\\.pdf$")
    expect_true(length(pdfs) > 0)
    expect_true(all(grepl("_CS40\\.pdf$", pdfs)))
    # Single CS only -- not one per CS level
    expect_false(any(grepl("_CS00\\.pdf$", pdfs)))

    # Invalid CS is rejected
    expect_error(taxa_resolution(project_dir = proj, CS = 999))

    unlink(proj, recursive = TRUE)
})

test_that("taxa_resolution default reads the final mosaic", {
    proj <- .kcs_setup_resolution_proj()
    # Mosaic with all ranks so the parent rank is present.
    suppressMessages(taxa_retention(project_dir = proj, tax_level = "Species"))
    suppressMessages(retrieve_selected_taxa(
        project_dir = proj, tax_level = NULL,
        CS_A = "auto", CS_B = "auto", CS_E = "auto", CS_V = "auto"
    ))

    expect_message(
        taxa_resolution(project_dir = proj),
        "SUCCESS: Resolution analysis completed."
    )
    pdfs <- list.files(
        file.path(proj, "004_taxa_resolution"),
        pattern = "\\.pdf$"
    )
    expect_true(any(grepl("Final_Mosaic\\.pdf$", pdfs)))

    unlink(proj, recursive = TRUE)
})

test_that("taxa_resolution errors clearly when no mosaic exists", {
    proj <- .kcs_setup_resolution_proj()
    expect_error(
        suppressMessages(taxa_resolution(project_dir = proj)),
        "No mosaic files"
    )
    unlink(proj, recursive = TRUE)
})
