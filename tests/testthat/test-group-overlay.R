# tests/testthat/test-group-overlay.R
# Tests for the shared group-overlay engine and the detail_samples opt-in
# wired into taxa_retention() (001) and reads_per_taxa() (003).

test_that(".grp_parse_group strips trailing digits", {
    g <- karioCaS:::.grp_parse_group
    expect_equal(g(c("SAMPLE33", "SAMPLE34")), c("SAMPLE", "SAMPLE"))
    expect_equal(g(c("CONTROL01", "TREATED02")), c("CONTROL", "TREATED"))
    expect_equal(g("PILO"), "PILO")
})

test_that(".grp_resolve_detail handles NULL, 'all', strings and vectors", {
    r <- karioCaS:::.grp_resolve_detail
    noop <- function(...) invisible(NULL)
    all_s <- c("SAMPLE33", "SAMPLE45", "SAMPLE50")
    expect_equal(r(NULL, all_s, noop), character(0))
    expect_equal(r("all", all_s, noop), all_s)
    expect_equal(r("SAMPLE33, SAMPLE45", all_s, noop), c("SAMPLE33", "SAMPLE45"))
    expect_equal(r(c("SAMPLE50"), all_s, noop), "SAMPLE50")
    # Unknown names are dropped
    expect_equal(r("SAMPLE33, NOPE", all_s, noop), "SAMPLE33")
})

test_that("taxa_retention writes a group overlay by default and detail on request", {
    temp_proj_dir <- tempfile(pattern = "kariocas_test_grp_")
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

    out_dir <- file.path(temp_proj_dir, "001_taxa_retention")

    # Default: group overlay only, no per_sample folder
    suppressMessages(taxa_retention(project_dir = temp_proj_dir))
    expect_true(file.exists(
        file.path(out_dir, "SAMPLE_Group_Retention_Species.pdf")
    ))
    expect_false(dir.exists(file.path(out_dir, "per_sample")))

    # Opt-in detail for a named sample -> per_sample subfolder
    suppressMessages(
        taxa_retention(project_dir = temp_proj_dir, detail_samples = "SAMPLE01")
    )
    detail_pdfs <- list.files(file.path(out_dir, "per_sample"), pattern = "\\.pdf$")
    expect_true(length(detail_pdfs) > 0)
    expect_true(all(grepl("^SAMPLE01", detail_pdfs)))

    unlink(temp_proj_dir, recursive = TRUE)
})
