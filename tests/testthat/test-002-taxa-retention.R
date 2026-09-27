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
        res <- taxa_retention(
            project_dir = temp_proj_dir, tax_level = "Species"
        ),
        "SUCCESS: 002_taxa_retention completed"
    )

    # Returns a kariocas_result holding the SI audit with a Primary SI
    expect_s3_class(res, "kariocas_result")
    audit_df <- res$data
    expect_s3_class(audit_df, "data.frame")
    expect_true("Primary_SI" %in% audit_df$SI_Type)
    expect_true(length(res$plots) > 0)
    expect_true(all(file.exists(res$paths)))

    # Audit files + group plot land in the 001 folder
    out_dir <- file.path(temp_proj_dir, "002_taxa_retention")
    expect_true(file.exists(file.path(out_dir, "SI_Audit_Species.rds")))
    expect_true(file.exists(file.path(out_dir, "SI_Audit_Species.tsv")))
    expect_true(length(list.files(out_dir, pattern = "\\.pdf$")) > 0)

    unlink(temp_proj_dir, recursive = TRUE)
})

test_that("taxa_retention with export = FALSE writes nothing and returns plots", {
    temp_proj_dir <- tempfile(pattern = "kariocas_test_si_noexp_")
    dir.create(temp_proj_dir)
    mock_data_src <- system.file(
        "extdata/your_project_name/000_mpa_original",
        package = "karioCaS"
    )
    dir.create(file.path(temp_proj_dir, "000_mpa_original"), recursive = TRUE)
    file.copy(
        list.files(mock_data_src, full.names = TRUE),
        file.path(temp_proj_dir, "000_mpa_original")
    )
    suppressMessages(import_karioCaS(project_dir = temp_proj_dir))

    before <- list.files(temp_proj_dir, recursive = TRUE)
    res <- suppressMessages(
        taxa_retention(temp_proj_dir, detail_samples = "SAMPLE01", export = FALSE)
    )
    after <- list.files(temp_proj_dir, recursive = TRUE)

    expect_identical(before, after)
    expect_false(dir.exists(file.path(temp_proj_dir, "002_taxa_retention")))
    expect_s3_class(res, "kariocas_result")
    expect_null(res$output_dir)
    expect_length(res$paths, 0)
    expect_true("SAMPLE_Group_Retention" %in% names(res$plots))
    expect_true("SAMPLE01_All_Levels" %in% names(res$plots))
    expect_s3_class(res$plots[["SAMPLE01_All_Levels"]], "ggplot")
    expect_output(print(res), "kariocas_result")

    # Unknown sample names are a hard error, not a silent fallback
    expect_error(
        suppressMessages(
            taxa_retention(temp_proj_dir, detail_samples = "NOPE", export = FALSE)
        ),
        "Unknown sample name"
    )
    unlink(temp_proj_dir, recursive = TRUE)
})
