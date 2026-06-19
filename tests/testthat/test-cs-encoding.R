# tests/testthat/test-cs-encoding.R
# Regression tests for the canonical Confidence Score (CS) unit system.
# Guards against the historical bug where CS = 1.0 (filename CS10) was
# misread as 0.1 and where CS = 1.0 arguments silently matched nothing.

test_that(".cs_filename_to_percent maps legacy, percent and decimal codes", {
    f <- karioCaS:::.cs_filename_to_percent
    # Legacy zero-padded tenths code (existing user data)
    expect_equal(f("00"), 0L)
    expect_equal(f("02"), 20L)
    expect_equal(f("09"), 90L)
    # The historical break point: CS10 must be 1.0 (100%), not 0.1 (10%)
    expect_equal(f("10"), 100L)
    # Explicit percentages
    expect_equal(f("50"), 50L)
    expect_equal(f("100"), 100L)
    # Explicit decimal fractions
    expect_equal(f("0.9"), 90L)
    expect_equal(f("0.05"), 5L)
    expect_equal(f("1.0"), 100L)
    # Out of range / unparseable
    expect_true(is.na(f("1.5")))
    expect_true(is.na(f("abc")))
})

test_that(".cs_arg_to_percent treats fractions and percentages correctly", {
    a <- karioCaS:::.cs_arg_to_percent
    expect_equal(a(1.0), 100L) # Kraken fraction for max stringency
    expect_equal(a(0.9), 90L)
    expect_equal(a(40), 40L) # percentage style (legacy argument usage)
    expect_equal(a(90), 90L)
    expect_equal(a(0), 0L)
    expect_true(is.na(a("auto")))
})

test_that("import_karioCaS reads a CS10 file as 1.0 (100%)", {
    temp_proj_dir <- tempfile(pattern = "kariocas_test_cs_")
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
    # Add a genuine maximum-stringency run named the natural way (CS10).
    file.copy(
        file.path(temp_proj_dir, "000_mpa_original", "SAMPLE01_CS09.mpa"),
        file.path(temp_proj_dir, "000_mpa_original", "SAMPLE01_CS10.mpa")
    )

    tse <- suppressMessages(import_karioCaS(project_dir = temp_proj_dir))
    cs_vals <- SummarizedExperiment::colData(tse)$Confidence_Score
    expect_true(100 %in% cs_vals)
    expect_false(10 %in% cs_vals) # must NOT be the old 0.1 misread

    unlink(temp_proj_dir, recursive = TRUE)
})
