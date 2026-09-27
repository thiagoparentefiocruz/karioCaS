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
        res <- upset_kariocas(project_dir = temp_proj_dir),
        "SUCCESS: 005_taxa_intersections_across_CS completed"
    )
    expect_s3_class(res, "kariocas_result")
    expect_true(length(res$plots) > 0)
    expect_true(all(file.exists(res$paths)))
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

test_that("upset_kariocas with export = FALSE returns UpSet objects and data", {
    temp_proj_dir <- tempfile(pattern = "kariocas_test_ups_noexp_")
    dir.create(file.path(temp_proj_dir, "000_mpa_original"), recursive = TRUE)
    file.copy(
        list.files(
            system.file("extdata/your_project_name/000_mpa_original",
                package = "karioCaS"
            ),
            full.names = TRUE
        ),
        file.path(temp_proj_dir, "000_mpa_original")
    )
    suppressMessages(import_karioCaS(project_dir = temp_proj_dir))

    res <- suppressMessages(
        upset_kariocas(temp_proj_dir, tax_level = "Genus", export = FALSE)
    )
    expect_false(dir.exists(
        file.path(temp_proj_dir, "005_taxa_intersections_across_CS")
    ))
    expect_s3_class(res, "kariocas_result")
    expect_length(res$paths, 0)
    expect_true("SAMPLE01_Bacteria_Genus" %in% names(res$plots))
    expect_s3_class(res$plots[["SAMPLE01_Bacteria_Genus"]], "upset")
    expect_true(all(c("sample", "Domain", "Taxon_Name", "N_CS") %in%
        colnames(res$data)))
    # N_CS counts the CS columns in which each taxon is present
    cs_cols <- grep("^CS[0-9]+$", colnames(res$data), value = TRUE)
    expect_true(length(cs_cols) >= 2)
    expect_equal(
        unname(res$data$N_CS),
        unname(rowSums(res$data[, cs_cols], na.rm = TRUE))
    )
    unlink(temp_proj_dir, recursive = TRUE)
})

test_that("a failed UpSet drawing is reported, not hidden behind SUCCESS", {
    log <- character(0)
    logger <- function(...) log <<- c(log, paste0(...))
    path <- file.path(tempdir(), "kariocas_bad_upset.pdf")
    out <- karioCaS:::.kcs_draw_upset_pdf(
        list(), path, "title", logger,
        draw = function(x) stop("simulated UpSetR failure")
    )
    expect_false(file.exists(path))
    expect_false(is.null(attr(out, "error")))
    expect_true(any(grepl("ERROR plotting", log)))
    expect_warning(
        karioCaS:::.kcs_warn_failed_plots("x: boom", function(...) NULL),
        "could not be drawn"
    )
})
