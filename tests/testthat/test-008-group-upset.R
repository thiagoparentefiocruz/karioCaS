# tests/testthat/test-006-group-upset.R
# group_upset(): cross-sample UpSet within a biological group. Default source is
# the final mosaic; CS = single Confidence Score.

# Build a project with TWO samples in the same group (SAMPLE01 + SAMPLE02) so a
# cross-sample comparison is possible (mock data ships a single sample).
.kcs_setup_two_sample_proj <- function() {
    proj <- tempfile(pattern = "kariocas_test_gup_")
    in_dir <- file.path(proj, "000_mpa_original")
    dir.create(in_dir, recursive = TRUE)
    src <- system.file(
        "extdata/your_project_name/000_mpa_original",
        package = "karioCaS"
    )
    if (src == "") {
        src <- file.path(
            "..", "..", "inst", "extdata",
            "your_project_name", "000_mpa_original"
        )
    }
    files <- list.files(src, full.names = TRUE)
    file.copy(files, in_dir)
    # Duplicate SAMPLE01_* as SAMPLE02_* (same group "SAMPLE").
    for (f in list.files(in_dir, full.names = TRUE)) {
        file.copy(f, sub("SAMPLE01", "SAMPLE02", f))
    }
    suppressMessages(import_karioCaS(project_dir = proj))
    proj
}

test_that("group_upset(CS = ...) builds membership with Core/Unique categories", {
    proj <- .kcs_setup_two_sample_proj()

    expect_message(
        res <- group_upset(project_dir = proj, CS = 40),
        "SUCCESS: 008_taxa_intersections_across_samples completed"
    )
    expect_s3_class(res, "kariocas_result")
    memb <- res$data
    expect_s3_class(memb, "data.frame")
    expect_true(all(
        c("Group", "Domain", "Taxon", "N_Samples", "Category") %in%
            colnames(memb)
    ))
    # Two identical samples -> every taxon is Core (present in both).
    expect_true("Core" %in% memb$Category)
    expect_equal(max(memb$N_Samples), 2)

    out_dir <- file.path(proj, "008_taxa_intersections_across_samples", "SAMPLE")
    expect_true(any(grepl("_CS40_SampleUpSet\\.pdf$", list.files(out_dir))))
    expect_true(any(grepl("_CS40_membership\\.tsv$", list.files(out_dir))))

    unlink(proj, recursive = TRUE)
})

test_that("group_upset default reads the final mosaic", {
    proj <- .kcs_setup_two_sample_proj()
    suppressMessages(taxa_retention(project_dir = proj, tax_level = "Species"))
    suppressMessages(retrieve_selected_taxa(
        project_dir = proj,
        CS_A = "auto", CS_B = "auto", CS_E = "auto", CS_V = "auto"
    ))

    expect_message(
        group_upset(project_dir = proj),
        "SUCCESS: 008_taxa_intersections_across_samples completed"
    )
    out_dir <- file.path(proj, "008_taxa_intersections_across_samples", "SAMPLE")
    expect_true(any(grepl("Final_Mosaic_SampleUpSet\\.pdf$", list.files(out_dir))))

    unlink(proj, recursive = TRUE)
})

test_that("group_upset rejects an invalid rank", {
    proj <- .kcs_setup_two_sample_proj()
    expect_error(
        suppressMessages(group_upset(project_dir = proj, tax_level = "Nope", CS = 40)),
        "Invalid 'tax_level'"
    )
    unlink(proj, recursive = TRUE)
})

test_that("group_upset warns instead of reporting success with one sample", {
    proj <- tempfile(pattern = "kariocas_test_gup_single_")
    dir.create(file.path(proj, "000_mpa_original"), recursive = TRUE)
    file.copy(
        list.files(
            system.file("extdata/your_project_name/000_mpa_original",
                package = "karioCaS"
            ),
            full.names = TRUE
        ),
        file.path(proj, "000_mpa_original")
    )
    suppressMessages(import_karioCaS(project_dir = proj))
    expect_warning(
        res <- suppressMessages(group_upset(proj, CS = 40, export = FALSE)),
        "no group comparisons"
    )
    expect_length(res$plots, 0)
    unlink(proj, recursive = TRUE)
})

test_that("group_upset with export = FALSE writes nothing", {
    proj <- .kcs_setup_two_sample_proj()
    res <- suppressMessages(group_upset(proj, CS = 40, export = FALSE))
    expect_false(dir.exists(
        file.path(proj, "008_taxa_intersections_across_samples")
    ))
    expect_true(length(res$plots) > 0)
    expect_s3_class(res$plots[[1]], "upset")
    unlink(proj, recursive = TRUE)
})
