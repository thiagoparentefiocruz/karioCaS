# tests/testthat/test-003-reads.R
# reads_per_taxa(): saturation-only group overlays + optimal-minimum-reads audit.

test_that(".rpt_cutoffs_for always includes low anchors and adapts the top", {
    f <- karioCaS:::.rpt_cutoffs_for
    anchors <- c(1, 2, 3, 4, 5, 7, 10)
    # Low anchors present even when max reads < 10
    expect_equal(f(3), anchors)
    expect_equal(f(8), anchors)
    # Above 10: log-spaced up to just past max, anchors still present
    out <- f(250)
    expect_true(all(anchors %in% out))
    expect_true(max(out) >= 250) # extends past the data max
    expect_true(all(out > 0))
    expect_false(is.unsorted(out))
})

test_that(".si_reads_elbow finds a low-count elbow on a singleton-heavy curve", {
    # Many taxa vanish by a few reads (rare tail), then a stable core persists.
    cutoffs <- c(1, 3, 5, 10, 30, 100, 1000)
    taxa <- c(1000, 300, 120, 60, 40, 35, 30) # steep early drop, then flat
    el <- karioCaS:::.si_reads_elbow(
        data.frame(Cutoff = cutoffs, Taxa_Count = taxa),
        method = "kneedle"
    )
    expect_false(is.null(el))
    # Elbow should sit in the low-count region, not out at 1000.
    expect_true(el$opt >= 3 && el$opt <= 100)
})

test_that("reads_per_taxa writes the optimal-reads audit and saturation plots", {
    temp_proj_dir <- tempfile(pattern = "kariocas_test_reads_")
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
        audit <- reads_per_taxa(
            project_dir = temp_proj_dir, analysis_level = "Species"
        ),
        "SUCCESS: Cutoff analysis completed."
    )
    expect_s3_class(audit, "data.frame")
    expect_true(all(
        c("Sample", "CS", "Domain", "Cutoff", "SI_Type") %in% colnames(audit)
    ))
    expect_true("Primary_SI" %in% audit$SI_Type)

    out_dir <- file.path(temp_proj_dir, "003_cutoffs")
    expect_true(file.exists(file.path(out_dir, "Reads_Audit_Species.tsv")))
    expect_true(file.exists(file.path(out_dir, "Reads_Audit_Species.rds")))
    # Saturation only -- no Rare_Taxa outputs anymore.
    expect_false(any(grepl("Rare", list.files(out_dir))))
    expect_true(any(grepl("Saturation\\.pdf$", list.files(out_dir))))

    unlink(temp_proj_dir, recursive = TRUE)
})
