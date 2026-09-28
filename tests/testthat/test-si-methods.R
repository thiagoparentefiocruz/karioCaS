# tests/testthat/test-si-methods.R
# Every Stability Index method on a synthetic decay curve whose answers are
# computed by hand, plus the alternative methods end-to-end.
#
# Taxa per CS:  CS  0   20   40   60   80   90
#               n 1000  400  300  280  270  265
# Pct_Retained:     100   40   30   28   27  26.5
# Step loss:          0   60   10    2    1   0.5
# Tail tolerance (step losses at CS >= 50: 2, 1, 0.5):
#   mean 1.1667 + 1.5 * sd 0.7638 = 2.31
# Kneedle: largest gap below the chord (0,100)-(90,26.5) is at CS 20.
# Post-cliff / dynamic: first CS (after the steepest drop) with loss <= 2.31
#   is CS 60. Manual (toll 1.5): first CS > 0 with loss <= 1.5 is CS 80.
# (Tolerances avoid exact ties: 0.3 * 100 is not exactly 30 in floating point.)

.kcs_synthetic_curve <- function() {
    n <- c(`0` = 1000, `20` = 400, `40` = 300, `60` = 280, `80` = 270,
        `90` = 265)
    do.call(rbind, lapply(names(n), function(cs) {
        data.frame(
            Domain = "Bacteria", CS = as.numeric(cs),
            Taxon_Name = paste0("t", seq_len(n[[cs]])),
            stringsAsFactors = FALSE
        )
    }))
}

.kcs_si <- function(method, manual_toll = 1) {
    a <- karioCaS:::.si_domain_audit(
        .kcs_synthetic_curve(), "Bacteria", method, manual_toll, "S1",
        function(...) invisible(NULL)
    )
    list(
        primary = a$CS[a$SI_Type %in% "Primary_SI"],
        secondary = a$CS[a$SI_Type %in% "Secondary_SI_1"],
        audit = a
    )
}

test_that("audit columns and percentages are computed as documented", {
    a <- .kcs_si("kneedle")$audit
    expect_equal(a$Taxa_Count, c(1000, 400, 300, 280, 270, 265))
    expect_equal(a$Pct_Retained, c(100, 40, 30, 28, 27, 26.5))
    expect_equal(a$Step_Loss_Pct, c(0, 60, 10, 2, 1, 0.5))
})

test_that("kneedle finds the elbow and the post-cliff floor as secondary", {
    r <- .kcs_si("kneedle")
    expect_equal(r$primary, 20)
    expect_equal(r$secondary, 60)
})

test_that("postcliff, dynamic and manual follow their definitions", {
    expect_equal(.kcs_si("postcliff")$primary, 60)
    expect_equal(.kcs_si("dynamic")$primary, 60)
    expect_equal(.kcs_si("manual", manual_toll = 1.5)$primary, 80)
    # Per-domain list for manual
    expect_equal(
        .kcs_si("manual", manual_toll = list(Bacteria = 10.5))$primary, 40
    )
})

test_that("segmented returns a breakpoint inside the curve", {
    p <- .kcs_si("segmented")$primary
    expect_length(p, 1)
    expect_true(p %in% c(40, 60))
})

test_that("domains without enough CS points or without loss are skipped", {
    short <- .kcs_synthetic_curve()
    short <- short[short$CS %in% c(0, 20), ]
    expect_null(karioCaS:::.si_domain_audit(
        short, "Bacteria", "kneedle", 1, "S1", function(...) NULL
    ))
    flat <- .kcs_synthetic_curve()
    flat <- flat[flat$Taxon_Name %in% paste0("t", 1:10), ]
    expect_null(karioCaS:::.si_domain_audit(
        flat, "Bacteria", "kneedle", 1, "S1", function(...) NULL
    ))
})

test_that("every method runs end-to-end on the example data", {
    proj <- tempfile(pattern = "kariocas_test_methods_")
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
    suppressMessages(import_karioCaS(proj))
    for (m in c("postcliff", "segmented", "dynamic", "manual")) {
        res <- suppressMessages(taxa_retention(proj, method = m, export = FALSE))
        expect_true("Primary_SI" %in% res$data$SI_Type, info = m)
    }
    for (m in c("postcliff", "segmented")) {
        res <- suppressMessages(reads_per_taxa(proj, method = m, export = FALSE))
        expect_true("Primary_SI" %in% res$data$SI_Type, info = m)
    }
    # Detailed per-sample saturation panels, one per CS
    res <- suppressMessages(reads_per_taxa(
        proj,
        detail_samples = "SAMPLE01", export = FALSE
    ))
    expect_true(any(grepl("^SAMPLE01_CS[0-9]+_Saturation$", names(res$plots))))
    unlink(proj, recursive = TRUE)
})
