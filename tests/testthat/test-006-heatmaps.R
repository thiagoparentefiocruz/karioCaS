# tests/testthat/test-006-heatmaps.R
# heatmaps_karioCaS(): survivors at a target CS, loss groups, and messages.

.kcs_setup_heatmap_proj <- function() {
    proj <- tempfile(pattern = "kariocas_test_hm_")
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
    proj
}

test_that("heatmaps_karioCaS exports one PDF per sample by default", {
    proj <- .kcs_setup_heatmap_proj()
    expect_message(
        res <- heatmaps_karioCaS(proj, confidence_score = 40, top_n = 10),
        "SUCCESS: 006_relative_abundance_across_CS completed"
    )
    expect_s3_class(res, "kariocas_result")
    out_dir <- file.path(proj, "006_relative_abundance_across_CS")
    pdfs <- list.files(out_dir, pattern = "\\.pdf$")
    expect_equal(pdfs, "SAMPLE01_Heatmap_Genus_CS40.pdf")
    expect_true(all(file.exists(res$paths)))
    unlink(proj, recursive = TRUE)
})

test_that("heatmaps_karioCaS data: top survivors, loss groups, relative abundance", {
    proj <- .kcs_setup_heatmap_proj()
    res <- suppressMessages(heatmaps_karioCaS(
        proj,
        analysis_rank = "Species", confidence_score = 40, top_n = 5,
        export = FALSE
    ))
    expect_false(dir.exists(file.path(proj, "006_relative_abundance_across_CS")))
    d <- res$data
    expect_true(all(c(
        "sample", "Domain", "Taxon_Name", "CS", "Counts", "Rel_Abund",
        "Target_CS"
    ) %in% colnames(d)))
    expect_true(all(d$Target_CS == 40))
    # Only CS levels up to the target are shown
    expect_true(all(d$CS <= 40))
    expect_true(all(d$Rel_Abund >= 0 & d$Rel_Abund <= 100))
    # At most top_n individual survivors per domain; the rest are loss groups
    individual <- d[!grepl("Recovered only in|Lowest abundance", d$Taxon_Name), ]
    per_dom <- tapply(individual$Taxon_Name, individual$Domain, function(x) {
        length(unique(x))
    })
    expect_true(all(per_dom <= 5))
    expect_true(any(grepl("Recovered only in", d$Taxon_Name)))
    expect_s3_class(res$plots[[1]], "ggplot")
    unlink(proj, recursive = TRUE)
})

test_that("heatmaps_karioCaS is explicit about the target CS it uses", {
    proj <- .kcs_setup_heatmap_proj()
    # NULL -> highest available CS, stated in a message
    expect_message(
        heatmaps_karioCaS(proj, export = FALSE),
        "highest available CS"
    )
    # Above the highest CS -> warning, then the maximum is used
    expect_warning(
        res <- suppressMessages(
            heatmaps_karioCaS(proj, confidence_score = 99, export = FALSE)
        ),
        "exceeds the highest available CS"
    )
    expect_true(all(res$data$Target_CS == max(res$data$CS)))
    # Invalid inputs
    expect_error(heatmaps_karioCaS(proj, confidence_score = "abc"), "Invalid")
    expect_error(heatmaps_karioCaS(proj, analysis_rank = NULL), "analysis_rank")
    unlink(proj, recursive = TRUE)
})

# MPA counts are cumulative (a genus row already includes its species' reads).
# Genus-level values must equal the genus's own MPA row, not genus + species.
.kcs_genus_rows <- function(cs_code) {
    f <- system.file(
        "extdata", "your_project_name", "000_mpa_original",
        paste0("SAMPLE01_CS", cs_code, ".mpa"),
        package = "karioCaS"
    )
    mpa <- utils::read.delim(
        f,
        header = FALSE, quote = "", comment.char = "",
        col.names = c("Taxonomy", "Counts"), stringsAsFactors = FALSE
    )
    g <- mpa[grepl("\\|g__[^|]+$", mpa$Taxonomy), ]
    data.frame(
        Domain = sub("^d__([^|]+).*", "\\1", g$Taxonomy),
        Genus = sub(".*\\|g__", "", g$Taxonomy),
        Counts = g$Counts,
        stringsAsFactors = FALSE
    )
}

test_that("genus heatmap counts equal the genus MPA rows (no double counting)", {
    proj <- .kcs_setup_heatmap_proj()
    res <- suppressMessages(heatmaps_karioCaS(
        proj,
        analysis_rank = "Genus", confidence_score = 40, top_n = 10000,
        export = FALSE
    ))
    got <- res$data[res$data$CS == 40 & res$data$Counts > 0 &
        !grepl("Recovered only in|Lowest abundance", res$data$Taxon_Name), ]
    # Sum homonymous genera (same name under different lineages), as the
    # heatmap groups by name within a domain.
    expected <- stats::aggregate(
        Counts ~ Domain + Genus,
        data = .kcs_genus_rows("04"), FUN = sum
    )
    m <- merge(
        data.frame(
            Domain = got$Domain, Genus = as.character(got$Taxon_Name),
            Got = got$Counts
        ),
        expected,
        by = c("Domain", "Genus")
    )
    expect_true(nrow(m) > 0)
    expect_equal(nrow(m), nrow(got))
    expect_equal(m$Got, m$Counts)
    unlink(proj, recursive = TRUE)
})
