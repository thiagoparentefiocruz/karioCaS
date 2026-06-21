# tests/testthat/test-1000_retrieve.R

test_that("retrieve_selected_taxa applies UX logic correctly", {
    # 1. Setup: Diret<U+00F3>rio tempor<U+00E1>rio
    temp_proj_dir <- tempfile(pattern = "kariocas_test_ret_")
    dir.create(temp_proj_dir)

    # 2. Localizar e copiar o Mock Data
    mock_data_src <- system.file("extdata/your_project_name/000_mpa_original", package = "karioCaS")
    if (mock_data_src == "") {
        mock_data_src <- file.path("..", "..", "inst", "extdata", "your_project_name", "000_mpa_original")
    }

    dir.create(file.path(temp_proj_dir, "000_mpa_original"), recursive = TRUE)
    file.copy(
        list.files(mock_data_src, full.names = TRUE),
        file.path(temp_proj_dir, "000_mpa_original")
    )

    # 3. Execu<U+00E7><U+00E3>o dos pr<U+00E9>-requisitos (Importa<U+00E7><U+00E3>o e SI) em modo silencioso.
    # taxa_retention() agora computa o SI audit (Step 001).
    suppressMessages(import_karioCaS(project_dir = temp_proj_dir))
    suppressMessages(taxa_retention(project_dir = temp_proj_dir, tax_level = "Species"))

    # 4. Rodar o Mosaico Final (Testando auto, secondary e manual num<U+00E9>rico simultaneamente)
    expect_message(
        retrieve_selected_taxa(
            project_dir = temp_proj_dir,
            tax_level = "Species",
            CS_B = "auto", # Deve puxar o Primary SI
            CS_A = "secondary", # Deve puxar o Secondary SI (se existir)
            CS_E = 40, # Override manual
            CS_V = 0 # Fallback num<U+00E9>rico
        ),
        "SUCCESS: Process completed."
    )

    # 5. Auditoria do Output Final
    out_dir <- file.path(temp_proj_dir, "004_final_mosaic")

    # A. Os arquivos finais foram gerados? (.mpa em mpa/, .tsv em tsv/)
    expect_true(file.exists(
        file.path(out_dir, "mpa", "SAMPLE01_karioCaS_Mosaic.mpa")
    ))
    expect_true(file.exists(
        file.path(out_dir, "tsv", "SAMPLE01_karioCaS_Mosaic.tsv")
    ))

    # B. The mosaic keeps ALL ranks even though tax_level = "Species" was given
    #    (tax_level no longer filters the output).
    m <- read.delim(
        file.path(out_dir, "tsv", "SAMPLE01_karioCaS_Mosaic.tsv"),
        check.names = FALSE
    )
    tax <- as.character(m[[1]])
    expect_true(any(!grepl("[|]s__", tax))) # parent-rank rows present

    # 6. Cleanup
    unlink(temp_proj_dir, recursive = TRUE)
})

test_that("retrieve_selected_taxa pulls optimal min-reads from Reads_Audit", {
    temp_proj_dir <- tempfile(pattern = "kariocas_test_autoreads_")
    dir.create(file.path(temp_proj_dir, "000_mpa_original"), recursive = TRUE)
    mock_data_src <- system.file(
        "extdata/your_project_name/000_mpa_original",
        package = "karioCaS"
    )
    if (mock_data_src == "") {
        mock_data_src <- file.path(
            "..", "..", "inst", "extdata", "your_project_name", "000_mpa_original"
        )
    }
    file.copy(
        list.files(mock_data_src, full.names = TRUE),
        file.path(temp_proj_dir, "000_mpa_original")
    )

    # Prerequisites: CS audit (001) and Reads audit (003)
    suppressMessages(import_karioCaS(project_dir = temp_proj_dir))
    suppressMessages(taxa_retention(project_dir = temp_proj_dir, tax_level = "Species"))
    suppressMessages(reads_per_taxa(project_dir = temp_proj_dir, analysis_level = "Species"))

    # Fully data-driven: auto CS + auto reads for every domain
    expect_message(
        retrieve_selected_taxa(
            project_dir = temp_proj_dir, tax_level = "Species",
            CS_A = "auto", reads_min_A = "auto",
            CS_B = "auto", reads_min_B = "auto",
            CS_E = "auto", reads_min_E = "auto",
            CS_V = "auto", reads_min_V = "auto"
        ),
        "SUCCESS: Process completed."
    )
    out_dir <- file.path(temp_proj_dir, "004_final_mosaic")
    expect_true(file.exists(
        file.path(out_dir, "mpa", "SAMPLE01_karioCaS_Mosaic.mpa")
    ))

    unlink(temp_proj_dir, recursive = TRUE)
})

test_that(".rst_resolve_reads resolves auto/secondary/manual correctly", {
    ra <- data.frame(
        Sample = "S1", CS = 40, Domain = "Bacteria",
        Cutoff = c(1, 10, 100),
        SI_Type = c(NA, "Primary_SI", "Secondary_SI_1"),
        stringsAsFactors = FALSE
    )
    noop <- function(...) invisible(NULL)
    expect_equal(
        karioCaS:::.rst_resolve_reads("auto", "Bacteria", "S1", 40, ra, noop)$val, 10
    )
    expect_equal(
        karioCaS:::.rst_resolve_reads("secondary", "Bacteria", "S1", 40, ra, noop)$val, 100
    )
    expect_equal(
        karioCaS:::.rst_resolve_reads(5, "Bacteria", "S1", 40, ra, noop)$val, 5
    )
    # Missing audit row -> fallback to 0
    expect_equal(
        karioCaS:::.rst_resolve_reads("auto", "Viruses", "S1", 40, ra, noop)$val, 0
    )
})
