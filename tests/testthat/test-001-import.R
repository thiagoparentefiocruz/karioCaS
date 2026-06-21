# tests/testthat/test-000_import.R

test_that("import_karioCaS generates TreeSummarizedExperiment successfully", {
    # 1. Setup: Criar diret<U+00F3>rio tempor<U+00E1>rio para n<U+00E3>o sujar o reposit<U+00F3>rio (Regra Bioconductor)
    temp_proj_dir <- tempfile(pattern = "kariocas_test_")
    dir.create(temp_proj_dir)

    # 2. Localizar o Mock Data
    # Quando o pacote <U+00E9> testado, system.file acha a pasta extdata
    mock_data_src <- system.file("extdata/your_project_name/000_mpa_original", package = "karioCaS")

    # Como estamos desenvolvendo localmente e o pacote pode n<U+00E3>o estar instalado,
    # fazemos um fallback seguro para a pasta local se system.file falhar.
    if (mock_data_src == "") {
        mock_data_src <- file.path("..", "..", "inst", "extdata", "your_project_name", "000_mpa_original")
    }

    # Copiar o Mock Data para o diret<U+00F3>rio tempor<U+00E1>rio
    dir.create(file.path(temp_proj_dir, "000_mpa_original"), recursive = TRUE)
    file.copy(
        list.files(mock_data_src, full.names = TRUE),
        file.path(temp_proj_dir, "000_mpa_original")
    )

    # 3. Execu<U+00E7><U+00E3>o: Rodar a fun<U+00E7><U+00E3>o
    # Usamos expect_message para garantir que o c<U+00F3>digo chega at<U+00E9> o final
    expect_message(
        tse <- import_karioCaS(project_dir = temp_proj_dir),
        "SUCCESS: Import completed."
    )

    # 4. Auditoria (Os Testes de Fato)

    # A. O objeto <U+00E9> da classe certa para o Bioconductor?
    expect_s4_class(tse, "TreeSummarizedExperiment")

    # B. A fun<U+00E7><U+00E3>o leu os 6 arquivos corretamente? (ncol do TSE deve ser 6)
    expect_equal(ncol(tse), 6)

    # C. Os metadados das colunas (CS e Sample) foram extra<U+00ED>dos corretamente?
    expect_true(all(c("Sample_ID", "Confidence_Score") %in% colnames(SummarizedExperiment::colData(tse))))

    # D. Os arquivos f<U+00ED>sicos foram gerados na pasta tempor<U+00E1>ria?
    expect_true(file.exists(file.path(temp_proj_dir, "001_imported_matrix", "karioCaS_TSE.rds")))
    expect_true(file.exists(file.path(temp_proj_dir, "001_imported_matrix", "karioCaS_matrix_audit.tsv")))

    # 5. Cleanup: Apagar diret<U+00F3>rio tempor<U+00E1>rio
    unlink(temp_proj_dir, recursive = TRUE)
})
