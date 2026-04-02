# tests/testthat/test-000_import.R

test_that("import_karioCaS generates TreeSummarizedExperiment successfully", {
  
  # 1. Setup: Criar diretório temporário para não sujar o repositório (Regra Bioconductor)
  temp_proj_dir <- tempfile(pattern = "kariocas_test_")
  dir.create(temp_proj_dir)
  
  # 2. Localizar o Mock Data
  # Quando o pacote é testado, system.file acha a pasta extdata
  mock_data_src <- system.file("extdata/your_project_name/000_mpa_original", package = "karioCaS")
  
  # Como estamos desenvolvendo localmente e o pacote pode não estar instalado,
  # fazemos um fallback seguro para a pasta local se system.file falhar.
  if (mock_data_src == "") {
    mock_data_src <- file.path("..", "..", "inst", "extdata", "your_project_name", "000_mpa_original")
  }
  
  # Copiar o Mock Data para o diretório temporário
  dir.create(file.path(temp_proj_dir, "000_mpa_original"), recursive = TRUE)
  file.copy(list.files(mock_data_src, full.names = TRUE), 
            file.path(temp_proj_dir, "000_mpa_original"))
  
  # 3. Execução: Rodar a função
  # Usamos expect_message para garantir que o código chega até o final
  expect_message(
    tse <- import_karioCaS(project_dir = temp_proj_dir),
    "SUCCESS: Import completed."
  )
  
  # 4. Auditoria (Os Testes de Fato)
  
  # A. O objeto é da classe certa para o Bioconductor?
  expect_s4_class(tse, "TreeSummarizedExperiment")
  
  # B. A função leu os 6 arquivos corretamente? (ncol do TSE deve ser 6)
  expect_equal(ncol(tse), 6)
  
  # C. Os metadados das colunas (CS e Sample) foram extraídos corretamente?
  expect_true(all(c("Sample_ID", "Confidence_Score") %in% colnames(SummarizedExperiment::colData(tse))))
  
  # D. Os arquivos físicos foram gerados na pasta temporária?
  expect_true(file.exists(file.path(temp_proj_dir, "000_karioCaS_input_matrix", "karioCaS_TSE.rds")))
  expect_true(file.exists(file.path(temp_proj_dir, "000_karioCaS_input_matrix", "karioCaS_matrix_audit.tsv")))
  
  # 5. Cleanup: Apagar diretório temporário
  unlink(temp_proj_dir, recursive = TRUE)
})
