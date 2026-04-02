# tests/testthat/test-006_optimize.R

test_that("optimize_CS calculates SI and generates outputs successfully", {
  
  # 1. Setup: Diretório temporário
  temp_proj_dir <- tempfile(pattern = "kariocas_test_opt_")
  dir.create(temp_proj_dir)
  
  # 2. Localizar e copiar o Mock Data
  mock_data_src <- system.file("extdata/your_project_name/000_mpa_original", package = "karioCaS")
  if (mock_data_src == "") {
    mock_data_src <- file.path("..", "..", "inst", "extdata", "your_project_name", "000_mpa_original")
  }
  
  dir.create(file.path(temp_proj_dir, "000_mpa_original"), recursive = TRUE)
  file.copy(list.files(mock_data_src, full.names = TRUE), 
            file.path(temp_proj_dir, "000_mpa_original"))
  
  # 3. Execução: Precisamos rodar a importação silenciosamente primeiro
  suppressMessages(import_karioCaS(project_dir = temp_proj_dir))
  
  # 4. Rodar a otimização (O teste principal)
  expect_message(
    audit_df <- optimize_CS(project_dir = temp_proj_dir, tax_level = "Species"),
    "SUCCESS: CS Optimization completed."
  )
  
  # 5. Auditoria do Output Numérico
  # A. Retornou um data frame?
  expect_s3_class(audit_df, "data.frame")
  
  # B. A matemática encontrou o Primary SI?
  expect_true("Primary_SI" %in% audit_df$SI_Type)
  
  # C. Os arquivos físicos de auditoria foram salvos na pasta certa?
  out_dir <- file.path(temp_proj_dir, "006_optimize_CS")
  expect_true(file.exists(file.path(out_dir, "SI_Audit_Species.rds")))
  expect_true(file.exists(file.path(out_dir, "SI_Audit_Species.tsv")))
  
  # D. O PDF foi gerado? (Verifica se existe pelo menos 1 arquivo PDF)
  pdfs <- list.files(out_dir, pattern = "\\.pdf$")
  expect_true(length(pdfs) > 0)
  
  # 6. Cleanup
  unlink(temp_proj_dir, recursive = TRUE)
})
