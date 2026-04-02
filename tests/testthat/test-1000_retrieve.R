# tests/testthat/test-1000_retrieve.R

test_that("retrieve_selected_taxa applies UX logic correctly", {
  
  # 1. Setup: Diretório temporário
  temp_proj_dir <- tempfile(pattern = "kariocas_test_ret_")
  dir.create(temp_proj_dir)
  
  # 2. Localizar e copiar o Mock Data
  mock_data_src <- system.file("extdata/your_project_name/000_mpa_original", package = "karioCaS")
  if (mock_data_src == "") {
    mock_data_src <- file.path("..", "..", "inst", "extdata", "your_project_name", "000_mpa_original")
  }
  
  dir.create(file.path(temp_proj_dir, "000_mpa_original"), recursive = TRUE)
  file.copy(list.files(mock_data_src, full.names = TRUE), 
            file.path(temp_proj_dir, "000_mpa_original"))
  
  # 3. Execução dos pré-requisitos (Importação e Otimização) em modo silencioso
  suppressMessages(import_karioCaS(project_dir = temp_proj_dir))
  suppressMessages(optimize_CS(project_dir = temp_proj_dir, tax_level = "Species"))
  
  # 4. Rodar o Mosaico Final (Testando auto, secondary e manual numérico simultaneamente)
  expect_message(
    retrieve_selected_taxa(
      project_dir = temp_proj_dir, 
      tax_level = "Species",
      CS_B = "auto",       # Deve puxar o Primary SI
      CS_A = "secondary",  # Deve puxar o Secondary SI (se existir)
      CS_E = 40,           # Override manual
      CS_V = 0             # Fallback numérico
    ),
    "SUCCESS: Process completed."
  )
  
  # 5. Auditoria do Output Final
  out_dir <- file.path(temp_proj_dir, "1000_final_selection")
  
  # A. Os arquivos finais foram gerados?
  expect_true(file.exists(file.path(out_dir, "SAMPLE01_karioCaS_Mosaic.mpa")))
  expect_true(file.exists(file.path(out_dir, "SAMPLE01_karioCaS_Mosaic.tsv")))
  
  # 6. Cleanup
  unlink(temp_proj_dir, recursive = TRUE)
})
