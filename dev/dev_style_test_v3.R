# ==============================================================================
# SCRIPT DE VALIDAÇÃO VISUAL V3 (Versão Final para Aprovação)
# ==============================================================================

# 1. Carregar Pacote e Dependências
devtools::load_all()
library(tidyverse)
library(patchwork)

# 2. Carregar Dados do Caminho Exato
input_path <- "/Users/thiagoparente/Desktop/2025_caetano_metagenomica_sample/INPE/000_karioCaS_input_matrix/karioCaS_input_matrix.tsv"

if(!file.exists(input_path)) stop("ERRO: Arquivo não encontrado no caminho: ", input_path)

message("Carregando dados de: ", input_path)

df_raw <- read_tsv(input_path, show_col_types = FALSE)
sample_cols <- names(df_raw)[grep("_CS", names(df_raw))]

df_long <- df_raw %>%
  pivot_longer(cols = all_of(sample_cols), names_to = "Sample_CS", values_to = "Counts") %>%
  separate(Sample_CS, into = c("Sample", "CS_Tag"), sep = "_CS", remove = FALSE) %>%
  mutate(CS = as.numeric(CS_Tag) * 10)

# ==============================================================================
# TESTE 1: Barras Horizontais (Resolution Style)
# ==============================================================================
df_bar <- df_long %>%
  filter(Domain == "Bacteria", Rank == "Genus", CS == 0) %>%
  group_by(Genus) %>% summarise(TotalReads = sum(Counts)) %>%
  arrange(desc(TotalReads)) %>% slice_head(n = 10)

p_bars <- ggplot(df_bar, aes(x = TotalReads, y = reorder(Genus, TotalReads))) +
  geom_col(fill = kariocas_colors$special["Level Reads"], width = 0.7) +
  scale_x_continuous(labels = label_k_number) +
  labs(
    title = "Teste 1: Resolução Taxonômica",
    subtitle = "Top 10 Genera | Bacteria | CS00", # SUBTÍTULO REINSERIDO
    x = "Reads",
    y = NULL
  ) +
  theme_kariocas() +
  theme(axis.text.y = element_text(face = "italic"))

# ==============================================================================
# TESTE 2: Linhas de Retenção (Retention Style)
# ==============================================================================
df_line <- data.frame(
  CS = seq(0, 100, 10),
  Retencao = c(100, 80, 50, 20, 10, 5, 1, 0.5, 0.1, 0.05, 0.01)
)

p_line <- ggplot(df_line, aes(x = CS, y = Retencao)) +
  geom_line(color = kariocas_colors$domains["Eukaryota"], linewidth = 1) +
  geom_point(size = 3, shape = kariocas_shapes["Species"]) +

  # AQUI ESTÁ A CORREÇÃO DOS DECIMAIS (Usa a nova função do Style)
  scale_y_kariocas_log10() +

  labs(
    title = "Teste 2: Curva de Retenção",
    subtitle = "Exemplo Eukaryota | Species", # SUBTÍTULO REINSERIDO
    x = kariocas_labels$y_confidence,
    y = kariocas_labels$y_log10_retained
  ) +
  theme_kariocas()

# ==============================================================================
# TESTE 3: Plot Vazio
# ==============================================================================
p_empty <- plot_kariocas_empty(
  title_text = "Viruses (Sem Dados)",
  subtitle_text = "Esqueleto mantido para alinhamento",
  x_label = "Confidence Score",
  y_label = "Reads"
)

# ==============================================================================
# TESTE 4: Heatmap (REINSERIDO)
# ==============================================================================
# Recorte pequeno para visualização
df_heat <- df_long %>%
  filter(Rank == "Class", Domain == "Archaea") %>%
  group_by(Class, CS) %>% summarise(Abund = sum(Counts), .groups = "drop")

p_heat <- ggplot(df_heat, aes(x = factor(CS), y = Class, fill = Abund)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colors = kariocas_colors$heatmap) +
  labs(
    title = "Teste 4: Heatmap",
    subtitle = "Archaea | Class | Abundância",
    x = "CS", y = NULL
  ) +
  theme_kariocas() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ==============================================================================
# MONTAGEM FINAL
# ==============================================================================
# Layout com 4 gráficos
layout_final <- (p_bars | p_line) / (p_empty | p_heat) +
  plot_annotation(
    title = "karioCaS Style Check V3",
    subtitle = "Validado: Inteiros no Log, Subtítulos, Heatmap presente e Caminho Correto",
    theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
  )

output_file <- "kariocas_style_check_v3.pdf"
ggsave(output_file, layout_final,
       width = kariocas_dims$width, height = kariocas_dims$height)

message("Validação V3 salva em: ", output_file)
system(paste("open", output_file))
