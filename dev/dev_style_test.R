# ==============================================================================
# SCRIPT DE VALIDAÇÃO VISUAL DO KARIOCAS
# ==============================================================================

# 1. Carregar Pacote Local e Dependências
devtools::load_all() # Carrega o kariocas_style.R
library(tidyverse)
library(patchwork)

# 2. Carregar Dados Reais (Ajuste o caminho se necessário)
# O script assume que o arquivo está na pasta '000_karioCaS_input_matrix'
input_file <- "/Users/thiagoparente/Desktop/2025_caetano_metagenomica_sample/INPE/000_karioCaS_input_matrix/karioCaS_input_matrix.tsv"

if(!file.exists(input_file)) {
  # Se não achar, tenta criar dados dummy só para o teste visual
  warning("Arquivo real não encontrado. Usando dados fictícios para teste visual.")
  df <- expand.grid(
    Domain = names(kariocas_colors$domains),
    Rank = names(kariocas_colors$ranks),
    sample = c("S1", "S2")
  ) %>% mutate(Counts = runif(n(), 10, 1000))
} else {
  df <- read_tsv(input_file, show_col_types = FALSE) %>%
    # Transforma de Wide para Long se necessário (o input matrix tsv é wide)
    pivot_longer(cols = starts_with("INPE_"), names_to = "Sample", values_to = "Counts")
}

# ==============================================================================
# TESTE A: Cores dos Ranks e Tipografia (Theme)
# ==============================================================================
p1 <- df %>%
  filter(Rank %in% names(kariocas_colors$ranks)) %>%
  group_by(Rank) %>% summarise(Total = sum(Counts)) %>%
  mutate(Rank = factor(Rank, levels = names(kariocas_colors$ranks))) %>%
  ggplot(aes(x = Rank, y = Total, fill = Rank)) +
  geom_col() +
  scale_fill_manual(values = kariocas_colors$ranks) +
  labs(title = "Teste A: Paleta de Ranks (Okabe-Ito)",
       subtitle = "Verificar cores, tamanho do título e eixos",
       y = "Total Counts (Log10)") +
  scale_y_log10() +
  theme_kariocas()

# ==============================================================================
# TESTE B: Cores dos Domínios e Facetas
# ==============================================================================
p2 <- df %>%
  filter(Domain %in% names(kariocas_colors$domains)) %>%
  group_by(Domain, Sample) %>% summarise(Total = sum(Counts)) %>%
  ggplot(aes(x = Domain, y = Total, fill = Domain)) +
  geom_boxplot(alpha = 0.8) +
  scale_fill_manual(values = kariocas_colors$domains) +
  labs(title = "Teste B: Paleta de Domínios (NPG)",
       subtitle = "Verificar cores e estilo das Facetas (Strip cinza)") +
  facet_wrap(~Domain, scales = "free") +
  theme_kariocas() +
  theme(axis.text.x = element_blank()) # Remove texto redundante

# ==============================================================================
# TESTE C: Shapes e Linhas (Scatter Plot)
# ==============================================================================
# Criando dados dummy para simular retenção (curva)
dummy_trend <- expand.grid(x = 1:10, Rank = names(kariocas_colors$ranks)) %>%
  mutate(y = (11-x) * runif(n(), 0.9, 1.1))

p3 <- ggplot(dummy_trend, aes(x = x, y = y, color = Rank, shape = Rank)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 3) +
  scale_color_manual(values = kariocas_colors$ranks) +
  scale_shape_manual(values = kariocas_shapes) +
  labs(title = "Teste C: Shapes e Linhas",
       subtitle = "Verificar consistência dos símbolos geométricos") +
  theme_kariocas()

# ==============================================================================
# TESTE D: Gradiente de Heatmap
# ==============================================================================
p4 <- df %>%
  head(100) %>% # Pegar apenas um subconjunto
  ggplot(aes(x = Sample, y = Taxonomy, fill = Counts)) +
  geom_tile() +
  scale_fill_gradientn(colors = kariocas_colors$heatmap) +
  labs(title = "Teste D: Gradiente Heatmap",
       subtitle = "Cyan (Baixo) -> Branco -> Vermelho (Alto)") +
  theme_kariocas() +
  theme(axis.text.y = element_blank(), axis.text.x = element_blank())

# ==============================================================================
# MONTAGEM FINAL E EXPORTAÇÃO
# ==============================================================================
layout_final <- (p1 | p2) / (p3 | p4) +
  plot_annotation(
    title = "karioCaS Visual Style Sheet Check",
    theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
  )

# Salvar usando as dimensões do pacote
output_file <- "kariocas_style_check.pdf"
ggsave(output_file, layout_final,
       width = kariocas_dims$width,
       height = kariocas_dims$height)

message("Style Check salvo em: ", output_file)
system(paste("open", output_file)) # Tenta abrir automaticamente no Mac
