# Vamos montar um esqueleto de script R com três abordagens:
# 1. K-means refinado com PCA
# 2. DBSCAN (clustering por densidade)
# 3. Regras baseadas em faixas de OR e pvalue

# Criação do script em formato de string para exportar
script_r = """
# ------------------------------
# Script de Refinamento de Risco com ML - Alzheimer e SNP rs429358
# ------------------------------


if (!require("dplyr")) install.packages("dplyr")
if (!require("dbscan")) install.packages("dbscan")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("factoextra")) install.packages("factoextra")

library(dplyr)
library(dbscan)
library(ggplot2)
library(factoextra)

# Assumimos que o dataframe `dados_ml` já está criado com:
# - p_log = -log10(pvalue)
# - or = OR or BETA
# - name = DISEASE/TRAIT

# ------------------------------
# Pré-processamento
# ------------------------------

# Padronizar os dados
dados_ml_std <- scale(dados_ml[, c("p_log", "or")])

# ------------------------------
# 1. K-Means com PCA
# ------------------------------

pca <- prcomp(dados_ml_std, scale. = TRUE)
pca_df <- as.data.frame(pca$x[, 1:2])  # Usar os dois primeiros componentes

set.seed(42)
kmeans_pca <- kmeans(pca_df, centers = 3)

dados_ml$cluster_kmeans_pca <- as.factor(kmeans_pca$cluster)

# ------------------------------
# 2. DBSCAN
# ------------------------------

dbscan_res <- dbscan(pca_df, eps = 0.5, minPts = 5)
dados_ml$cluster_dbscan <- as.factor(dbscan_res$cluster)

# ------------------------------
# 3. Classificação por regras clínicas
# ------------------------------

dados_ml <- dados_ml %>%
  mutate(
    risk_manual = case_when(
      or < 0.8 & 10^(-p_log) < 5e-8 ~ "low_risk",
      or > 1.2 & 10^(-p_log) < 5e-8 ~ "high_risk",
      TRUE ~ "neutral"
    )
  )

# ------------------------------
# Visualizações (opcional)
# ------------------------------

# K-means com PCA
fviz_cluster(kmeans_pca, data = pca_df, geom = "point", main = "K-means + PCA")

# DBSCAN
ggplot(pca_df, aes(x = PC1, y = PC2, color = dados_ml$cluster_dbscan)) +
  geom_point(size = 2) +
  labs(title = "DBSCAN Clustering", color = "Cluster") +
  theme_minimal()

# ------------------------------
# Exportar ou continuar com integração ao Neo4j...
# ------------------------------
"""

# Exemplo: traduzir clusters para rótulos de risco
dados_ml <- dados_ml %>%
  mutate(
    risk_kmeans = case_when(
      cluster_kmeans_pca == 1 ~ "neutral",
      cluster_kmeans_pca == 2 ~ "high_risk",
      cluster_kmeans_pca == 3 ~ "low_risk"
    ),
    risk_dbscan = case_when(
      cluster_dbscan == 1 ~ "low_risk",
      cluster_dbscan == 2 ~ "high_risk",
      TRUE ~ "neutral"  # 0 ou NA vira neutro
    )
  )

dados_ml <- dados_ml %>%
  mutate(risco_consistente = ifelse(
    risk_manual == risk_kmeans & risk_manual == risk_dbscan,
    "CONSISTENTE", "VARIÁVEL"
  ))

consistentes <- dados_ml %>%
  filter(risco_consistente == "CONSISTENTE")

View(consistentes)

dados_ml <- dados_ml %>%
  mutate(risco_consistente = ifelse(`DISEASE/TRAIT` %in% consistentes$`DISEASE/TRAIT`, "SIM", "NAO"))

library(glue)
library(purrr)
library(neo4r)

# Conexão (ajuste se necessário)
con <- neo4j_api$new(
  url = "http://neo4j:7474",
  user = "neo4j",
  password = "admin123"  
)

# Lista com os valores
consistencia_list <- dados_ml %>%
  transmute(name = `DISEASE/TRAIT`, risco_consistente = risco_consistente) %>%
  purrr::transpose()

walk(consistencia_list, function(row) {
  query <- glue('
    MATCH (d:Disease_Trait {{name: "{row$name}"}})
    SET d.risco_consistente = "{row$risco_consistente}"
  ')
  call_neo4j(query, con)
})

# ----------------------------
# Consulta para alto risco
# ----------------------------
query_high <- "
MATCH (d:Disease_Trait)
WHERE d.risk = 'high_risk' AND d.risco_consistente = 'SIM'
RETURN d.name AS name, d.or AS or, d.pvalue AS pvalue, d.p_log AS p_log
"

resultado_alto <- call_neo4j(query_high, con, type = "row")

# Montar dataframe corretamente
high_risk_df <- data.frame(
  name = resultado_alto$name[[1]],
  or = resultado_alto$or[[1]],
  pvalue = resultado_alto$pvalue[[1]],
  p_log = resultado_alto$p_log[[1]],
  stringsAsFactors = FALSE
)

write.csv(high_risk_df, "high_risk_confirmed.csv", row.names = FALSE)

# ----------------------------
# Consulta para baixo risco
# ----------------------------
query_low <- "
MATCH (d:Disease_Trait)
WHERE d.risk = 'low_risk' AND d.risco_consistente = 'SIM'
RETURN d.name AS name, d.or AS or, d.pvalue AS pvalue, d.p_log AS p_log
"

resultado_baixo <- call_neo4j(query_low, con, type = "row")

# Montar dataframe corretamente
low_risk_df <- data.frame(
  name = resultado_baixo$name[[1]],
  or = resultado_baixo$or[[1]],
  pvalue = resultado_baixo$pvalue[[1]],
  p_log = resultado_baixo$p_log[[1]],
  stringsAsFactors = FALSE
)

write.csv(low_risk_df, "low_risk_confirmed.csv", row.names = FALSE)

resultado_bruto <- call_neo4j(query_high, con, type = "row")
str(resultado_bruto)

library(ggplot2)

# Adiciona coluna de grupo
high_risk_df$risk <- "High"
low_risk_df$risk <- "Low"

# Junta tudo
combined_df <- rbind(high_risk_df, low_risk_df)

library(patchwork)

p1 <- ggplot(low_risk_df, aes(x = p_log, y = or)) +
  geom_point(color = "#1ABC9C", size = 2) +
  labs(title = "Low Risk", x = "-log10(p-valor)", y = "OR") +
  theme_minimal()

p2 <- ggplot(high_risk_df, aes(x = p_log, y = or)) +
  geom_point(color = "#E74C3C", size = 2) +
  labs(title = "High Risk", x = "-log10(p-valor)", y = "OR") +
  theme_minimal()

p1 + p2  # lado a lado



