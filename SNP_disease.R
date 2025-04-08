# Instalar pacotes se necessário
if (!require("dplyr")) install.packages("dplyr")
if (!require("readr")) install.packages("readr")
if (!require("neo4r")) install.packages("neo4r")
if (!require("ggraph")) install.packages("ggraph")
if (!require("tidygraph")) install.packages("tidygraph")
if (!require("glue")) install.packages("glue")
if (!require("purrr")) install.packages("purrr")

# Carregar pacotes
library(dplyr)
library(readr)
library(neo4r)
library(ggraph)
library(tidygraph)
library(glue)
library(purrr)

# Baixar o catálogo GWAS
url <- "https://www.ebi.ac.uk/gwas/api/search/downloads/full"
dest <- "~/dados/gwas_catalog.tsv"
download.file(url, destfile = dest)

# Carregar os dados
gwas <- read_tsv("/home/rstudio/dados/gwas_catalog.tsv", show_col_types = FALSE)

# Encontrar SNP mais frequente relacionado a Alzheimer
alzheimer <- gwas %>% filter(grepl("Alzheimer", `DISEASE/TRAIT`, ignore.case = TRUE))
top_snp <- alzheimer %>% count(SNPS, sort = TRUE) %>% slice(1) %>% pull(SNPS)

# Filtrar todos os dados do SNP mais frequente (ex: rs429358)
snp_data <- gwas %>%
  filter(SNPS == top_snp) %>%
  select(SNPS, `DISEASE/TRAIT`, `P-VALUE`, `OR or BETA`, `MAPPED_GENE`, CHR_ID, CHR_POS)%>%
  filter(!is.na(`MAPPED_GENE`)) %>%
  mutate(
    SNPS = gsub("\"", "", SNPS),
    `MAPPED_GENE` = toupper(gsub("\"", "", `MAPPED_GENE`)),
    `P-VALUE` = as.numeric(`P-VALUE`),
    `OR or BETA` = as.numeric(`OR or BETA`)
  )

# Conectar ao Neo4j
con <- neo4j_api$new(
  url = "http://neo4j:7474",
  user = "neo4j",
  password = "admin123"
)

# Criar lista para envio
snp_batch <- snp_data %>%
  pmap(function(SNPS, `DISEASE/TRAIT`, `P-VALUE`, `OR or BETA`, MAPPED_GENE, CHR_ID, CHR_POS) {
    list(
      snp = SNPS,
      gene = MAPPED_GENE,
      pval = `P-VALUE`,
      or = `OR or BETA`,
      disease = `DISEASE/TRAIT`,
      chr = CHR_ID,
      pos = CHR_POS
    )
  })


# Enviar para o Neo4j

walk(snp_batch, function(row) {
  query <- glue('
    MERGE (s:SNP {{id: "{row$snp}"}})
    SET s.pvalue = {row$pval},
        s.effect_size = {row$or},
        s.chr = "{row$chr}",
        s.position = "{row$pos}"
        
    MERGE (d:Disease_Trait {{name: "{row$disease}"}})
    MERGE (g:Gene {{name: "{row$gene}"}})
    
    MERGE (s)-[:ASSOCIATED_WITH]->(d)
    MERGE (d)-[:TARGETS]->(g)
  ')
  
  call_neo4j(query, con)
})

results <- call_neo4j(
  "MATCH (s:SNP {id: 'rs429358'})-[:ASSOCIATED_WITH]->(d:Disease_Trait)
   RETURN s, d",
  con,
  type = "graph"
)

library(dplyr)

dados_ml <- snp_completo %>%
  filter(SNPS == "rs429358", !is.na(`P-VALUE`), !is.na(`OR or BETA`)) %>%
  mutate(
    p_log = -log10(`P-VALUE`),
    or = `OR or BETA`
  ) %>%
  select(`DISEASE/TRAIT`, p_log, or)

# Padronizar
dados_ml_std <- scale(dados_ml[, c("p_log", "or")])
rownames(dados_ml_std) <- dados_ml$`DISEASE/TRAIT`

set.seed(42)
kmeans_result <- kmeans(dados_ml_std, centers = 3)

# Adicionar cluster aos dados
dados_ml$cluster <- as.factor(kmeans_result$cluster)

cluster_list <- dados_ml %>%
  transmute(name = `DISEASE/TRAIT`, cluster = as.character(cluster)) %>%
  purrr::transpose()

library(glue)

walk(cluster_list, function(row) {
  query <- glue('
    MATCH (d:Disease_Trait {{name: "{row$name}"}})
    SET d.cluster = "{row$cluster}"
  ')
  
  call_neo4j(query, con)
})

dados_ml$risk <- case_when(
  dados_ml$cluster == 2 ~ "high_risk",   # 🟢 verde – risco alto
  dados_ml$cluster == 3 ~ "low_risk",    # 🔵 azul – efeito protetor
  TRUE ~ "neutral"                       # 🔴 vermelho – efeito leve/neutro
)

risk_list <- dados_ml %>%
  transmute(name = `DISEASE/TRAIT`, risk = risk) %>%
  purrr::transpose()

walk(risk_list, function(row) {
  query <- glue('
    MATCH (d:Disease_Trait {{name: "{row$name}"}})
    SET d.risk = "{row$risk}"
  ')
  call_neo4j(query, con)
})

dados_ml <- dados_ml %>%
  mutate(pvalue = 10^(-p_log))


# Criar lista para envio
or_pval_list <- dados_ml %>%
  transmute(name = `DISEASE/TRAIT`, or = or, pvalue = pvalue) %>%
  purrr::transpose()

# Enviar ao Neo4j
walk(or_pval_list, function(row) {
  query <- glue('
    MATCH (d:Disease_Trait {{name: "{row$name}"}})
    SET d.or = {row$or},
        d.pvalue = {row$pvalue}
  ')
  call_neo4j(query, con)
})

p_log_list <- dados_ml %>%
  transmute(name = `DISEASE/TRAIT`, p_log = p_log) %>%
  purrr::transpose()

walk(p_log_list, function(row) {
  query <- glue('
    MATCH (d:Disease_Trait {{name: "{row$name}"}})
    SET d.p_log = {row$p_log}
  ')
  call_neo4j(query, con)
})


