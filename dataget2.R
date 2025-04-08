# Instalar pacotes se necessário
if (!require("dplyr")) install.packages("dplyr")
if (!require("readr")) install.packages("readr")
if (!require("neo4r")) install.packages("neo4r")
if (!require("ggraph")) install.packages("ggraph")
if (!require("tidygraph")) install.packages("tidygraph")

library(dplyr)
library(readr)
library(neo4r)
library(ggraph)
library(tidygraph)
library(glue)
library(purrr)


url <- "https://www.ebi.ac.uk/gwas/api/search/downloads/full"
dest <- "~/dados/gwas_catalog.tsv"
download.file(url, destfile = dest)


gwas <- read_tsv("/home/rstudio/dados/gwas_catalog.tsv", show_col_types = FALSE)

# 2. Ler e filtrar dados de Alzheimer
alzheimer <- gwas %>% filter(grepl("Alzheimer", `DISEASE/TRAIT`, ignore.case = TRUE))

cat("Total de registros sobre Alzheimer: ", nrow(alzheimer), "\n")

# 3. Análise simples
cat("SNPs únicos: ", length(unique(alzheimer$SNPS)), "\n")
cat("Top 10 genes associados:\n")
print(sort(table(alzheimer$`MAPPED_GENE`), decreasing = TRUE)[1:10])

# 4. Conectar ao Neo4j
con <- neo4j_api$new(
  url = "http://neo4j:7474",
  user = "neo4j",
  password = "admin123"  
)

batch_df <- alzheimer %>%
  select(SNPS, `MAPPED_GENE`, `P-VALUE`, `OR or BETA`, CHR_ID, CHR_POS) %>%
  filter(!is.na(SNPS), !is.na(`MAPPED_GENE`)) %>%
  mutate(
    SNPS = gsub("\"", "", SNPS),
    MAPPED_GENE = toupper(gsub("\"", "", `MAPPED_GENE`)),
    `P-VALUE` = as.numeric(`P-VALUE`),
    `OR or BETA` = as.numeric(`OR or BETA`),
    CHR_ID = as.character(CHR_ID),
    CHR_POS = as.character(CHR_POS)
  )



batch_data <- batch_df %>%
  transmute(
    snp = SNPS,
    gene = `MAPPED_GENE`,
    pval = `P-VALUE`,
    or = `OR or BETA`,
    chr = CHR_ID,
    pos = CHR_POS
  ) %>%
  purrr::transpose()

# Pegamos um subconjunto para teste, com 100 registros por exemplo
batch_data <- batch_df

# Usar os dados direto (sem sobrescrever nada depois)
batch_list <- batch_df %>%
  pmap(function(SNPS, MAPPED_GENE, `P-VALUE`, `OR or BETA`, CHR_ID, CHR_POS) {
    list(
      snp = SNPS,
      gene = MAPPED_GENE,
      pval = `P-VALUE`,
      or = `OR or BETA`,
      chr = CHR_ID,
      pos = CHR_POS
    )
  })

cypher <- "
UNWIND $batch AS row
MERGE (d:Disease_Trait {name: 'Alzheimer'})
MERGE (g:Gene {name: row.gene})
MERGE (s:SNP {id: row.snp})
SET s.pvalue = row.pval,
    s.effect_size = row.or,
    s.chr = row.chr,
    s.position = row.pos
MERGE (s)-[:ASSOCIATED_WITH]->(d)
MERGE (s)-[:LOCATED_IN]->(g)
"

library(glue)

walk(batch_list, function(row) {
  query <- glue('
    MERGE (d:Disease_Trait {{name: "Alzheimer"}})
    MERGE (g:Gene {{name: "{row$gene}"}})
    MERGE (s:SNP {{id: "{row$snp}"}})
    SET s.pvalue = {row$pval},
        s.effect_size = {row$or},
        s.chr = "{row$chr}",
        s.position = "{row$pos}"
    MERGE (s)-[:ASSOCIATED_WITH]->(d)
    MERGE (s)-[:LOCATED_IN]->(g)
  ')
  
  call_neo4j(
    query = query,
    con = con
  )
})



results <- call_neo4j(
  "MATCH (s:SNP)-[:ASSOCIATED_WITH]->(d:Disease_Trait) RETURN s, d",
  con,
  type = "graph"
)

library(tidyr)


# Corrigir rótulos
results$nodes$label <- purrr::map_chr(results$nodes$label, ~ .x[[1]])

# (opcional) Limpar colunas que causam conflito
results$nodes <- results$nodes[, c("id...1", "label", "name")]

# Usar nodes e rels diretamente
graph <- tbl_graph(nodes = results$nodes, edges = results$rels, directed = TRUE)

library(ggraph)

ggraph(graph, layout = "fr") +
  geom_edge_link(alpha = 0.5, color = "gray") +
  geom_node_point(aes(color = label), size = 5) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  theme_void() +
  ggtitle("SNPs associados ao Alzheimer")



