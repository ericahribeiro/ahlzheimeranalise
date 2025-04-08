# Carregar pacotes necessários
library(dplyr)

# 🔹 1. Encontrar o SNP mais frequente no dataset filtrado (Alzheimer)
top_snp <- batch_data %>%
  count(SNPS, sort = TRUE) %>%
  slice(1) %>%
  pull(SNPS)

cat("🔍 SNP mais frequente no Alzheimer: ", top_snp, "\n")

# 🔹 2. Filtrar o SNP no catálogo GWAS completo
snp_completo <- gwas %>%
  filter(SNPS == top_snp) %>%
  select(SNPS, `DISEASE/TRAIT`, `P-VALUE`, `OR or BETA`, `MAPPED_GENE`)

cat("\n📊 Total de registros encontrados no GWAS para esse SNP:", nrow(snp_completo), "\n")

# 🔹 3. Mostrar resumo estatístico
cat("\n📈 Resumo dos P-VALUES:\n")
print(summary(snp_completo$`P-VALUE`))

cat("\n📈 Resumo dos OR or BETA:\n")
print(summary(snp_completo$`OR or BETA`))

# 🔹 4. Ver todas as doenças associadas a esse SNP
cat("\n🧠 Doenças associadas a esse SNP:\n")
print(unique(snp_completo$`DISEASE/TRAIT`))

# 🔹 5. Ver os genes mapeados para esse SNP
cat("\n🧬 Genes mapeados para esse SNP:\n")
print(unique(snp_completo$`MAPPED_GENE`))

# (Opcional) Visualizar os dados em tabela
# View(snp_completo)
