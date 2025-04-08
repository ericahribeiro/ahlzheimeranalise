library(dplyr)

# ------------------------------------------
# 1. Verificar valores ausentes em cada coluna
# ------------------------------------------
cat("\n🔍 1. Valores ausentes por coluna:\n")
print(colSums(is.na(batch_data)))

# ------------------------------------------
# 2. Verificar estrutura e tipo de dados
# ------------------------------------------
cat("\n🔍 2. Estrutura do dataset:\n")
str(batch_data)

# ------------------------------------------
# 3. Resumo estatístico dos dados
# ------------------------------------------
cat("\n🔍 3. Resumo estatístico:\n")
print(summary(batch_data))

# ------------------------------------------
# 4. Verificar SNPs duplicados
# ------------------------------------------
cat("\n🔍 4. Número de SNPs duplicados:\n")
print(sum(duplicated(batch_data$SNPS)))

# ------------------------------------------
# 5. Checar se todos os P-VALUES estão no intervalo válido [0, 1]
# ------------------------------------------
cat("\n🔍 5. P-VALUES fora do intervalo [0, 1]:\n")
print(batch_data %>% filter(`P-VALUE` < 0 | `P-VALUE` > 1))

# ------------------------------------------
# 6. Ver genes com múltiplos nomes (vírgula ou hífen)
# ------------------------------------------
cat("\n🔍 6. Genes com múltiplos nomes (vírgula ou hífen):\n")
print(batch_data %>% filter(grepl(",", MAPPED_GENE) | grepl("-", MAPPED_GENE)))

# ------------------------------------------
# 7. Verificar caracteres estranhos em nomes de genes
# (qualquer coisa fora de letras maiúsculas, números, vírgulas, hífens ou espaços)
# ------------------------------------------
cat("\n🔍 7. Genes com caracteres incomuns:\n")
print(table(grepl("[^A-Z0-9,\\- ]", batch_data$MAPPED_GENE)))
