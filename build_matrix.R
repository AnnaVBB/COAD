library(dplyr)

# Carregar manifesto
manifesto <- read.delim("C:\\Users\\Usuário\\Documents\\Projetos em R\\Bioinformática\\COAD\\manifest_coad_filtered.txt", header = TRUE)
manifest_ids <- manifesto$file_name  # Ajuste para a coluna que contém os IDs no manifesto
print(manifest_ids)  # Verifique os IDs do manifesto

# Listar todos os arquivos .tsv diretamente
rnaseq_files <- list.files("C:\\Users\\Usuário\\Documents\\Projetos em R\\Bioinformática\\COAD\\tsv", pattern = ".tsv$", full.names = TRUE)

# Filtrar os arquivos que estão no manifesto
filtered_files <- rnaseq_files[basename(rnaseq_files) %in% manifest_ids]

# Inicializar com o primeiro
raw_counts <- read.delim(filtered_files[1], skip = 1)
raw_counts <- raw_counts[, c("gene_id", "unstranded")]
colnames(raw_counts)[2] <- basename(filtered_files[1])

# Loop para os demais
for (f in filtered_files[-1]) {
  tmp <- read.delim(f, skip = 1)
  tmp <- tmp[, c("gene_id", "unstranded")]
  colnames(tmp)[2] <- basename(f)
  raw_counts <- merge(raw_counts, tmp, by = "gene_id")
}

# Ajustar nomes das linhas e converter para matriz numérica
rownames(raw_counts) <- raw_counts$gene_id
raw_counts <- raw_counts[, -1]

# Transformar matriz para numérica
raw_counts <- as.matrix(raw_counts)
raw_counts <- apply(raw_counts, 2, as.numeric)

# Verificar se a matriz é numérica
if (!is.numeric(raw_counts)) {
  stop("Erro: Valores na matriz 'raw_counts' não são numéricos.")
}

# Salvar matriz gerada
saveRDS(raw_counts, file = "C:\\Users\\Usuário\\Documents\\Projetos em R\\Bioinformática\\COAD\\TCGA_COAD_RNASEQ_ORDERED.rds")

# Verificar dimensões da matriz
print(dim(raw_counts))
print(head(raw_counts))
##################################
manifesto <- read.delim("C:\\Users\\Usuário\\Documents\\Projetos em R\\Bioinformática\\COAD\\manifest_coad_filtered.txt", header = TRUE, check.names= FALSE)
matrix_ids <- setdiff(colnames(raw_counts), "gene_id")
uuid_to_tcga <- manifesto[, c("file_name", "submitter_sample")]
metadata_ids <- uuid_to_tcga$file_name

missing_ids <- setdiff(matrix_ids, metadata_ids)

if (length(missing_ids) > 0) {
  warning("IDs ausentes nos metadados: ", paste(missing_ids, collapse = ", "))
  stop("Erro: Existem IDs na matriz de contagem que não têm correspondência nos metadados.")
}

#se quiser verificar quais matrix_ids deram match:
for (i in seq_along(matrix_ids)) {
  file <- matrix_ids[i]
  match_idx <- which(manifesto$file_name == file)
  if (length(match_idx) == 0) {
    cat("❌ Sem match exato para:", shQuote(file), "\n")
  } else {
    cat("✅ Match com:", shQuote(file), "\n")
  }
}

# Adicionar condição com base no código de tipo de tecido no submitter_sample
# Sufixos "-11A", "-11B" indicam amostras normais; tudo o mais é considerado tumoral
uuid_to_tcga$condition <- ifelse(
  grepl("-11[A-Z]?$", uuid_to_tcga$submitter_sample), "Normal", "Tumor"
)

# Verificar o mapeamento com condições
head(uuid_to_tcga)

# Criar o dataframe de condições somente com IDs correspondentes
sample_conditions <- data.frame(
  row.names = matrix_ids,
  condition = uuid_to_tcga$condition[match(matrix_ids, uuid_to_tcga$file_name)]
)
head(sample_conditions)
# Verificar se todas as amostras têm condição associada
if (any(is.na(sample_conditions$condition))) {
  stop("Erro: Algumas amostras não têm condição associada. Verifique o mapeamento.")
}
print(missing_ids)



# --- Processar Matriz de Contagem ---
# Remover coluna de IDs dos genes, se presente
if ("gene_id" %in% colnames(raw_counts)) {
  row.names(raw_counts) <- raw_counts$gene_id
  raw_counts <- raw_counts[, -1]
}

# Garantir que a matriz é numérica
raw_counts <- as.matrix(raw_counts)
mode(raw_counts) <- "numeric"

# Verificar se há valores negativos
if (any(raw_counts < 0)) {
  warning("A matriz contém valores negativos. Substituindo por 0.")
  raw_counts[raw_counts < 0] <- 0
}

# Remover genes com valores ausentes
raw_counts <- raw_counts[complete.cases(raw_counts), ]

# Verificar estrutura final da matriz de contagem
summary(as.vector(raw_counts))
all(raw_counts == floor(raw_counts))  # Deve retornar TRUE

# --- Verificar Dados ---
print(dim(raw_counts))
print(table(sample_conditions$condition))

# Converter condição para fator
sample_conditions$condition <- as.factor(sample_conditions$condition)


# --- Criar Objeto DESeqDataSet ---
dds <- DESeqDataSetFromMatrix(
  countData = raw_counts,
  colData = sample_conditions,
  design = ~ condition
)
#Verificar se os nomes batem (deve retornar TRUE): 
all(colnames(raw_counts) == rownames(sample_conditions)) 

# --- Análise com DESeq ---
dds <- DESeq(dds)

# Extrair resultados.  Ordenar por p-valor ajustado.Visualizar os primeiros genes
res <- results(dds, contrast = c("condition", "Tumor", "Normal"))
res <- res[order(res$padj), ]
head(res)

#saveRDS(dds, file = 'C:/Users/Usuário/Downloads/TCGA/TCGA_COAD_DESeq.rds')

res_sig <- res[!is.na(res$padj) & res$padj < 0.05 & abs(res$log2FoldChange) > 1, ]
dim(res_sig)  # número de genes significativos
head(res_sig[order(res_sig$padj), ])  # top genes

#MA plot (expressão x log2FC)
plotMA(res, ylim = c(-5, 5))

#PCA (agrupamento das amostras) 
vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = "condition")

#Heatmap dos 30 genes mais expressos
library(pheatmap)
res_clean <- res[!is.na(res$padj), ] #remover os genes com padj NA
top_genes <- head(order(res$padj), 30)
mat_raw <- assay(vsd)[top_genes, ] # Matriz normalizada (sem escalonamento ainda)
gene_sd <- apply(mat_raw, 1, sd)# Calcular desvio padrão por gene
zero_sd_genes <- names(gene_sd[gene_sd == 0])# Verificar quais são 0 (causa de NaN)
print(zero_sd_genes)
mat_filtered <- mat_raw[gene_sd > 0, ]# Filtrar genes com variação
mat <- t(scale(t(mat_filtered)))# Escalonar com segurança
pheatmap(mat,
         annotation_col = annotation_col,
         show_rownames = TRUE,
         show_colnames = TRUE,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         fontsize_row = 8)

#write.csv(as.data.frame(res), file = "C://Users//Usuário//Downloads//DESeq2_all_results.csv")
#write.csv(as.data.frame(res_sig), file = "DESeq2_significant_results.csv")

#Estatísticas descritivas
# Contagem total de reads por amostra (soma das contagens por coluna)
sample_sums <- colSums(counts(dds))

summary(sample_sums)  # Mostra min, 1Q, mediana, média, etc.
boxplot(sample_sums, main = "Total de reads por amostra", ylab = "Total counts")

# Média e mediana por gene
gene_means <- rowMeans(counts(dds))
gene_medians <- apply(counts(dds), 1, median)

summary(gene_means)
 
