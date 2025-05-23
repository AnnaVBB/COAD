# Pacotes necessários
packages <- c("tidyverse", "DESeq2", "pheatmap", "EnhancedVolcano", "RColorBrewer", "ggplot2", "readr")
installed <- rownames(installed.packages())
lapply(setdiff(packages, installed), install.packages)
lapply(packages, library, character.only = TRUE)
library(purrr)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
library(DESeq2)

manifesto <- read.delim("C:/Users/Usuário/Documents/Projetos em R/Bioinformática/COAD/manifest_coad_filtered.txt", check.names = FALSE)
map_df <- read.delim("C:/Users/Usuário/Documents/Projetos em R/Bioinformática/COAD/gdc_manifest_coad_one.txt.map2submitterID.tsv")
sample_ids <- colnames(raw_counts)[colnames(raw_counts) != "gene_id"]   #Extrai os sample IDs da matriz

#Subset do mapa para os sample_ids da matriz de contagem
samples_metadata <- map_df[match(sample_ids, map_df$file_name), c("file_name", "submitter_id", "submitter_sample")]
# Adiciona a condição com base no sufixo do código da amostra (Normal vs Tumor)
samples_metadata$condition <- ifelse(grepl("-11[A-Z]?$", samples_metadata$submitter_sample), "Normal", "Tumor")

# Ler os dados clínicos brutos e limpar
clinical_data <- read.delim("C:/Users/Usuário/Documents/Projetos em R/Bioinformática/COAD/clinical.tsv")
clinical_data <- na.omit(clinical_data[, c(
  "cases.case_id",
  "cases.submitter_id",
  "demographic.age_at_index",
  "demographic.gender",
  "demographic.vital_status",
  "diagnoses.ajcc_pathologic_stage"
)])

# Juntar samples_metadata com clinical_dt via submitter_id
colData_full <- merge(
  samples_metadata,
  clinical_dt,
  by.x = "submitter_id",
  by.y = "cases.submitter_id",
  all.x = TRUE
)

# Garantir a ordem das amostras igual à da matriz
colData_full <- colData_full[match(sample_ids, colData_full$file_name), ]
rownames(colData_full) <- colData_full$file_name

stopifnot(all(rownames(colData_full) == sample_ids))

# Garantir que 'condition' está presente (build_matrix.R)
if (!"condition" %in% colnames(colData_full)) {
  stop("A coluna 'condition' está ausente de colData_full.")
}

# Transformar 'condition' explicitamente em fator e definir a ordem correta
colData_full$condition <- factor(colData_full$condition, levels = c("Normal", "Tumor"))

# Verificar os níveis
levels(colData_full$condition)

all(colnames(raw_counts[, sample_ids]) == rownames(colData_full))

# Inicializar com o primeiro arquivo
raw_counts <- read.delim(filtered_files[1], skip = 1)
raw_counts <- raw_counts[, c("gene_id", "unstranded")]
colnames(raw_counts)[2] <- basename(filtered_files[1])

# Loop para os demais arquivos
for (f in filtered_files[-1]) {
  tmp <- read.delim(f, skip = 1)
  tmp <- tmp[, c("gene_id", "unstranded")]
  colnames(tmp)[2] <- basename(f)
  raw_counts <- merge(raw_counts, tmp, by = "gene_id")
}

# Verificar se a coluna gene_id ainda existe
if (!"gene_id" %in% colnames(raw_counts)) {
  stop("A coluna 'gene_id' não foi encontrada após o merge. Verifique os arquivos de entrada.")
}

# Salvar os gene_ids antes de remover a coluna
gene_ids <- raw_counts$gene_id

# Atribuir gene_ids como rownames
rownames(raw_counts) <- gene_ids

# Remover coluna gene_id
raw_counts <- raw_counts[, -1]

# Verificar se rownames estão corretos
if (is.null(rownames(raw_counts))) {
  stop("Erro: rownames da matriz 'raw_counts' ainda estão NULL após atribuição.")
}

# Garantir matriz numérica
raw_counts <- as.matrix(raw_counts)
mode(raw_counts) <- "numeric"

dds <- DESeqDataSetFromMatrix(
  countData = raw_counts,
  colData = sample_conditions,
  design = ~ condition
)

# Análise
dds <- DESeq(dds)
vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = "condition")

head(rownames(raw_counts))

# --- Análise com DESeq em build_matrix.R---
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "Tumor", "Normal"))
res <- res[order(res$padj), ]
head(res)
# ----------------------------
##Estatísticas dados clínicos 
### Idade ao diagnóstico (média, mediana, desvio padrão, 1° e 3° quartil) ###
colData_full$demographic.age_at_index <- as.numeric(colData_full$demographic.age_at_index)
summary(colData_full$demographic.age_at_index)

### Gênero (Frequência absoluta/relativa) ###
table(colData_full$`demographic.gender`)
prop.table(table(colData_full$`demographic.gender`)) * 100

### Tecido- Normal vs Tumoral ###
table(colData_full$condition)
prop.table(table(colData_full$condition)) * 100

### Estágio do câncer (distribuição entre os estágios) ###
table(colData_full$`diagnoses.ajcc_pathologic_stage`)
prop.table(table(colData_full$`diagnoses.ajcc_pathologic_stage`)) * 100

stage_colors <- c(
  "--" = "#8da0cb",
  "Stage I" = "#FFD700",
  "Stage IA" = "#FFFF00",
  "Stage II" = "#9ACD32",
  "Stage IIA" = "#ADFF2F",
  "Stage IIB" = "#3CB371",
  "Stage IIC" = "#008000",
  "Stage III" = "#87CEEB",
  "Stage IIIA" = "#1E90FF",
  "Stage IIIB" = "#4682B4",
  "Stage IIIC" = "#4169E1",
  "Stage IV" = "#FF69B4",
  "Stage IVA" = "#DA70D6",
  "Stage IVB" = "#C71585"
)
colData_full_stage <- colData_full[!is.na(colData_full$`diagnoses.ajcc_pathologic_stage`), ]
colData_full_stage$diagnoses.ajcc_pathologic_stage <- factor(
  colData_full_stage$`diagnoses.ajcc_pathologic_stage`,
  levels = names(stage_colors)
)

ggplot(colData_full_stage, aes(x = diagnoses.ajcc_pathologic_stage, fill = diagnoses.ajcc_pathologic_stage)) +
  geom_bar() +
  scale_fill_manual(values = stage_colors, name = "Estágio do Câncer") +
  labs(
    title = "Distribuição de Estágios Patológicos do Câncer",
    x = "Estágio",
    y = "Frequência Absoluta"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

### Status vital (proporções) ###
table(colData_full$`demographic.vital_status`)
prop.table(table(colData_full$`demographic.vital_status`)) * 100

##Estatísticas inferenciais 
### Status vital x Gênero x Idade ###
ggplot(colData_full, aes(x = demographic.age_at_index,
                         fill = demographic.vital_status)) +
  geom_histogram(binwidth = 5, color = "black", alpha = 0.8) +
  facet_wrap(~ demographic.gender) +
  scale_fill_manual(values = c("Alive" = "#66c2a5", "Dead" = "#fc8d62")) +
  labs(title = "Distribuição da idade por status vital e gênero",
       x = "Idade", y = "Frequência") +
  theme_minimal()

# Tecido × Estágio do câncer
ggplot(colData_full, aes(x = diagnoses.ajcc_pathologic_stage, fill = condition)) +
  geom_bar(position = "dodge") +
  scale_fill_manual(values = c("Normal" = "#8da0cb", "Tumor" = "#e78ac3")) +
  labs(title = "Estágio do câncer por tipo de tecido",
       x = "Estágio AJCC", y = "Frequência") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Tecido × Status vital
ggplot(colData_full, aes(x = demographic.vital_status, fill = condition)) +
  geom_bar(position = "dodge") +
  scale_fill_manual(values = c("Normal" = "#8da0cb", "Tumor" = "#e78ac3")) +
  labs(title = "Status vital por tipo de tecido",
       x = "Status Vital", y = "Frequência") +
  theme_minimal()


# ----------------------------
### Estatísticas dados clínicos + expressão gênica
head(colnames(vsd))
head(colData_full$file_name)
colData_full <- colData_full[match(colnames(vsd), colData_full$file_name), ]
all(colData_full$file_name == colnames(vsd))  # Deve retornar TRUE

colnames(colData_full)[colnames(colData_full) == "demographic.age_at_index"] <- "age"
colnames(colData_full)[colnames(colData_full) == "demographic.gender"] <- "gender"
colnames(colData_full)[colnames(colData_full) == "demographic.vital_status"] <- "vital_status"
colnames(colData_full)[colnames(colData_full) == "diagnoses.ajcc_pathologic_stage"] <- "tumor_stage"

colData(vsd)$condition <- factor(colData_full$condition)
colData(vsd)$age <- as.numeric(colData_full$age)
colData(vsd)$gender <- factor(colData_full$gender)
colData(vsd)$vital_status <- factor(colData_full$vital_status)
colData(vsd)$tumor_stage <- factor(colData_full$tumor_stage)

colnames(colData(vsd))  # Deve listar: condition, age, gender, vital_status, tumor_stage
summary(colData(vsd))   # Estatísticas descritivas
##############Genes mais expressos
gene_means <- rowMeans(assay(vsd))
top20_genes <- names(sort(gene_means, decreasing = TRUE)[1:20])

# Barplot com as expressões médias
barplot(gene_means[top20_genes],
        las = 2,           # rotação dos rótulos
        col = "steelblue",
        main = "Top 20 genes mais expressos (média)",
        ylab = "Expressão média (log)",
        cex.names = 0.5)

########
gene_vars <- apply(assay(vsd), 1, var)
top20_genes <- names(sort(gene_vars, decreasing = TRUE))[1:20]

# Plotando idade × expressão para os top 20
par(mfrow = c(4, 5))  # grid 4 linhas × 5 colunas

#por conta do tamanho, gerar um pdf 
pdf("C:/Users/Usuário/Documents/Projetos em R/idade_vs_expressao_top20.pdf", width = 14, height = 10)  # salva em PDF grande 
par(mfrow = c(4, 5))  # 4 linhas × 5 colunas
for (gene in top20_genes) {
  expr <- assay(vsd)[gene, ]
  idade <- colData(vsd)$age
  plot(idade, expr,
       main = gene, xlab = "Idade", ylab = "Expressão (log)",
       pch = 19, col = "darkblue")
  abline(lm(expr ~ idade), col = "red")
}
dev.off()

#PCA 
plotPCA(vsd, intgroup = "condition")
plotPCA(vsd, intgroup = "tumor_stage")
plotPCA(vsd, intgroup = "vital_status")
print(plotPCA(vsd, intgroup = "tumor_stage"))
print(plotPCA(vsd, intgroup = c("condition", "gender")))
print(plotPCA(vsd, intgroup = c("condition", "vital_status","gender")))


### Estágio do câncer x Expressão gênica DO GENE MAIS EXPRESSO ###
boxplot(assay(vsd)["ENSG00000198804.2", ] ~ colData(vsd)$tumor_stage,
        main = "MT-CO1 por estágio tumoral", ylab = "Expressão (log)", xlab = "Estágio")

anova_res <- aov(assay(vsd)["ENSG00000198804.2", ] ~ colData(vsd)$tumor_stage)
summary(anova_res)

### Status vital x Expressão gênica (a expressão gênica afeta as chances de sobreviver?) ###
boxplot(assay(vsd)["ENSG00000198804.2", ] ~ colData(vsd)$vital_status,
        main = "MT-CO1 por status vital", ylab = "Expressão (log)", xlab = "Status Vital")
t.test(assay(vsd)["ENSG00000198804.2", ] ~ colData(vsd)$vital_status)

### Sobrevivência (kaplan-Meier)
library(survival)
library(survminer)

# Supondo que você tenha colData(vsd)$days_to_death ou $days_to_last_follow_up
surv_time <- ifelse(colData(vsd)$vital_status == "Dead", 
                    colData(vsd)$days_to_death, 
                    colData(vsd)$days_to_last_follow_up)

surv_obj <- Surv(time = surv_time, event = colData(vsd)$vital_status == "Dead")

fit <- survfit(surv_obj ~ colData(vsd)$condition)
ggsurvplot(fit, data = colData(vsd), pval = TRUE, risk.table = TRUE)

