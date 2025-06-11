## ----------------------------------------
## Installing Packages
## ----------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

BiocManager::install(c(
  "TCGAbiolinks", "DESeq2", "SummarizedExperiment",
  "EnhancedVolcano", "pheatmap", "biomaRt"
))
install.packages(c("dplyr", "ggplot2"))

## ----------------------------------------
## Loading Libraries
## ----------------------------------------

library(TCGAbiolinks)
library(DESeq2)
library(SummarizedExperiment)
library(EnhancedVolcano)
library(pheatmap)
library(dplyr)
library(ggplot2)
library(biomaRt)

## ---------------------------------------------------------------------
## Query TCGA ( The Cancer Genomic Atlas) for BRCA Gene Expression Data
## ---------------------------------------------------------------------

query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor", "Solid Tissue Normal")
)

## --------------------------------------------------------
## 113 patients with both tumor tissues and normal tissues
## --------------------------------------------------------

meta <- query$results[[1]]
meta$patient_id <- substr(meta$cases, 1, 12)

## -------------------------------------------
## Identify 113 Matched Tumor-Normal Pairs
## -------------------------------------------

paired_patients <- meta %>%
  group_by(patient_id) %>%
  summarise(types = paste(unique(sample_type), collapse = ",")) %>%
  filter(grepl("Primary Tumor", types) & grepl("Solid Tissue Normal", types)) %>%
  pull(patient_id)

paired_patient_ids <- paired_patients[1:113]
meta_matched <- meta %>% filter(patient_id %in% paired_patient_ids)

#Extracting one tumor sample and one normal sample per patient

tumor_samples <- meta_matched %>%
  filter(sample_type == "Primary Tumor") %>%
  group_by(patient_id) %>%
  slice(1)

normal_samples <- meta_matched %>%
  filter(sample_type == "Solid Tissue Normal") %>%
  group_by(patient_id) %>%
  slice(1)

meta_matched_clean <- bind_rows(tumor_samples, normal_samples)
barcodes_matched <- meta_matched_clean$cases

## ----------------------------------------
## Download and Prepare Gene Expression Data
## ----------------------------------------

query_final <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  barcode = barcodes_matched
)

GDCdownload(query_final)
data <- GDCprepare(query_final)

## ----------------------------------------
## Labeling samples as normal or tumor data
## ----------------------------------------

meta <- as.data.frame(colData(data))
meta$condition <- ifelse(meta$sample_type == "Solid Tissue Normal", "Normal", "Tumor")
meta$condition <- factor(meta$condition, levels = c("Normal", "Tumor"))
colData(data)$condition <- meta$condition

## ----------------------------------------
## Running Differential Expression Analysis (DESeq2)
## ----------------------------------------
dds <- DESeqDataSetFromMatrix(countData = assay(data),
                              colData = colData(data),
                              design = ~ condition)

dds <- dds[rowSums(counts(dds)) >= 10, ]  # removing low-expression genes
dds <- DESeq(dds)

res <- lfcShrink(dds, coef = "condition_Tumor_vs_Normal", type = "apeglm")
resOrdered <- res[order(res$padj), ]
write.csv(as.data.frame(resOrdered), "DESeq2_TCGA_BRCA_results.csv")

## ----------------------------------------
## PCA Plot
## ----------------------------------------

vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = "condition") + ggtitle("PCA: Tumor vs Normal")

## ----------------------------------------
## Heatmap: Top 20 Genes with Gene Symbols
## ----------------------------------------

# Select top 20 most significant genes
topGeneIndices <- head(order(resOrdered$padj), 20)
vsd_mat <- assay(vsd)[topGeneIndices, ]

# Map Ensembl IDs to gene symbols
mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
clean_ids <- gsub("\\..*", "", rownames(vsd_mat))
genes <- getBM(
  filters = "ensembl_gene_id",
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  values = clean_ids,
  mart = mart
)
matched <- match(clean_ids, genes$ensembl_gene_id)
gene_symbols <- ifelse(
  is.na(matched) | genes$hgnc_symbol[matched] == "",
  clean_ids,
  genes$hgnc_symbol[matched]
)
rownames(vsd_mat) <- gene_symbols

# Annotate sample conditions (tumor vs normal)
annotation_col <- data.frame(Condition = colData(vsd)$condition)
rownames(annotation_col) <- colnames(vsd_mat)

# Generate heatmap (sample IDs hidden for clarity)
pheatmap(vsd_mat,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         show_colnames = FALSE,  # <--- Hide messy bottom labels
         fontsize_row = 10,
         annotation_col = annotation_col,
         main = "Top 20 Differentially Expressed Genes: Tumor vs Normal")

# Step 1: Select top 20 genes
topGeneIndices <- head(order(resOrdered$padj), 20)
vsd_mat <- assay(vsd)[topGeneIndices, ]

# Step 2: Replace Ensembl IDs with gene symbols
library(biomaRt)
mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
clean_ids <- gsub("\\..*", "", rownames(vsd_mat))
genes <- getBM(
  filters = "ensembl_gene_id",
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  values = clean_ids,
  mart = mart
)
matched <- match(clean_ids, genes$ensembl_gene_id)
gene_symbols <- ifelse(
  is.na(matched) | genes$hgnc_symbol[matched] == "",
  clean_ids,
  genes$hgnc_symbol[matched]
)
rownames(vsd_mat) <- gene_symbols

# Step 3: Calculate distances and cluster genes
gene_dist <- dist(vsd_mat)                  # Distance between gene expression profiles
gene_clust <- hclust(gene_dist, method = "complete")  # Hierarchical clustering

# Step 4: Plot the dendrogram
plot(gene_clust,
     main = "Dendrogram of Top 20 Differentially Expressed Genes",
     xlab = "Genes",
     sub = "",
     cex = 0.8,
     hang = -1)




