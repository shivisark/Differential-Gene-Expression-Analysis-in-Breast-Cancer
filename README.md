# TCGA-BRCA Differential Expression Analysis

This project performs differential gene expression analysis on TCGA Breast Cancer (BRCA) data using RNA-Seq count data from tumor and matched normal tissue samples. It identifies the top 20 significantly altered genes and visualizes the data using PCA, heatmaps, and clustering.

---

## Installation

Install the necessary packages:

```r
# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "TCGAbiolinks", "DESeq2", "SummarizedExperiment",
  "EnhancedVolcano", "pheatmap", "biomaRt"
))

# CRAN packages
install.packages(c("dplyr", "ggplot2"))
