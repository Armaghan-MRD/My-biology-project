# GSE154109 Pseudobulk DESeq2 Analysis

**Pseudo-bulk differential expression analysis** (Disease vs Control)  
This script uses **Seurat v5** to aggregate single-cell 10x Genomics data into pseudo-bulk counts per sample, and then performs **DESeq2** analysis.


## Overview
- Input: 10x Genomics matrices from GSE154109 (stored in `Control/` and `Disease/` subfolders).
- Processing:
  1. Load each sample as a Seurat object.
  2. Perform QC filtering (min/max features, mitochondrial percentage).
  3. Merge samples into one Seurat object.
  4. Aggregate counts by sample (pseudo-bulk).
  5. Run DESeq2 on pseudo-bulk counts.
- Output:
  - `DEG_Disease_vs_Control.csv` → all genes with DESeq2 results.
  - `DEG_Disease_vs_Control_sig.csv` → significant genes (`padj < 0.05` and `|log2FC| > 1`).

---

## Requirements
- R (>= 4.0)
- Packages:
  ```r
  install.packages(c("Seurat", "Matrix", "dplyr", "readr", "tibble"))
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("DESeq2")
