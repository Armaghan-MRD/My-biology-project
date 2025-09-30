# GSE116256_DEG
**limma-voom bulk DEG analysis (summed single-cell counts per sample)**

## Overview
This repository contains an R script that performs differential expression analysis between **Control** and **Disease** samples using a limma-voom pipeline on summed single-cell counts (bulk-like). The analysis saves full DE results and a table of significant DEGs.

## Files
- `GSE116256_DEG_analysis.R` — main R script (reads per-sample `.dem.txt` files, builds DGEList, voom, limma, saves results).
- `DEG_Disease_vs_Control.csv` — (output) full DE results (created when you run the script).
- `DEG_Disease_vs_Control_sig.csv` — (output) significant DEGs filtered by `adj.P.Val < 0.05` and `|logFC| > 1`.

## Requirements
- R (>= 4.0 recommended)
- R packages: `edgeR`, `limma`, `dplyr`, `readr`, `tibble`
  ```r
  install.packages(c("dplyr","readr","tibble"))
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install(c("edgeR","limma"))
