# My-Project: Multi-omics Analysis & Visualization

This repository contains scripts for:
1. **Pseudo-bulk Differential Expression (R, Seurat + DESeq2)**
2. **Gene Expression Bar Plots (Python, Matplotlib + Seaborn)**

---

##  Project Structure

- `R_pipeline/` → R script for pseudo-bulk DEG analysis (Disease vs Control)
- `Python_barplot/` → Python script for hub gene bar plots
- `data/` → Example input files (optional, if you want to share demo data)

---

## Usage

### 1. R Pipeline
- Navigate to `R_pipeline/`
- Run `pseudo_bulk_DEG.R` in RStudio
- Outputs:
  - `DEG_Disease_vs_Control.csv`
  - `DEG_Disease_vs_Control_sig.csv`

Details: [R_pipeline/README.md](R_pipeline/README.md)

---

### 2. Python Bar Plots
- Navigate to `Python_barplot/`
- Run `barplot_genes.py` in Spyder/Jupyter
- Outputs:
  - PDF bar plots in your chosen folder

Details: [Python_barplot/README.md](Python_barplot/README.md)

---

## 📦 Dependencies

- **R**: Seurat v5, DESeq2, dplyr, readr, tibble
- **Python**: pandas, matplotlib, seaborn

---

## ✍️ Author
Armaghan_MRD
