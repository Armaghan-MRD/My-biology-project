"""
Bar plot visualization of gene expression
-----------------------------------------

This script loads a CSV file, filters for hub genes,
and creates a bar plot of expression for a selected sample.

Author: Armaghan_MRD
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# -------------------------------
# Load data
# -------------------------------
# Replace "X" with the path to your CSV file
stagei = pd.read_csv("X")

# First column as row names (like row.names in R)
stagei.index = stagei.iloc[:, 0]
stagei = stagei.iloc[:, 1:]   # drop the first column

# -------------------------------
# Filter for hub genes
# -------------------------------
# Replace "X" with your hub gene names
hubgenes = ["X"]
stagei = stagei.loc[hubgenes]

# Transpose (samples as rows, genes as columns)
stagei = stagei.T
stagei["sample"] = stagei.index
stagei["sample"] = stagei["sample"].str.slice(0, 15)

stagei.index = stagei["sample"]
stagei = stagei.drop(columns=["sample"])

# Transpose back (genes as rows)
stagei = stagei.T
stagei["Gene"] = stagei.index

# -------------------------------
# Select one sample for barplot
# -------------------------------
# Replace "X" with your sample ID
sample_id = "X"
y = stagei[["Gene", sample_id]]

print("Available samples:", stagei.columns.tolist())

# -------------------------------
# Bar plot
# -------------------------------
plt.figure(figsize=(10, 5))
sns.barplot(x="Gene", y=sample_id, data=y, palette="viridis")

plt.xticks(rotation=45, ha="right", fontsize=8)
plt.ylabel(f"Expression of genes in {sample_id} (log2)")
plt.xlabel("Gene")
plt.tight_layout()

# -------------------------------
# Save plot
# -------------------------------
save_dir = r"E:\R\RStudio\R_CLASS\Plots"   # change if needed
os.makedirs(save_dir, exist_ok=True)

plt.savefig(os.path.join(save_dir, "arma_first_bar_chart_in_python.pdf"))
plt.show()
