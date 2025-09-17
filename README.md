# My-project
data analysis 
Bar plots
# Load data Python code in Spyder

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# -------------------------------
# Load data Python code in Spyder
# -------------------------------
stagei = pd.read_csv("X")

# First column as row names (like row.names in R)
stagei.index = stagei.iloc[:, 0]
stagei = stagei.iloc[:, 1:]   # drop the first column

# Hub genes list
hubgenes = ["X"]

# Keep only hubgenes (rows)
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
sample_id = "TCGA-E2-A1L7-11"
y = stagei[["Gene", sample_id]]
print(stagei.columns.tolist())

# -------------------------------
# Bar plots
# -------------------------------
plt.figure(figsize=(10,5))
sns.barplot(x="Gene", y=sample_id, data=y, palette="viridis")

plt.xticks(rotation=45, ha="right", fontsize=8)
plt.ylabel(f"Expression of genes in {sample_id} (log2)")
plt.xlabel("Gene")
plt.tight_layout()

# Save to PDF
import os
import matplotlib.pyplot as plt

# Folder to save into
save_dir = r"E:\R\RStudio\R CLASS\Plots"
os.makedirs(save_dir, exist_ok=True)

# Save and show
plt.savefig(os.path.join(save_dir, "arma_first_bar_chart_in_python.pdf"))
plt.show()
