
# Pseudo-bulk Differential Expression (Disease vs Control, multiple samples)
# Seurat v5 compatible: handles Assay5 multilayer assay
# Input: 10x Genomics matrices for Control and Disease subfolders
# Output:
#   DEG_Disease_vs_Control.csv       (all genes)
#   DEG_Disease_vs_Control_sig.csv   (padj < 0.05 & |log2FC| > 1)
################################################################################

# ---- 0. Install & load packages ----
required <- c("Seurat", "Matrix", "dplyr", "DESeq2", "readr", "tibble")
for (pkg in required) if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)

if (!requireNamespace("DESeq2", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("DESeq2")
}

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(DESeq2)
  library(readr)
  library(tibble)
})

# ---- 1. Define sample directories ----
control_dirs <- c("the/location/of/your/file")

disease_dirs <- c("the/location/of/your/file")

# ---- 2. Load datasets into Seurat ----
seurat_list <- list()

load_sample <- function(subdir, group) {
  sample_id <- basename(subdir)
  message("Loading sample: ", sample_id, " (", group, ")")
  counts <- Read10X(data.dir = subdir)
  so <- CreateSeuratObject(counts = counts, project = group, min.cells = 3, min.features = 200)
  so$group <- group
  so$sample <- sample_id
  seurat_list[[sample_id]] <<- so
}

for (subdir in control_dirs) load_sample(subdir, "Control")
for (subdir in disease_dirs) load_sample(subdir, "Disease")

# ---- 3. QC filtering ----
for (s in names(seurat_list)) {
  so <- seurat_list[[s]]
  so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")
  so <- subset(so, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)
  seurat_list[[s]] <- so
  message("Sample ", s, ": retained ", ncol(so), " cells after QC")
}

# QC summary
qc_summary <- data.frame(
  sample = names(seurat_list),
  group = sapply(seurat_list, function(x) x$group[1]),
  cells_after_QC = sapply(seurat_list, ncol)
)
print(qc_summary)
print(table(qc_summary$group))

# ---- 4. Merge all Seurat objects ----
combined <- merge(seurat_list[[1]], y = seurat_list[-1])

# ---- 5. Pseudo-bulk: sum counts per sample ----
pb_list <- AggregateExpression(combined, group.by = "sample", assays = "RNA", slot = "counts")
pb_mat <- pb_list$RNA

# ---- 6. Metadata for DESeq2 ----

sample_table <- data.frame(
  sample = colnames(pb_mat),
  group  = factor(
    sapply(colnames(pb_mat), function(sid) seurat_list[[sid]]$group[1]),
    levels = c("Control", "Disease")
  ),
  stringsAsFactors = FALSE
)

# make sure rownames match columns in pb_mat
rownames(sample_table) <- sample_table$sample

# check!
all(colnames(pb_mat) == rownames(sample_table))

# ---- 7. Run DESeq2 ----
dds <- DESeqDataSetFromMatrix(countData = pb_mat,
                              colData = sample_table,
                              design = ~ group)

dds <- dds[rowSums(counts(dds)) >= 10, ]   # keep expressed genes

dds <- DESeq(dds)
res <- results(dds, contrast = c("group", "Disease", "Control"))

# ---- 8. Save outputs ----
res_df <- as.data.frame(res) %>%
  tibble::rownames_to_column("gene") %>%
  dplyr::select(gene, log2FoldChange, lfcSE, stat, pvalue, padj)

out_csv <- file.path("E:/R/RStudio/R CLASS/Data/GSE154109", "DEG_Disease_vs_Control.csv")
write.csv(res_df, out_csv, row.names = FALSE)
message("Full DE results saved to: ", out_csv)

sig <- res_df %>%
  filter(!is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1)

sig_csv <- file.path("E:/R/RStudio/R CLASS/Data/GSE154109", "DEG_Disease_vs_Control_sig.csv")
write.csv(sig, sig_csv, row.names = FALSE)
message("Significant DEGs saved to: ", sig_csv)

################################################################################
# End
#Armaghan_MRD
