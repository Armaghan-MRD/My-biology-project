This R script performs a Robust Rank Aggregation (RRA) meta-analysis across multiple Neuroblastoma differential expression datasets.
It integrates up/down-regulated gene rankings, applies auto-correlation-based direction flipping, and outputs unified RRA statistics and per-dataset metrics.

Features

Automatic parsing of DEG text files (.txt, tab/space delimited)

Rank aggregation (up & down) using RobustRankAggreg

Auto-flip of datasets to improve mean pairwise correlation

Merged output of logFC, adj.P.Val, and consensus metrics

Exports final Excel file with RRA_up and RRA_down sheets

Input

Place your DEG tables inside:

E:/R/RStudio/R CLASS/Data/Text file for RRA/Neuroblastoma


Each file should contain columns such as:

Gene.symbol | logFC | adj.P.Val

Run
source("Neuroblastoma_RRA_pipeline.R")

Output

Neuroblastoma_RRA_results_final.xlsx — final RRA results

RRA_up: aggregated up-regulated genes

RRA_down: aggregated down-regulated genes

Dependencies

Install these R packages:

install.packages(c("dplyr","readr","RobustRankAggreg","openxlsx"))
