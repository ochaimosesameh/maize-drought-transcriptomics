# ==============================================================================
# Title: Maize (Zea mays L.) Drought Stress Transcriptomic Analysis
# Author: Ameh Ochai Moses
# Description: Analysis of RNA-seq expression data identifying upregulated 
#              antioxidant enzyme genes during severe water deficit.
# ==============================================================================

library(DESeq2)
library(tidyverse)
library(pheatmap)
library(EnhancedVolcano)

# ------------------------------------------------------------------------------
# STEP 1: Create Simulated RNA-seq Gene Expression Dataset
# (Models real NCBI GEO transcriptomic data for maize leaves: Control vs Drought)
# ------------------------------------------------------------------------------
set.seed(42)

# Define 20 genes (including key antioxidant enzyme gene markers)
gene_names <- c(
  "ZmCAT1", "ZmCAT2", "ZmSOD1", "ZmSOD2", "ZmAPX1", "ZmAPX2", "ZmPPO1", "ZmPPO2",
  "ZmGST1", "ZmGST2", "ZmGR1", "ZmPRX1", "ZmACT1", "ZmTUB1", "ZmGAPDH",
  "ZmMYB1", "ZmWRKY1", "ZmHSP70", "ZmLEA1", "ZmP5CS1"
)

# 6 biological samples: 3 Well-Watered (Control), 3 Water Deficit (Drought)
samples <- c("Control_1", "Control_2", "Control_3", "Drought_1", "Drought_2", "Drought_3")

# Generate synthetic raw count data
counts_matrix <- matrix(
  rpois(120, lambda = 500), 
  nrow = 20, 
  ncol = 6, 
  dimnames = list(gene_names, samples)
)

# Induce strong upregulation in antioxidant genes under drought conditions
antioxidant_genes <- c("ZmCAT1", "ZmCAT2", "ZmSOD1", "ZmAPX1", "ZmPPO1", "ZmLEA1")
counts_matrix[antioxidant_genes, 4:6] <- counts_matrix[antioxidant_genes, 4:6] * 4

# Define sample metadata
col_data <- data.frame(
  condition = factor(c("Control", "Control", "Control", "Drought", "Drought", "Drought")),
  row.names = samples
)

# ------------------------------------------------------------------------------
# STEP 2: Differential Expression Analysis with DESeq2
# ------------------------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = counts_matrix,
  colData = col_data,
  design = ~ condition
)

# Run DESeq2 statistical pipeline
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "Drought", "Control"))

# Convert results to clean data frame
res_df <- as.data.frame(res) %>%
  rownames_to_column(var = "Gene") %>%
  arrange(padj)

print("Top Differentially Expressed Genes (Drought vs Control):")
print(head(res_df))

# Save output table
if(!dir.exists("results")) dir.create("results")
write.csv(res_df, "results/differential_expression_results.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# STEP 3: Data Visualization (Volcano Plot & Heatmap)
# ------------------------------------------------------------------------------

# 1. Volcano Plot
pdf("results/volcano_plot.pdf", width = 8, height = 6)
EnhancedVolcano(
  res_df,
  lab = res_df$Gene,
  x = 'log2FoldChange',
  y = 'pvalue',
  pCutoff = 0.05,
  FCcutoff = 1.5,
  pointSize = 3.0,
  labSize = 4.0,
  title = 'Maize Leaf Gene Expression: Drought vs Control',
  subtitle = 'Highlighting Antioxidant Response Genes',
  legendPosition = 'right'
)
dev.off()

# 2. Heatmap of Antioxidant Genes
norm_counts <- counts(dds, normalized = TRUE)
antioxidant_counts <- norm_counts[antioxidant_genes, ]

pdf("results/antioxidant_heatmap.pdf", width = 7, height = 5)
pheatmap(
  log2(antioxidant_counts + 1),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_col = col_data["condition"],
  main = "Expression Profile of Key Antioxidant Genes under Water Stress"
)
dev.off()

cat("Analysis complete! Plots and results saved to the 'results/' folder.\n")
