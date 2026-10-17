# Setup required bioinformatics and visualization packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install Bioconductor packages for RNA-seq analysis
BiocManager::install(c("DESeq2", "GEOquery", "EnhancedVolcano"), update = FALSE)

# Install standard CRAN packages for visualization and data manipulation
install.packages(c("tidyverse", "pheatmap"))

cat("All bioinformatics dependencies installed successfully!\n")
