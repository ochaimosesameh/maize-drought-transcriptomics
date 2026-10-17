# Maize (*Zea mays* L.) Drought Stress Transcriptomics Pipeline

[![R](https://img.shields.io/badge/Language-R-blue.svg)](https://www.r-project.org/)
[![Bioconductor](https://img.shields.io/badge/Bioconductor-DESeq2-green.svg)](https://bioconductor.org/)

## Overview
This repository provides an R computational pipeline for analyzing transcriptomic responses (RNA-seq) in maize (*Zea mays* L.) leaves subjected to severe water deficit (drought stress). 

Building upon wet-lab biochemical observations of enzymatic antioxidant responses (Catalase, Superoxide Dismutase, Ascorbate Peroxidase, and Polyphenol Oxidase), this dry-lab workflow evaluates the underlying transcriptional mechanisms and gene expression changes regulating abiotic stress tolerance in local crop varieties.

---

## Workflow Steps
1. **Data Ingestion & Preprocessing:** Structuring raw gene counts and sample metadata (Control vs. Drought).
2. **Differential Expression Analysis:** Utilizing **DESeq2** to calculate $\log_2$ fold changes and adjusted p-values ($p_{adj} < 0.05$).
3. **Targeted Biomarker Identification:** Filtering specifically for antioxidant enzyme genes (`ZmCAT1`, `ZmSOD1`, `ZmAPX1`, `ZmPPO1`).
4. **Genomic Visualization:** Generating volcano plots and expression heatmaps using **EnhancedVolcano** and **pheatmap**.

---

## Key Results
* Identified significant transcriptional upregulation of primary antioxidant defense genes (`ZmCAT1`, `ZmSOD1`, `ZmAPX1`) during water deficit conditions.
* Established a reproducible framework for screening drought-resilient gene expression signatures in indigenous tropical crop landraces.

---

## How to Run

1. Clone the repository:
   ```bash
   git clone [https://github.com/ochaimosesameh/maize-drought-transcriptomics.git](https://github.com/ochaimosesameh/maize-drought-transcriptomics.git)
   cd maize-drought-transcriptomics
   ```

2. Install dependencies:
   ```R
   source("requirements.R")
   ```

3. Run the analysis:
   ```R
   source("maize_drought_analysis.R")
   ```

Outputs will be saved automatically in the `results/` directory.

---
**Author:** Ameh Ochai Moses  
**Field:** Biochemistry / Computational Biology & Bioinformatics
