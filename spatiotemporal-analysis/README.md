# Yeast Colony Transcriptomics: Spatiotemporal & Genotypic RNA-seq Analysis

[![R](https://img.shields.io/badge/R-276DC3?style=flat&logo=r&logoColor=white)]()
[![Bash](https://img.shields.io/badge/Bash-4EAA25?style=flat&logo=gnu-bash&logoColor=white)]()
[![RNA-seq](https://img.shields.io/badge/Analysis-RNA--seq-blue)]()

## Overview
This repository contains a comprehensive RNA-seq analysis pipeline investigating how gene expression in yeast colonies is influenced by intersecting genotypic, spatial, and temporal variables. 

**Primary Reference:**
Cromie, G. A., Tan, Z., Hays, M., Sirr, A., Dudley, A. M. (2024). *Spatiotemporal patterns of gene expression during development of yeast colonies*. PLOS One. [DOI: 10.1371/journal.pone.0311061](https://doi.org/10.1371/journal.pone.0311061)

## Dataset Specifications
The workflow processes 36 raw RNA-seq samples sourced from the NCBI Sequence Read Archive (SRA), structured across two primary dimensions:
* **Developmental & Spatial Conditions:** Early developmental stage (Day 2), Mature Outside layer (Day 5), and Mature Inside layer (Day 5).
* **Genotypes:** Wild-Type (F13) alongside three transcription factor knockouts (`tec1Δ`, `sfl1Δ`, `dig1Δ`).

## Pipeline Architecture
The analysis is executed using a standard Bash/R command-line toolchain, taking raw data through read alignment, statistical modeling, and final visualization.

1. **Data Acquisition:** Retrieval of the 36 raw FASTQ files using the `SRA Toolkit`.
2. **Quality Control:** Comprehensive read evaluation using `FastQC` and `MultiQC`.
3. **Alignment:** Mapping sequencing reads against the yeast reference genome utilizing `Bowtie2`.
4. **Quantification:** Feature counting and count matrix generation via `featureCounts`.
5. **Statistical Modeling:** Two-factor ANOVA and differential expression analysis conducted in R using `edgeR`.
6. **Data Visualization:** Generation of expression heatmaps and cluster profiles in R to map localized spatiotemporal gradients.

## Key Biological Findings
The statistical analysis revealed clear hierarchies regarding the regulatory drivers of biofilm development:
* **Spatiotemporal Dominance:** Physical location within the colony (Inside vs. Outside) and chronological age (Day 2 vs. Day 5) overwhelmingly dominate gene expression variance across the biofilm.
* **Localized Genotypic Effects:** The genotypic effects of the specific knockouts (`tec1Δ`, `sfl1Δ`, `dig1Δ`) are highly specific and limited in scope compared to the massive shifts driven by spatial positioning.
* **Additive Regulation:** Most genes follow an additive regulatory pattern. Very few genes exhibit complex interaction effects between genotype and spatiotemporal context, indicating that these knockouts generally shift baseline expression uniformly rather than causing location-specific disruptions.

## Repository Structure
* `scripts/` — Bash processing workflows and R statistical scripts.
* `results/` — Output tables, count matrices, and differential expression logs.
* `plots/` — Generated heatmaps, cluster profiles, and visualizations.
* `R_project/` — The structured R environment and analysis code.

---
**Romdhane MRAD NAMJI**  
MSc Bioinformatics Candidate | Pázmány Péter Catholic University  
