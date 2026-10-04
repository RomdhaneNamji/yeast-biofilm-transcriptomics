# Yeast Colony Biofilm Transcriptomics: Spatiotemporal Dynamics & Metabolic Differentiation

[![R](https://img.shields.io/badge/R-276DC3?style=flat&logo=r&logoColor=white)]()
[![Bash](https://img.shields.io/badge/Bash-4EAA25?style=flat&logo=gnu-bash&logoColor=white)]()
[![RNA-seq](https://img.shields.io/badge/Analysis-RNA--seq-blue)]()

## Overview
This repository consolidates two bioinformatics workflows for analyzing RNA-seq data derived from yeast colony biofilms. The projects investigate how gene expression in yeast colonies is driven by spatial localization, temporal development, and genotype, highlighting the complex metabolic adaptations that occur during biofilm maturation.

## Part I: Metabolic Differentiation (2017 Study Reproduction)
An independent computational reproduction of a major yeast biofilm transcriptomics study, validating the distinct metabolic signatures of localized cellular sub-populations.

**Primary Reference:** 
Maršíková, J., Wilkinson, D., Hlaváček, O., et al. (2017). *Metabolic differentiation of surface and invasive cells of yeast colony biofilms revealed by gene expression profiling*. BMC Genomics 18, 814. [DOI: 10.1186/s12864-017-4214-4](https://doi.org/10.1186/s12864-017-4214-4)

### Biological Objectives & Findings
* **Goal:** Reproduce the main transcriptional differences between aerial (surface) and root (invasive) cells.
* **Findings:** 
  * Successfully validated the original study's biological patterns.
  * **Aerial Cells:** Exhibited significant upregulation in signatures related to stress response, glucose starvation, and sporulation.
  * **Root Cells:** Demonstrated strong signatures related to active translation and nutrient transport.

### Pipeline Implementation
* **Primary Processing:** Quality control via `FastQC` and `MultiQC`.
* **Alignment & Quantification:** Read mapping to the reference genome using `HISAT2`, manipulated with `samtools`, and counted via `featureCounts`.
* **Statistical Analysis:** Differential expression and Gene Ontology (GO) enrichment analysis executed in R utilizing `DESeq2`.

---

## Part II: Spatiotemporal & Genotypic Expression Patterns
A comprehensive analysis of 36 RNA-seq samples investigating the interplay between spatial positioning, temporal aging, and genetic knockouts in yeast colonies.

**Primary Reference:**
Cromie, G. A., Tan, Z., Hays, M., Sirr, A., Dudley, A. M. (2024). *Spatiotemporal patterns of gene expression during development of yeast colonies*. PLOS One. [DOI: 10.1371/journal.pone.0311061](https://doi.org/10.1371/journal.pone.0311061)

### Biological Objectives & Findings
* **Dataset Context:** Analyzed across temporal/spatial variables (Day 2 vs. Day 5 Outside vs. Day 5 Inside) and genotypic variables (Wild-Type F13, tec1Δ, sfl1Δ, dig1Δ).
* **Goal:** Model how spatiotemporal context and specific genotypes influence the global transcriptional landscape.
* **Findings:**
  * Spatiotemporal effects heavily dominate gene expression variance across the biofilm.
  * Genotype effects are highly specific and limited in scope compared to spatial positioning.
  * Most genes follow an additive regulatory pattern, with very few exhibiting complex interaction effects between time, space, and genotype.

### Pipeline Implementation
* **Data Acquisition & QC:** Raw data retrieval via `SRA Toolkit`, followed by `FastQC` and `MultiQC`.
* **Alignment & Quantification:** Read alignment utilizing `Bowtie2`, followed by expression quantification with `featureCounts`.
* **Statistical Modeling:** Two-factor ANOVA and differential expression analysis conducted in R using `edgeR`.
* **Visualization:** Generation of expression heatmaps and cluster profiles to map spatiotemporal gradients.


### Reproducibility & Setup
This project utilizes a Conda environment to ensure seamless reproducibility across local workstations and HPC clusters.

**1. Clone the repository:**
\`\`\`bash
git clone https://github.com/RomdhaneNamji/yeast-biofilm-transcriptomics.git
cd yeast-biofilm-transcriptomics
\`\`\`

**2. Create and activate the environment:**
\`\`\`bash
conda env create -f environment.yml
conda activate yeast-transcriptomics
\`\`\`
*(Note: This environment automatically installs R, necessary Bioconductor packages, and all required Bash utilities including HISAT2, Bowtie2, and samtools).*



---
## Repository Structure
```text
📦 yeast-biofilm-transcriptomics
 ┣ 📂 workflows
 ┃ ┣ 📂 reproduction-2017       # Complete pipeline (Bash/R), results, and plots for the 2017 study
 ┃ ┗ 📂 spatiotemporal          # Complete pipeline (Bash/R), results, and plots for the 2024 study
 ┣ 📜 .gitignore
 ┣ 📜 LICENSE
 ┗ 📜 README.md                 # Project documentation
```

**Romdhane MRAD NAMJI**
MSc Bioinformatics Candidate | Pázmány Péter Catholic University
