# Yeast Colony Biofilm Transcriptomics: Metabolic Differentiation Reproduction

[![R](https://img.shields.io/badge/R-276DC3?style=flat&logo=r&logoColor=white)]()
[![Bash](https://img.shields.io/badge/Bash-4EAA25?style=flat&logo=gnu-bash&logoColor=white)]()
[![RNA-seq](https://img.shields.io/badge/Analysis-RNA--seq-blue)]()

## Overview
This repository contains a computational reproduction workflow for a foundational yeast colony biofilm RNA-seq study. The project validates the localized transcriptomic signatures and metabolic reprogramming that distinguish the aerial (surface) layer from the root (invasive) layer of mature biofilms.

**Primary Reference:**
Maršíková, J., Wilkinson, D., Hlaváček, O., et al. (2017). *Metabolic differentiation of surface and invasive cells of yeast colony biofilms revealed by gene expression profiling*. BMC Genomics 18, 814. [DOI: 10.1186/s12864-017-4214-4](https://doi.org/10.1186/s12864-017-4214-4)

## Project Goals
- **Independent Validation:** Reproduce the primary transcriptional differences between aerial and root cellular sub-populations.
- **End-to-End Pipeline Execution:** Perform comprehensive quality control, read alignment, transcript counting, differential expression modeling, and Gene Ontology (GO) analysis.
- **Methodological Benchmarking:** Compare the independently reproduced computational results with the published study to confirm biological consistency.

## Computational Workflow
The analysis is executed using a standard command-line and R toolchain, taking raw data through to functional enrichment.

1. **Data Acquisition:** Download raw RNA-seq datasets.
2. **Quality Control:** Run `FastQC` and aggregate reports with `MultiQC`.
3. **Alignment:** Map sequencing reads to the reference genome utilizing `HISAT2`.
4. **Quantification:** Count mapped reads against genomic features using `featureCounts`.
5. **Differential Expression:** Perform statistical modeling and expression analysis using `DESeq2`.
6. **Functional Analysis:** Run GO enrichment analysis to interpret biological pathway activation.

## Main Findings
The reproduced pipeline successfully validated the core biological patterns reported in the published study, revealing a stark metabolic dichotomy driven by physical localization within the colony:

* **Aerial Cells (Surface):** Exhibited significant upregulation in signatures related to generalized stress responses, glucose starvation, and sporulation. 
* **Root Cells (Invasive):** Demonstrated strong transcriptional signatures related to active translation and transmembrane nutrient transport.
* **Conclusion:** The independently reproduced biological pattern was entirely consistent with the original published findings.

## Main Tools
* **Command Line Utilities:** `HISAT2`, `samtools`, `featureCounts` (Subread)
* **Statistical Environment:** `R`, `DESeq2`

## Notes
Raw `FASTQ` and aligned `BAM` files are excluded from this repository due to file size constraints. 

---
**Romdhane MRAD NAMJI**  
MSc Bioinformatics Candidate | Pázmány Péter Catholic University  
