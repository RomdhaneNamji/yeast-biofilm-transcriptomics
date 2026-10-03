# RNA-seq Analysis of Yeast Colony Spatiotemporal Patterns

## Overview
This project reproduces an RNA-seq study investigating how gene expression in yeast colonies is influenced by genotype and spatial/temporal context.

## Dataset
- Source: SRA (36 samples)
- Conditions:
  - Day 2
  - Day 5 outside
  - Day 5 inside
- Genotypes:
  - WT (F13)
  - tec1Δ
  - sfl1Δ
  - dig1Δ

## Pipeline

1. Data download (SRA Toolkit)
2. Quality control (FastQC, MultiQC)
3. Alignment (Bowtie2)
4. Read counting (featureCounts)
5. Differential analysis (edgeR)
6. Statistical modeling (two-factor ANOVA)
7. Visualization (heatmaps, cluster profiles)

## Key Findings

- Spatiotemporal effects dominate gene expression
- Genotype effects are limited and specific
- Most genes follow additive regulation
- Few genes show interaction effects

## Structure
scripts/ # bash + R scripts
results/ # output tables
plots/ # figures
R_project/ # analysis code
