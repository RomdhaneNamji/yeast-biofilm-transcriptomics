#!/bin/bash

#Base project directory
BASE = "/mnt/d/GIS/RNAseq_yeast_2024"


RAW="$BASE/raw_fastq" 
QC="$BASE/qc"
FASTQC_OUT="$QC/fastqc_raw"
LOG="$BASE/logs/qc.logs"

mkdir -p "$FASTQC_OUT" 


echo "==== QC start ==== " | tee -a "$LOG"
date | tee -a "$LOG"


#Run FastQC

echo " === Running fastQC ===" | tee -a "$LOG"
fastqc "$RAW/*.fastq.gz" -o "$FASTQC_OUT" --threads 10 >> "$LOG" 2>$1

#Run MultiQC

echo " === Running MultiQC ===" | tee -a "$LOG"
multiqc "$FASTQC_OUT" -o "$QC" >> "$LOG" 2>$1

echo " === QC DONE ===" | tee -a "$LOG"
date | tee -a "$LOG"
