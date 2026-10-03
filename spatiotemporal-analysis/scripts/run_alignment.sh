#!/bin/bash

set -euo pipefail
BASE="/mnt/d/GIS/RNAseq_yeast_2024"
RAW_EXT="$BASE/raw_fastq"
IDX="$BASE/reference/index/yeast_ext"
BAM_EXT="$BASE/alignments"
LOG_EXT="$BASE/logs"

TMP_BASE="$HOME/yeast_tmp"
TMP_RAW="$TMP_BASE/raw"
TMP_WORK="$TMP_BASE/work"

BOWTIE2_THREADS=8
SAMTOOLS_THREADS=4

mkdir -p "$TMP_RAW" "$TMP_WORK" "$BAM_EXT" "$LOG_EXT"

find "$RAW_EXT" -name "*_1.fastq.gz" | sort > "$LOG_EXT/samples.txt"

while read -r R1_EXT; do

    SAMPLE=$(basename "$R1_EXT" _1.fastq.gz)
    R2_EXT="$RAW_EXT/${SAMPLE}_2.fastq.gz"

    LOG="$LOG_EXT/${SAMPLE}.align.log"
    FINAL_BAM="$BAM_EXT/${SAMPLE}.sorted.bam"

    if [[ -f "$FINAL_BAM" ]]; then
        echo "[SKIP] $SAMPLE"
        continue
    fi

    echo "=== START $SAMPLE ===" | tee "$LOG"

    TMP_R1="$TMP_RAW/${SAMPLE}_1.fastq.gz"
    TMP_R2="$TMP_RAW/${SAMPLE}_2.fastq.gz"

    TMP_BAM="$TMP_WORK/${SAMPLE}.sorted.bam"

    # Copy FASTQ locally
    cp "$R1_EXT" "$TMP_R1"
    cp "$R2_EXT" "$TMP_R2"

    # Alignment + sorting (SAFE)
    bowtie2 \
      -x "$IDX" \
      -1 "$TMP_R1" \
      -2 "$TMP_R2" \
      -N 1 \
      -I 50 \
      -X 450 \
      --reorder \
      -p "$BOWTIE2_THREADS" \
      2>> "$LOG" \
    | samtools view -@ "$SAMTOOLS_THREADS" -b - \
    | samtools sort -@ "$SAMTOOLS_THREADS" -o "$TMP_BAM" 2>> "$LOG"

    samtools index "$TMP_BAM" 2>> "$LOG"

    # Move result to external drive
    mv "$TMP_BAM" "$FINAL_BAM"
    mv "${TMP_BAM}.bai" "$FINAL_BAM.bai"

    # Cleanup
    rm -f "$TMP_R1" "$TMP_R2"

    echo "=== DONE $SAMPLE ===" | tee -a "$LOG"

done < "$LOG_EXT/samples.txt"
