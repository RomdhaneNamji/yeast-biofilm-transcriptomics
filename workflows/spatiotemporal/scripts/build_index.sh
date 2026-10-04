#!/bin/bash
set -euo pipefail


BASE="/mnt/d/GIS/RNAseq_yeast_2024"
REF="$BASE/reference"
OUT="$REF/index"
LOG="$BASE/logs/build_index.log"

FASTA="$REF/extended_reference.fa"

mkdir -p  "$OUT" "$(dirname "$LOG")"

echo "=== Building Bowtie2 index =echo "=== BUILD INDEX START ===" | tee -a "$LOG"
date | tee -a "$LOG"
echo "Reference FASTA: $FASTA" | tee -a "$LOG"
echo "Output prefix: $OUT/yeast_ext" | tee -a "$LOG"

if [[ ! -f "$FASTA" ]]; then
    echo "ERROR: FASTA file not found: $FASTA" | tee -a "$LOG"
    exit 1
fi

echo "Running bowtie2-build..." | tee -a "$LOG"
bowtie2-build "$FASTA" "$OUT/yeast_ext" >> "$LOG" 2>&1

echo "=== BUILD INDEX DONE ===" | tee -a "$LOG"
date | tee -a "$LOG"

echo "Index files created in: $OUT"=="

bowtie2-build \
	"$REF/extended_reference.fa" \
	"$OUT/yeast_ext"

echo "Index built at: $OUT"


