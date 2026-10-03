#!/bin/bash
SECONDS=0

set -euo pipefail

BASE="/mnt/d/GIS/RNAseq_yeast_2024"

BAMDIR="$BASE/alignments"
GFF="$BASE/reference/extended_annotation.gtf"
OUTDIR="$BASE/counts"
LOG="$BASE/logs/featurecounts.log"

mkdir -p "$OUTDIR"

echo "=== FEATURECOUNTS START ===" | tee "$LOG"
date | tee -a "$LOG"

featureCounts \
  -F GTF \
  -T 8 \
  -a "$GFF" \
  -o "$OUTDIR/gene_counts.txt" \
  -t gene \
  -g gene_id \
  -s 2 \
  -p -P -B \
  -d 50 \
  -D 450 \
  "$BAMDIR"/*.sorted.bam \
  >> "$LOG" 2>&1

echo "=== FEATURECOUNTS DONE ===" | tee -a "$LOG"
date | tee -a "$LOG"

duration=$SECONDS
echo "$((duration/60)) minutes and $((minutes%60)) seconds elapsed"
