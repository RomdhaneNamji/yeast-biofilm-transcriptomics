#!/bin/bash

INPUT="./metadata/srr_runs.txt"
OUTDIR="./raw_fastq"
LOGDIR="./logs"

mkdir -p "$OUTDIR" "$LOGDIR"

while read -r srr; do
  echo "Processing $srr"

  # Skip if already downloaded
  if [[ -f "$OUTDIR/${srr}_1.fastq.gz" ]]; then
    echo "$srr already exists, skipping"
    continue
  fi

  # Download
  prefetch "$srr" >> "$LOGDIR/prefetch.log" 2>&1

  # Convert to FASTQ
  fasterq-dump "$srr" --split-files --threads 6 -O "$OUTDIR" >> "$LOGDIR/fasterq.log" 2>&1

  # Compress
  gzip "$OUTDIR/${srr}_1.fastq"
  gzip "$OUTDIR/${srr}_2.fastq"

  echo "$srr DONE"
done < "$INPUT"
