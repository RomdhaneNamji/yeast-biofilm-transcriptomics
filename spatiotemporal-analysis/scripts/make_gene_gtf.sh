#!/bin/bash
set -euo pipefail

IN="/mnt/d/GIS/RNAseq_yeast_2024/reference/extended_annotation.gff"
OUT="/mnt/d/GIS/RNAseq_yeast_2024/reference/extended_annotation.gtf"

awk -F '\t' '
BEGIN { OFS="\t" }
$0 ~ /^#/ { next }
NF < 9 { next }
$3 != "gene" { next }

{
    attr = $9
    gene_id = ""

    if (match(attr, /ID[ ="]+([^";]+)/, m)) {
        gene_id = m[1]
    } else if (match(attr, /Name=([^;]+)/, m)) {
        gene_id = m[1]
    } else if (match(attr, /Name[ ="]+([^";]+)/, m)) {
        gene_id = m[1]
    }

    if (gene_id != "") {
        print $1, $2, $3, $4, $5, $6, $7, $8, "gene_id \"" gene_id "\"; gene_name \"" gene_id "\";"
    }
}
' "$IN" > "$OUT"

echo "Wrote: $OUT"
echo "Preview:"
head "$OUT"
