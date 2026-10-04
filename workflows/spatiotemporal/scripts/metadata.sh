grep -v "^#" /mnt/d/GIS/RNAseq_yeast_2024/counts/gene_counts.txt \
| head -1 \
| awk '{for(i=7;i<=NF;i++) print $i}' \
> /mnt/d/GIS/RNAseq_yeast_2024/R/data/sample_names.txt
