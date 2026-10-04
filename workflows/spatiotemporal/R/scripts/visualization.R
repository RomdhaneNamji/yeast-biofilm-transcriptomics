# =========================================
# Visualization script for yeast colony RNA-seq study
# =========================================

setwd("D:/GIS/RNAseq_yeast_2024/R/")

library(edgeR)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

# ----------------------------
# 1. Paths
# ----------------------------
counts_file   <- "data/gene_counts.txt"
metadata_file <- "data/sample_metadata.csv"
anova_file    <- "results/ANOVA_results_all_genes.csv"

plots_dir     <- "plots"
results_dir   <- "results"

dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# 2. Read counts
# ----------------------------
counts <- read.delim(
  counts_file,
  comment.char = "#",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

count_matrix <- counts[, 7:ncol(counts), drop = FALSE]
rownames(count_matrix) <- counts$Geneid

sample_names <- colnames(count_matrix)
sample_names <- basename(sample_names)
sample_names <- sub("\\.sorted\\.bam$", "", sample_names)
colnames(count_matrix) <- sample_names

count_matrix <- as.matrix(count_matrix)
storage.mode(count_matrix) <- "integer"

# ----------------------------
# 3. Read metadata
# ----------------------------
meta <- read.csv(metadata_file, stringsAsFactors = FALSE)

# support either 'condition' or 'spatiotemporal'
if ("condition" %in% colnames(meta) && !("spatiotemporal" %in% colnames(meta))) {
  meta$spatiotemporal <- meta$condition
}

required_cols <- c("sample", "genotype", "spatiotemporal")
missing_cols <- setdiff(required_cols, colnames(meta))
if (length(missing_cols) > 0) {
  stop("Missing required metadata columns: ", paste(missing_cols, collapse = ", "))
}

meta <- meta[match(colnames(count_matrix), meta$sample), ]

if (any(is.na(meta$sample))) {
  stop("Some count columns are missing from sample_metadata.csv")
}
if (!all(meta$sample == colnames(count_matrix))) {
  stop("Metadata order mismatch after matching")
}

meta$genotype <- factor(meta$genotype, levels = c("tec1", "wt", "sfl1", "dig1"))
meta$spatiotemporal <- factor(
  meta$spatiotemporal,
  levels = c("day2", "day5_outside", "day5_inside")
)

# ----------------------------
# 4. Normalize counts
# ----------------------------
dge <- DGEList(counts = count_matrix)
dge <- calcNormFactors(dge, method = "TMM")

# Keep only ORFs beginning with Y
keep_orf <- grepl("^Y", rownames(dge))
dge <- dge[keep_orf, , keep.lib.sizes = FALSE]

logCPM <- cpm(dge, log = TRUE, prior.count = 20)

# ----------------------------
# 5. Read ANOVA results
# ----------------------------
anova_res <- read.csv(anova_file, stringsAsFactors = FALSE)
rownames(anova_res) <- anova_res$Geneid

# Keep only genes present in logCPM
anova_res <- anova_res[rownames(logCPM), , drop = FALSE]

# ----------------------------
# 6. Plot colors and annotations
# ----------------------------
ann_col <- data.frame(
  Genotype = meta$genotype,
  Spatiotemporal = meta$spatiotemporal
)
rownames(ann_col) <- meta$sample

ann_colors <- list(
  Genotype = c(
    tec1 = "#E69F00",
    wt   = "#999999",
    sfl1 = "#009E73",
    dig1 = "#56B4E9"
  ),
  Spatiotemporal = c(
    day2         = "#000000",
    day5_outside = "#E76F51",
    day5_inside  = "#4CAF50"
  )
)

# ----------------------------
# 7. Refined MDS plot
# ----------------------------
mds <- plotMDS(dge, plot = FALSE)

# var.explained may not exist in some edgeR versions
xlab_text <- "Leading logFC dim 1"
ylab_text <- "Leading logFC dim 2"

if (!is.null(mds$var.explained)) {
  xlab_text <- paste0("Leading logFC dim 1 (", round(mds$var.explained[1], 1), "%)")
  ylab_text <- paste0("Leading logFC dim 2 (", round(mds$var.explained[2], 1), "%)")
}

mds_df <- data.frame(
  Dim1 = mds$x,
  Dim2 = mds$y,
  sample = meta$sample,
  genotype = meta$genotype,
  spatiotemporal = meta$spatiotemporal
)

condition_cols <- c(
  day2 = "black",
  day5_outside = "#E76F51",
  day5_inside = "#4CAF50"
)

genotype_pch <- c(
  tec1 = 15,   # square
  wt   = 16,   # circle
  sfl1 = 17,   # triangle
  dig1 = 18    # diamond
)

pdf(file.path(plots_dir, "01_MDS_plot_refined.pdf"), width = 9, height = 7)

plot(
  mds_df$Dim1,
  mds_df$Dim2,
  type = "n",
  xlab = xlab_text,
  ylab = ylab_text,
  main = "MDS plot of RNA-seq samples"
)

points(
  mds_df$Dim1,
  mds_df$Dim2,
  pch = genotype_pch[as.character(mds_df$genotype)],
  col = condition_cols[as.character(mds_df$spatiotemporal)],
  cex = 1.5
)

legend(
  "topright",
  legend = c("day2", "day5_outside", "day5_inside"),
  col = condition_cols[c("day2", "day5_outside", "day5_inside")],
  pch = 16,
  title = "Spatiotemporal",
  bty = "n"
)

legend(
  "bottomleft",
  legend = c("tec1", "wt", "sfl1", "dig1"),
  pch = genotype_pch[c("tec1", "wt", "sfl1", "dig1")],
  col = "black",
  title = "Genotype",
  bty = "n"
)

dev.off()

# ----------------------------
# 8. Sample distance heatmap
# ----------------------------
sample_cor <- cor(logCPM, method = "pearson")

pdf(file.path(plots_dir, "02_sample_correlation_heatmap.pdf"), width = 10, height = 8)
pheatmap(
  sample_cor,
  annotation_col = ann_col,
  annotation_row = ann_col,
  annotation_colors = ann_colors,
  main = "Sample correlation heatmap"
)
dev.off()

# ----------------------------
# 9. P-value histograms
# ----------------------------
plot_pval_hist <- function(pvals, title, outfile) {
  pvals <- pvals[is.finite(pvals) & !is.na(pvals) & pvals > 0]
  df <- data.frame(log10p = -log10(pvals))
  
  p <- ggplot(df, aes(x = log10p)) +
    geom_histogram(bins = 30, fill = "#d9e7e7", color = "gray50") +
    labs(
      title = title,
      x = expression(-log[10](P-value)),
      y = "Frequency"
    ) +
    theme_bw(base_size = 12)
  
  ggsave(outfile, p, width = 6, height = 5)
}

plot_pval_hist(
  anova_res$add_p_genotype,
  "Genotype factor",
  file.path(plots_dir, "03_pvalue_hist_genotype.pdf")
)

plot_pval_hist(
  anova_res$add_p_spatiotemporal,
  "Spatiotemporal factor",
  file.path(plots_dir, "04_pvalue_hist_spatiotemporal.pdf")
)

plot_pval_hist(
  anova_res$int_p_interaction,
  "Interaction factor",
  file.path(plots_dir, "05_pvalue_hist_interaction.pdf")
)

# ----------------------------
# 10. Helper: normalize expression by factor
# ----------------------------
normalize_for_factor <- function(mat, meta, keep_factor = c("genotype", "spatiotemporal")) {
  keep_factor <- match.arg(keep_factor)
  out <- mat
  
  if (keep_factor == "genotype") {
    # remove spatiotemporal mean effect
    for (lvl in levels(meta$spatiotemporal)) {
      idx <- meta$spatiotemporal == lvl
      out[, idx] <- sweep(out[, idx, drop = FALSE], 1, rowMeans(out[, idx, drop = FALSE]), "-")
    }
  } else {
    # remove genotype mean effect
    for (lvl in levels(meta$genotype)) {
      idx <- meta$genotype == lvl
      out[, idx] <- sweep(out[, idx, drop = FALSE], 1, rowMeans(out[, idx, drop = FALSE]), "-")
    }
  }
  
  out
}

scale_rows <- function(mat) {
  scaled <- t(scale(t(mat)))
  scaled[is.na(scaled)] <- 0
  scaled
}

# ----------------------------
# 11. Heatmap: genotype-responsive genes
# ----------------------------
geno_genes <- rownames(anova_res)[anova_res$add_padj_genotype < 0.01 & anova_res$add_rsq > 0.7]
geno_genes <- intersect(geno_genes, rownames(logCPM))

if (length(geno_genes) > 2) {
  mat_geno <- normalize_for_factor(logCPM[geno_genes, , drop = FALSE], meta, keep_factor = "genotype")
  mat_geno <- scale_rows(mat_geno)
  
  pdf(file.path(plots_dir, "06_heatmap_genotype_responsive.pdf"), width = 12, height = 10)
  pheatmap(
    mat_geno,
    show_rownames = FALSE,
    annotation_col = ann_col,
    annotation_colors = ann_colors,
    clustering_method = "complete",
    main = "Genes responsive to genotype"
  )
  dev.off()
}

# ----------------------------
# 12. Heatmap: spatiotemporal-responsive genes
# ----------------------------
spatio_genes <- rownames(anova_res)[anova_res$add_padj_spatiotemporal < 0.01 & anova_res$add_rsq > 0.7]
spatio_genes <- intersect(spatio_genes, rownames(logCPM))

if (length(spatio_genes) > 2) {
  mat_spatio <- normalize_for_factor(logCPM[spatio_genes, , drop = FALSE], meta, keep_factor = "spatiotemporal")
  mat_spatio <- scale_rows(mat_spatio)
  
  pdf(file.path(plots_dir, "07_heatmap_spatiotemporal_responsive.pdf"), width = 12, height = 10)
  pheatmap(
    mat_spatio,
    show_rownames = FALSE,
    annotation_col = ann_col,
    annotation_colors = ann_colors,
    clustering_method = "complete",
    main = "Genes responsive to spatiotemporal factor"
  )
  dev.off()
}

# ----------------------------
# 13. Heatmap: strong interaction genes
# ----------------------------
int_genes <- rownames(anova_res)[anova_res$strong_interaction]
int_genes <- intersect(int_genes, rownames(logCPM))

if (length(int_genes) > 2) {
  mat_int <- scale_rows(logCPM[int_genes, , drop = FALSE])
  
  pdf(file.path(plots_dir, "08_heatmap_interaction_genes.pdf"), width = 12, height = 10)
  pheatmap(
    mat_int,
    show_rownames = FALSE,
    annotation_col = ann_col,
    annotation_colors = ann_colors,
    clustering_method = "complete",
    main = "Genes with strong genotype × spatiotemporal interaction"
  )
  dev.off()
}


cat("Refined visualization plots generated in:", plots_dir, "\n")