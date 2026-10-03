# =========================================
# Paper-style cluster plots
# =========================================

setwd("D:/GIS/RNAseq_yeast_2024/R/")

library(edgeR)
library(pheatmap)


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

sample_names <- basename(colnames(count_matrix))
sample_names <- sub("\\.sorted\\.bam$", "", sample_names)
colnames(count_matrix) <- sample_names

count_matrix <- as.matrix(count_matrix)
storage.mode(count_matrix) <- "integer"

# ----------------------------
# 3. Read metadata
# ----------------------------
meta <- read.csv(metadata_file, stringsAsFactors = FALSE)

if ("condition" %in% colnames(meta) && !("spatiotemporal" %in% colnames(meta))) {
  meta$spatiotemporal <- meta$condition
}

meta$genotype <- factor(meta$genotype, levels = c("tec1", "wt", "sfl1", "dig1"))
meta$spatiotemporal <- factor(
  meta$spatiotemporal,
  levels = c("day2", "day5_outside", "day5_inside")
)

meta <- meta[match(colnames(count_matrix), meta$sample), ]

if (any(is.na(meta$sample))) stop("Metadata mismatch")

# ----------------------------
# 4. Normalize counts
# ----------------------------
dge <- DGEList(counts = count_matrix)
dge <- calcNormFactors(dge)

keep_orf <- grepl("^Y", rownames(dge))
dge <- dge[keep_orf, , keep.lib.sizes = FALSE]

logCPM <- cpm(dge, log = TRUE, prior.count = 20)

# ----------------------------
# 5. Read ANOVA results
# ----------------------------
anova_res <- read.csv(anova_file, stringsAsFactors = FALSE)
rownames(anova_res) <- anova_res$Geneid
anova_res <- anova_res[rownames(logCPM), , drop = FALSE]

# ----------------------------
# 6. Reorder columns like the paper
# Order:
# day2 -> day5_outside -> day5_inside
# and within each block:
# tec1 -> wt -> sfl1 -> dig1
# ----------------------------
meta$replicate_index <- ave(seq_len(nrow(meta)),
                            meta$genotype, meta$spatiotemporal,
                            FUN = seq_along)

ord <- order(meta$spatiotemporal, meta$replicate_index, meta$genotype)

meta_ord <- meta[ord, ]
logCPM_ord <- logCPM[, ord, drop = FALSE]

# ----------------------------
# 7. Helper functions
# ----------------------------
normalize_for_spatiotemporal <- function(mat, meta) {
  out <- mat
  for (lvl in levels(meta$genotype)) {
    idx <- meta$genotype == lvl
    out[, idx] <- sweep(out[, idx, drop = FALSE], 1, rowMeans(out[, idx, drop = FALSE]), "-")
  }
  out
}

scale_rows <- function(mat) {
  z <- t(scale(t(mat)))
  z[is.na(z)] <- 0
  z
}

geno_letters <- c(tec1 = "T", wt = "F", sfl1 = "S", dig1 = "D")
geno_colors  <- c(tec1 = "#E69F00", wt = "#666666", sfl1 = "#009E73", dig1 = "#56B4E9")

draw_paper_profile <- function(submat, meta, title, outfile) {
  mean_profile <- colMeans(submat)
  
  pdf(outfile, width = 7.5, height = 4.2)
  
  y_lim <- range(submat, na.rm = TRUE)
  
  plot(
    seq_len(ncol(submat)),
    mean_profile,
    type = "n",
    ylim = y_lim,
    xaxt = "n",
    xlab = "",
    ylab = "Normalized expression",
    main = title
  )
  
  # shaded condition blocks
  block_sizes <- table(meta$spatiotemporal)
  starts <- c(1, cumsum(block_sizes)[-length(block_sizes)] + 1)
  ends   <- cumsum(block_sizes)
  
  rect(starts[1] - 0.5, y_lim[1], ends[1] + 0.5, y_lim[2], col = "#f0f0f0", border = NA)
  rect(starts[2] - 0.5, y_lim[1], ends[2] + 0.5, y_lim[2], col = "#dddddd", border = NA)
  rect(starts[3] - 0.5, y_lim[1], ends[3] + 0.5, y_lim[2], col = "#f0f0f0", border = NA)
  
  # subset of genes in gray
  set.seed(1)
  idx <- sample(seq_len(nrow(submat)), min(80, nrow(submat)))
  for (i in idx) {
    lines(seq_len(ncol(submat)), submat[i, ], col = "gray70")
  }
  
  # mean profile
  lines(seq_len(ncol(submat)), mean_profile, col = "red", lwd = 2)
  
  # separators
  abline(v = ends + 0.5, lty = 2, col = "gray40")
  
  # bottom labels
  mids <- c((starts[1] + ends[1]) / 2,
            (starts[2] + ends[2]) / 2,
            (starts[3] + ends[3]) / 2)
  
  axis(1, at = mids, labels = c("Day 2", "Day 5 outside", "Day 5 inside"), tick = FALSE, line = 1)
  
  # top genotype letters
  y_top <- par("usr")[4]
  text(
    x = seq_len(ncol(submat)),
    y = y_top + 0.04 * diff(par("usr")[3:4]),
    labels = geno_letters[as.character(meta$genotype)],
    col = geno_colors[as.character(meta$genotype)],
    xpd = NA,
    cex = 0.9,
    font = 2
  )
  
  box()
  dev.off()
}

# ----------------------------
# 8. Build spatiotemporal-responsive gene set
# ----------------------------
spatio_genes <- rownames(anova_res)[anova_res$add_padj_spatiotemporal < 0.01 & anova_res$add_rsq > 0.7]
spatio_genes <- intersect(spatio_genes, rownames(logCPM_ord))

mat_spatio <- normalize_for_spatiotemporal(logCPM_ord[spatio_genes, , drop = FALSE], meta_ord)
mat_spatio_z <- scale_rows(mat_spatio)

# ----------------------------
# 9. Paper-style heatmap
# ----------------------------
hc <- hclust(dist(mat_spatio_z), method = "complete")
clusters <- cutree(hc, k = 8)

write.csv(
  data.frame(Geneid = names(clusters), cluster = clusters),
  file.path(results_dir, "spatiotemporal_clusters.csv"),
  row.names = FALSE
)

pdf(file.path(plots_dir, "spatiotemporal_heatmap.pdf"), width = 10, height = 8)

# order rows by cluster, then dendrogram order
row_order <- order(clusters)

heatmap_mat <- mat_spatio_z[row_order, , drop = FALSE]

pheatmap(
  heatmap_mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  main = "Spatiotemporal-responsive genes"
)

dev.off()

# ----------------------------
# 10. Paper-style profile plots for each cluster
# ----------------------------
for (cl in sort(unique(clusters))) {
  genes <- names(clusters)[clusters == cl]
  if (length(genes) < 5) next
  
  submat <- mat_spatio_z[genes, , drop = FALSE]
  
  draw_paper_profile(
    submat = submat,
    meta = meta_ord,
    title = paste("Spatiotemporal profiles - cluster", cl),
    outfile = file.path(plots_dir, paste0("spatiotemporal_cluster_", cl, ".pdf"))
  )
}

cat("Paper-style plots generated.\n")