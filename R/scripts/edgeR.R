## Reproduction of Cromie et al.paper

setwd("D:/GIS/RNAseq_yeast_2024/R/")

library (edgeR)


# ---------------- Paths -----------

counts_file <- "data/gene_counts.txt" 
metadata_file <- "data/sample_metadata.csv"
results_dir <- "results"
plots_dir <- "plots"

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

#------------------Read featureCounts table ---

counts <- read.delim(
  counts_file,
  comment.char = "#",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# first 6 columns are feature annotation from featureCounts
count_matrix <- counts[,7:ncol(counts), drop = FALSE]
rownames(count_matrix) <- counts$Geneid

# clean sample names
sample_names <- colnames(count_matrix)
sample_names <- basename(sample_names)
sample_names <- sub("\\.sorted\\.bam$", "", sample_names)
colnames(count_matrix) <- sample_names

#convert to integer matrix
count_matrix <- as.matrix(count_matrix)
storage.mode <- "integer" 


# ---------- Read metadata -------------

meta <- read.csv(metadata_file, stringsAsFactors = FALSE)
required_cols <- c("sample", "genotype", "spatiotemporal")
missing_cols <- setdiff(required_cols,colnames(meta))

if (length(missing_cols) > 0) {
  stop("Missing required metadata columns: ", paste(missing_cols, collapse = ", "))
}

# Reorder metadata to match counts 
meta <- meta[match(colnames(count_matrix),meta$sample), ]
if (any(is.na(meta$sample))) {
  stop("Some count columns are missing from sample_metadata.csv")
}
if (!all(meta$sample == colnames(count_matrix))) {
  stop("Metadata order mismatch after matching")
}

# Set factor levels
meta$genotype <- factor(meta$genotype, levels = c("tec1","wt","sfl1","dig1"))
meta$condition <- factor(
  meta$condition, levels = c("day2","day5_inside","day5_outside")
)

if (any(is.na(meta$genotype))) {
  stop("Invalid genotype values found in metadata")
}
if (any(is.na(meta$condition))) {
  stop("Invalid spatiotemporal values found in metadata")
}

write.csv(meta, file.path(results_dir, "metadata_reordered.csv"), row.names = FALSE)


# --------------- Create DGElist and normaize -------
dge <- DGEList(counts = count_matrix)
dge <- calcNormFactors(dge, method = "TMM")

# -------------- Filter S288c ORFs beginning with Y
keep_orf <- grepl("^Y", rownames(dge))
dge <- dge[keep_orf, keep.lib.sizes = FALSE]


#--------------- Convert to log2 CPM prior count 20--------

logCPM <- cpm(dge, log = TRUE, prior.count = 20)
write.csv(logCPM, file.path(results_dir,"logCPM_matrix.csv"))


# ----------------------------
# 8) Per-gene ANOVA models
# Paper used:
# - additive model: genotype + spatiotemporal
# - interaction model: genotype * spatiotemporal
# and Holm correction on relevant p-values
# ----------------------------

n_genes <- nrow(logCPM)

add_p_genotype         <- numeric(n_genes)
add_p_spatiotemporal   <- numeric(n_genes)
add_p_vs_null          <- numeric(n_genes)
add_rsq                <- numeric(n_genes)

int_p_interaction      <- numeric(n_genes)
int_p_vs_null          <- numeric(n_genes)
int_rsq                <- numeric(n_genes)

gene_ids <- rownames(logCPM)

for (i in seq_len(n_genes)) {
  y <- as.numeric(logCPM[i, ])
  
  df <- data.frame(
    y = y,
    genotype = meta$genotype,
    spatiotemporal = meta$condition
  )
  
  fit_null <- lm(y ~ 1, data = df)
  fit_add  <- lm(y ~ genotype + spatiotemporal, data = df)
  fit_int  <- lm(y ~ genotype * spatiotemporal, data = df)
  
  # Additive model term p-values
  d_add <- drop1(fit_add, test = "F")
  add_p_genotype[i]       <- d_add["genotype", "Pr(>F)"]
  add_p_spatiotemporal[i] <- d_add["spatiotemporal", "Pr(>F)"]
  
  # Model vs null
  add_p_vs_null[i] <- anova(fit_null, fit_add)$`Pr(>F)`[2]
  int_p_vs_null[i] <- anova(fit_null, fit_int)$`Pr(>F)`[2]
  
  # Interaction term p-value
  d_int <- drop1(fit_int, test = "F")
  int_p_interaction[i] <- d_int["genotype:spatiotemporal", "Pr(>F)"]
  
  # R-squared
  add_rsq[i] <- summary(fit_add)$r.squared
  int_rsq[i] <- summary(fit_int)$r.squared
}

# ----------------------------
# 9) Multiple testing correction (Holm)
# Done separately for each family, matching the paper description
# ----------------------------
res <- data.frame(
  Geneid = gene_ids,
  variance_logCPM = apply(logCPM, 1, var),
  
  add_p_genotype = add_p_genotype,
  add_p_spatiotemporal = add_p_spatiotemporal,
  add_p_vs_null = add_p_vs_null,
  add_rsq = add_rsq,
  
  int_p_interaction = int_p_interaction,
  int_p_vs_null = int_p_vs_null,
  int_rsq = int_rsq,
  
  stringsAsFactors = FALSE
)

res$add_padj_genotype       <- p.adjust(res$add_p_genotype, method = "holm")
res$add_padj_spatiotemporal <- p.adjust(res$add_p_spatiotemporal, method = "holm")
res$add_padj_vs_null        <- p.adjust(res$add_p_vs_null, method = "holm")

res$int_padj_interaction    <- p.adjust(res$int_p_interaction, method = "holm")
res$int_padj_vs_null        <- p.adjust(res$int_p_vs_null, method = "holm")

# ----------------------------
# 9) Multiple testing correction (Holm)
# Done separately for each family, matching the paper description
# ----------------------------
res <- data.frame(
  Geneid = gene_ids,
  variance_logCPM = apply(logCPM, 1, var),

  add_p_genotype = add_p_genotype,
  add_p_spatiotemporal = add_p_spatiotemporal,
  add_p_vs_null = add_p_vs_null,
  add_rsq = add_rsq,

  int_p_interaction = int_p_interaction,
  int_p_vs_null = int_p_vs_null,
  int_rsq = int_rsq,

  stringsAsFactors = FALSE
)

res$add_padj_genotype       <- p.adjust(res$add_p_genotype, method = "holm")
res$add_padj_spatiotemporal <- p.adjust(res$add_p_spatiotemporal, method = "holm")
res$add_padj_vs_null        <- p.adjust(res$add_p_vs_null, method = "holm")

res$int_padj_interaction    <- p.adjust(res$int_p_interaction, method = "holm")
res$int_padj_vs_null        <- p.adjust(res$int_p_vs_null, method = "holm")

# ----------------------------
# 10) Paper-style subsets
# ----------------------------

# Genes significantly responding to additive model factors
res$signif_genotype_additive <- res$add_padj_genotype < 0.01
res$signif_spatiotemporal_additive <- res$add_padj_spatiotemporal < 0.01

# Strong additive responses
res$strong_additive <- (
  res$variance_logCPM > 0.2 &
    res$add_rsq > 0.7 &
    res$add_padj_vs_null < 0.01
)

# Strong non-additive responses
res$strong_interaction <- (
  res$variance_logCPM > 0.2 &
    res$int_rsq > 0.9 &
    res$int_padj_interaction <= 0.01
)

# Shared responses
res$responds_to_both_additive <- (
  res$signif_genotype_additive &
    res$signif_spatiotemporal_additive
)

# ----------------------------
# 11) Save main results
# ----------------------------
write.csv(res, file.path(results_dir, "ANOVA_results_all_genes.csv"), row.names = FALSE)
write.csv(
  subset(res, signif_genotype_additive),
  file.path(results_dir, "genes_significant_genotype_additive.csv"),
  row.names = FALSE
)
write.csv(
  subset(res, signif_spatiotemporal_additive),
  file.path(results_dir, "genes_significant_spatiotemporal_additive.csv"),
  row.names = FALSE
)
write.csv(
  subset(res, responds_to_both_additive),
  file.path(results_dir, "genes_significant_both_additive.csv"),
  row.names = FALSE
)
write.csv(
  subset(res, strong_additive),
  file.path(results_dir, "genes_strong_additive.csv"),
  row.names = FALSE
)
write.csv(
  subset(res, strong_interaction),
  file.path(results_dir, "genes_strong_interaction.csv"),
  row.names = FALSE
)

# ----------------------------
# 12) Summary file
# ----------------------------
summary_lines <- c(
  paste("Genes after ORF filter (^Y):", nrow(logCPM)),
  paste("Significant genotype additive (Holm < 0.01):", sum(res$signif_genotype_additive)),
  paste("Significant spatiotemporal additive (Holm < 0.01):", sum(res$signif_spatiotemporal_additive)),
  paste("Significant for both additive factors:", sum(res$responds_to_both_additive)),
  paste("Strong additive genes:", sum(res$strong_additive)),
  paste("Strong interaction genes:", sum(res$strong_interaction))
)

writeLines(summary_lines, con = file.path(results_dir, "analysis_summary.txt"))

cat(paste(summary_lines, collapse = "\n"), "\n")