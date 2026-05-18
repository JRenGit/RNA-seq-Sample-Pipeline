# Script to run DESeq2
# Accesses paths and parameters via the snakemake object

library(DESeq2)
library(ggplot2)
library(dplyr)
library(readr)
library(EnhancedVolcano)
library(pheatmap)

# Reroute R warnings to the Snakemake log
log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "message")

# Load featureCounts output
count_data <- read.delim(snakemake@input[["counts"]], comment.char = "#", row.names = 1)

# Strip featureCounts metadata column and extract sample name
count_data <- count_data[, -(1:5)]
colnames(count_data) <- gsub(".*\\/|\\.Aligned.*", "", colnames(count_data))

# Load sample metadata and align with count matrix
sample_info <- read.delim(snakemake@input[["samples"]])
sample_info <- sample_info[match(colnames(count_data), sample_info$sample), ]
if (!all(colnames(count_data) == sample_info$sample)) {
  stop("Sample names in count matrix and sample sheet do not match.")
}

# Set reference level for contrast
contrast <- snakemake@params[["contrast"]]
sample_info[[contrast[1]]] <- factor(
  sample_info[[contrast[1]]],
  levels = c(contrast[3], contrast[2])
)

# Build DESeq DataSet
design_formula <- as.formula(snakemake@params[["design"]])

dds <- DESeqDataSetFromMatrix(
  countData = count_data,
  colData = sample_info,
  design = design_formula
)

# Prefilter low-count genes
keep <- rowSums(counts(dds)) >= snakemake@params[["min_count"]]
dds <- dds[keep, ]

message(paste("Genes after filtering:", sum(keep)))                        

# Run DESeq2
dds <- DESeq(dds)
message("Model coefficients available:")
message(paste(resultsNames(dds), collapse = ", "))

res <- results(dds,
    contrast = contrast,
    alpha = snakemake@params[["alpha"]]
)

message("Results Summary:")
message(capture.output(summary(res)))

# log2FC shrinkage using apeglm - match coefficient name to contrast
coef_name <- grep(
    paste0(contrast[1], "_", contrast[2], "_vs_", contrast[3]),
    resultsNames(dds),
    value = TRUE
)

if (length(coef_name) == 1) {
    res_shrunk <- lfcShrink(dds, coef = coef_name, type = "apeglm")
    message("LFC Shrinkage applied (apeglm).")
} else {
    res_shrunk <- lfcShrink(dds, contrast = contrast, type = "ashr")
    message("LFC shrinkage applied (ashr fallback - apeglm coef not found).")
}

# Generate output tables
res_df <- as.data.frame(res_shrunk) %>%
    mutate(gene_id = rownames(.)) %>%
    filter(!is.na(padj)) %>%
    arrange(padj) %>%
    relocate(gene_id)

write_csv(res_df, snakemake@output[["results_table"]])

norm_counts <- as.data.frame(counts(dds, normalized = TRUE)) %>%
  mutate(gene_id = rownames(.)) %>%
  relocate(gene_id)

write_csv(norm_counts, snakemake@output[["normalized_counts"]])

# Generate QC plots

# PCA to ensure samples cluster by condition, not batch
vst_counts <- vst(dds, blind = FALSE)

pca_plot <- plotPCA(vst_counts, intgroup = contrast[1]) + 
    theme_minimal() +
    ggtitle("PCA: Variance-Stabilized Counts")
ggsave(snakemake@output[["pca_plot"]], pca_plot, width = 8, height = 6)

# Generate volcano plot - statistical significance and biologically meaningful
volcano_plot <- EnhancedVolcano(
    res_df,
    lab = rownames(res_df),
    x = "log2FoldChange",
    y = "padj",
    pCutoff = snakemake@params[["alpha"]],
    FCcutoff = 1,
    title = paste(contrast[2], "vs", contrast[3])
)
ggsave(snakemake@output[["volcano_plot"]], volcano_plot, width = 10, height = 8)

# MA Plot - relate fold change to expression level
pdf(snakemake@output[["ma_plot"]], width = 8, height = 6)
plotMA(res_shrunk, ylim = c(-4,4), main = "MA Plot (shrunken LFC)")
dev.off()

message("DESeq2 analysis complete")
sink(type = "message")
