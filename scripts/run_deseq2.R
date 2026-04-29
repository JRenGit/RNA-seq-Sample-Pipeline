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
sample_info <- read.delim(snakemake@input[["samples"]]
sample_info <- sample_info[math(colnames(count_data), sample_info$sample), ]
if (!all(colnames(count_data) == sample_info$sample)) {
  stop("Sample names in count matrix and sample sheet do not match.")
}

                          
