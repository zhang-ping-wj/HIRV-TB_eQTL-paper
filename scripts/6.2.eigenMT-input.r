#!/usr/bin/env Rscript
# eigenMT-input.r - prepare eigenMT input
# usage: Rscript ~/eigenMT-input.r --CHR chr${SLURM_ARRAY_TASK_ID}  --output QTL.chr${SLURM_ARRAY_TASK_ID}.txt --verbose


suppressPackageStartupMessages({
  library(argparse)
  library(data.table)
})

# ---- Argument parsing ----
parser <- ArgumentParser(description = "Prepare QTL file for eigenMT")
parser$add_argument("--CHR", required = TRUE, help = "Chromosome identifier (e.g. chr1 or 1)")
parser$add_argument("--output", required = TRUE, help = "Output QTL file path")
args <- parser$parse_args()

# ---- Input / Output paths ----
chr <- gsub("^chr", "", args$CHR)
input_path <- sprintf("input/LMM-output-chr%s_PC47.txt.gz", chr)
output_path <- args$output

# ---- Read and prepare data ----
QTL <- fread(input_path, stringsAsFactors = FALSE)
cols <- c("snp", "gene_name", "Estimate", "t.value", "p")

if (!all(cols %in% colnames(QTL))) {
  stop("Input file missing required columns: ", paste(setdiff(cols, colnames(QTL)), collapse = ", "))
}

QTL <- QTL[, ..cols]
setnames(QTL, old = cols, new = c("SNP", "gene", "beta", "t-stat", "p-value"))

# ---- Write output ----
fwrite(QTL, output_path, sep = "\t", quote = FALSE, na = "NA")

# ---- Done ----
cat("Saved eigenMT QTL file:", output_path, "\n")




