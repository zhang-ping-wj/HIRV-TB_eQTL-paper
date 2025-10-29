#!/usr/bin/env Rscript
# run_matrix_eqtl.R — run MatrixEQTL using different expression PCs

suppressPackageStartupMessages({
  library(MatrixEQTL)
  library(data.table)
})

# ---- user-editable settings ----
input_dir <- "./1.TB.day2.input"
snps_file <- file.path(input_dir, "SNP_rsid_213.Latent.TB.txt.gz")
genes_file <- file.path(input_dir, "GE.txt")
snpPos_file <- file.path(input_dir, "snpsloc_rsid_213.Latent.TB.txt.gz")
genePos_file <- file.path(input_dir, "geneloc.txt")
covariates_file <- file.path(input_dir, "Covariates.all.PCs.txt")
out_summary <- "zz_day2_latent.csv"
max_PCs <- 100

# ---- load data once ----
snps.X <- SlicedData$new(); snps.X$LoadFile(snps_file)
genes  <- SlicedData$new(); genes$LoadFile(genes_file)

snpPos  <- fread(snpPos_file, stringsAsFactors = FALSE)
genePos <- fread(genePos_file, stringsAsFactors = FALSE)

cov_all <- fread(covariates_file, header = TRUE, row.names = 1)
# ensure covariates are in rows as original script expected
# cov_all rows = PCs, cols = samples; we'll subset rows 1:i and transpose to matrix for SlicedData
# If covariates are samples in rows, adjust accordingly.

zz <- vector("list", max_PCs)

for (i in seq_len(max_PCs)) {
  cov_subset <- cov_all[1:i, , drop = FALSE]
  cvrt <- SlicedData$new()
  cvrt$CreateFromMatrix(as.matrix(cov_subset))
  
  eQTL <- Matrix_eQTL_main(
    snps = snps.X,
    gene = genes,
    cvrt = cvrt,
    useModel = modelLINEAR,
    pvOutputThreshold = 0,
    pvOutputThreshold.cis = 1e-3,
    snpspos = as.data.frame(snpPos),
    genepos = as.data.frame(genePos)
  )
  
  cis <- as.data.frame(eQTL$cis$eqtls)
  sig <- cis[cis$FDR < 0.05, , drop = FALSE]
  
  zz[[i]] <- data.frame(
    ene.PC = i,
    paire = nrow(sig),
    uniq.snp = length(unique(as.character(sig$snps))),
    uniq.gene = length(unique(as.character(sig$gene))),
    stringsAsFactors = FALSE
  )
}

zz_df <- do.call(rbind, zz)
fwrite(zz_df, out_summary)
