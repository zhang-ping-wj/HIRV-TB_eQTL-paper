#!/usr/bin/env Rscript
# lmm.R — run LMM model for cis-eQTL mapping
# Usage: Rscript run_LMM_chr.R <CHR> <NUM_PCS>

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(lme4)
  library(lmerTest)
  library(foreach)
  library(doMC)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript run_LMM_chr.R <CHR> <NUM_PCS>")
}
CHR <- args[1]          # chromosome identifier (string or number)
NUM_PCS <- as.integer(args[2])

# ----- user-editable base paths (edit for your repo) -----
input_dir   <- "input"   # directory containing prepared input files
out_prefix  <- "LMM-output" 

# input files (constructed from CHR / input_dir)
rna_file    <- file.path(input_dir, "GE_415.txt.gz")                   # expression matrix
meta_file   <- file.path(input_dir, "meta.data_415.txt")              # metadata
cov_file    <- file.path(input_dir, "Covariates.PCs_415.txt")         # PCs file
dosage_file <- file.path(input_dir, paste0("CHR.", CHR, ".SNP_rsid_ALL.267.TB.txt")) # genotype dosage (tab-delim)
pairs_file  <- file.path(input_dir, paste0("snp.gene.pairs.CHR.", CHR, ".txt.gz"))   # SNP-gene pairs

# ----- load inputs -----
RNAmatrixLM <- as.data.frame(fread(rna_file))
rownames(RNAmatrixLM) <- RNAmatrixLM$gene_name
RNAmatrixLM$gene_name <- NULL

MetaTable <- as.data.frame(fread(meta_file))
donor <- as.factor(MetaTable$Sample_ID)

covariatesPCsLM <- fread(cov_file)

# Build fixed effects string from PCs (columns 2..NUM_PCS)
pc_cols <- colnames(covariatesPCsLM)[2:NUM_PCS]
coFixedLM <- paste(paste0("covariatesPCsLM$", pc_cols), collapse = " + ")

# dosage/genotype: rows = SNP IDs in column 'ID', columns = samples
dosagefiltered_DR <- as.data.frame(fread(dosage_file, stringsAsFactors = FALSE))
rownames(dosagefiltered_DR) <- dosagefiltered_DR$ID
dosagefiltered_DR$ID <- NULL

# transpose genotype so rows = samples, cols = SNPs, then merge with metadata to align sample IDs
genotype <- as.data.frame(t(dosagefiltered_DR))
genotype$genotype.ID <- rownames(genotype)

# Merge metadata and genotype using the fields present in your meta file.
# Here we expect MetaTable to have a column named 'genotype.ID_corrected' and 'RNA.Sequencing_ID'
xxxx <- merge(MetaTable, genotype, by.x = "genotype.ID_corrected", by.y = "genotype.ID", all = TRUE)
rownames(xxxx) <- xxxx$RNA.Sequencing_ID

# Reconstruct genotype matrix with RNAseq sample order
genotype2 <- as.data.frame(t(xxxx[ , -seq_len(16), drop = FALSE]))  # keep behavior from original script
reorder_idx <- match(colnames(RNAmatrixLM), colnames(genotype2))
genotype2_reorder <- genotype2[, reorder_idx, drop = FALSE]

# load SNP-gene pairs
Pairs <- as.data.frame(fread(pairs_file))

# simple sanity checks (silent unless fails)
stopifnot(ncol(RNAmatrixLM) == ncol(genotype2_reorder))
stopifnot(all(colnames(RNAmatrixLM) == covariatesPCsLM$RNA.Sequencing_ID))

# ----- parallel setup -----
n.cpus_env <- Sys.getenv("SLURM_CPUS_PER_TASK", unset = "")
if (nchar(n.cpus_env)) {
  n.cpus <- as.numeric(n.cpus_env)
} else {
  n.cpus <- parallel::detectCores(logical = FALSE)
  if (is.na(n.cpus) || n.cpus < 1) n.cpus <- 1
}
registerDoMC(cores = n.cpus)

# ----- run LMM per SNP-gene pair in parallel -----
start.time <- Sys.time()
t <- nrow(Pairs)

# Export only required objects to workers
exports <- c("RNAmatrixLM", "genotype2_reorder", "Pairs", "coFixedLM", "donor")

dopar <- foreach(i = seq_len(t), .combine = rbind,
                 .packages = c("lme4", "lmerTest"),
                 .export = exports) %dopar% {
                   gene_vec <- as.numeric(RNAmatrixLM[Pairs[i, 2], ])
                   snp_vec  <- as.numeric(genotype2_reorder[Pairs[i, 1], ])
                   form1 <- as.formula(paste("gene_vec ~ snp_vec +", coFixedLM, "+ (1|donor)"))
                   
                   # fit model (suppress messages), extract SNP effect row
                   fit <- suppressMessages(lmer(form1))
                   s   <- summary(fit)
                   coef_row <- s$coefficients[2, , drop = TRUE]  # [Estimate, Std.Error, t value, Pr(>|t|)]
                   data.frame(
                     gene_name = Pairs[i, 2],
                     snp       = Pairs[i, 1],
                     Estimate  = as.numeric(coef_row[1]),
                     se        = as.numeric(coef_row[2]),
                     t.value   = as.numeric(coef_row[3]),
                     p         = as.numeric(coef_row[4]),
                     stringsAsFactors = FALSE
                   )
                 }

end.time <- Sys.time()
message("Time taken: ", signif(as.numeric(difftime(end.time, start.time, units = "secs")), 4), " secs")

# ----- write results -----
out_file <- gzfile(paste0(out_prefix, "-chr", CHR, "_PC", NUM_PCS, ".txt.gz"), "w")
write.table(dopar, file = out_file, quote = FALSE, row.names = FALSE, sep = "\t")
close(out_file)
