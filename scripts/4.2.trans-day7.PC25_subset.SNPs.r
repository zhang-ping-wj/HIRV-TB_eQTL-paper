#!/usr/bin/env Rscript
# run_matrixeqtl_day7.R — trans-eQTL mapping using day7 lead cis-eQTLs
# Usage:
#   Rscript run_matrixeqtl_day7.R

suppressPackageStartupMessages({
  library(MatrixEQTL)
  library(data.table)
  library(dplyr)
})

# ---------- user-editable paths/params ----------
input_dir        <- "./1.TB.day7.input"
snps_file        <- file.path(input_dir, "SNP_rsid_day7.latent.TB.txt.gz")
genes_file       <- file.path(input_dir, "GE.txt")
snpPos_file      <- file.path(input_dir, "snpsloc_rsid_day7.latent.TB.txt.gz")
genePos_file     <- file.path(input_dir, "geneloc.txt")
cov_file         <- file.path(input_dir, "Covariates.all.PCs.txt")

# external resources (set to local copies; do NOT commit private /well paths)
egene_lead_file  <- "permute_full_d7.txt.gz"   # path to fastQTL permutation results (variant_id, ...)
gene_name_file   <- "gencode.v31lift37_gene.txt" # gencode gene metadata (no pseudogenes)
tst_file         <- "./AllTSTResponseGenes.csv"   # TST gene list

# outputs
out_trans_all    <- "Day7_trans_subset.day7cis_txt.gz"
out_trans_local  <- "Day7_trans_subset.day7cis_locus.sig.txt.gz"

# analysis parameters
pcs_to_use       <- 25
cis_distance     <- 5e6

# ---------- load input data ----------
snps_all <- fread(snps_file)
genes    <- SlicedData$new(); genes$LoadFile(genes_file)

snpPos   <- fread(snpPos_file, stringsAsFactors = FALSE)
genePos  <- fread(genePos_file, stringsAsFactors = FALSE)

# ---------- select candidate SNPs from eGene leads ----------
eGene.lead <- fread(egene_lead_file)
setnames(eGene.lead, c("gene_name","num_var","beta_shape1","beta_shape2","Dummy",
                       "variant_id","tss_distance","pval_nominal","slope","pval_perm","pval_beta"))
eGene.lead[, p_fdr := p.adjust(pval_beta, method = "fdr")]
eGene.lead <- eGene.lead[!is.na(p_fdr) & p_fdr < 0.05]

# filter to TST genes
TST <- fread(tst_file)
TST <- TST[MaxFC >= 0.58]
eGene.lead.TST <- eGene.lead[gene_name %in% TST$HGNC]

# subset SNPs to lead variants
filtered_snps <- snps_all[ID %in% eGene.lead.TST$variant_id]

# ---------- build SlicedData for SNPs ----------
snps.X <- SlicedData$new()
snps.X$CreateFromMatrix(as.matrix(filtered_snps[ , -1, with = FALSE]))
rownames(snps.X) <- filtered_snps$ID

# ---------- filter genes (remove pseudogenes) ----------
gene.name <- fread(gene_name_file)
gene.name <- gene.name[!V6 %like% "pseudogene", ]
keep_genes <- gene.name$V7

keep_indices <- which(rownames(genes) %in% keep_genes)
genes$RowReorder(keep_indices)

# ---------- covariates ----------
Covariates_values <- read.table(cov_file, row.names = 1, header = TRUE)
Covariates_values <- Covariates_values[1:pcs_to_use, , drop = FALSE]
cvrt <- SlicedData$new(); cvrt$CreateFromMatrix(as.matrix(Covariates_values))

# ---------- Matrix eQTL run ----------
eQTL <- Matrix_eQTL_main(
  snps = snps.X,
  gene = genes,
  cvrt = cvrt,
  useModel = modelLINEAR,
  pvOutputThreshold = 1,
  pvOutputThreshold.cis = 1,
  cisDist = cis_distance,
  snpspos = as.data.frame(snpPos),
  genepos = as.data.frame(genePos)
)

# ---------- post-processing ----------
eQTL$all$eqtls$beta_se <- eQTL$all$eqtls$beta / eQTL$all$eqtls$statistic
trans <- as.data.frame(eQTL$trans$eqtls)

# write full trans results
fwrite(trans, out_trans_all, compress = "gzip", sep = "\t")

# ---------- SNP-level (local) BH correction ----------
colnames(trans) <- c("snp","Gene","statistic","pvalue","FDR","beta")
transQTLs <- trans %>%
  group_by(snp) %>%
  mutate(adjusted.P.local = p.adjust(pvalue, method = "BH")) %>%
  ungroup()

transQTLs_adj.p <- transQTLs[transQTLs$adjusted.P.local <= 0.05, ]

# write local-significant hits
fwrite(transQTLs_adj.p, out_trans_local, compress = "gzip", sep = "\t")
