#!/usr/bin/env Rscript
# interaction.analysis-Day7vsDay2.r - test SNP x TB_Status_day7 vs day2 interaction with lmerTest
# Usage: Rscript interaction.analysis-Day7vsDay2.r <chunk_number>

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(lmerTest)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript run_interaction_lmm_chunk.R <chunk_number>")

chunk_number <- as.integer(args[1])

# ---- user-editable paths (adjust for your repo) ----
input_dir   <- "../input"
leads_file  <- "../conditional_analysis/backwards_threshod.pairs_lmerTest_ALL.txt"
rna_file    <- file.path(input_dir, "GE_415.txt.gz")
meta_file   <- file.path(input_dir, "meta.data_415.txt")
pcs_file    <- file.path(input_dir, "Covariates.PCs_415.txt")
dosage_dir  <- input_dir  # contains CHR.<chr>.SNP_rsid_ALL.267.TB.txt
out_prefix  <- paste0("interaction_chunk", chunk_number)
no_of_chunks <- 100

# ---- load inputs ----
RNAmatrixLM <- as.data.frame(fread(rna_file, stringsAsFactors = FALSE))
rownames(RNAmatrixLM) <- RNAmatrixLM$gene_name
RNAmatrixLM$gene_name <- NULL

MetaTable <- as.data.frame(fread(meta_file, stringsAsFactors = FALSE))

covariatesPCsLM <- fread(pcs_file, stringsAsFactors = FALSE)
# build covariate names (columns 2..47 as before)
pc_cols <- colnames(covariatesPCsLM)[2:47]
# Create covariate formula part (we will use column names, but put covariates into the model data frame)
cov_formula_terms <- pc_cols

# ---- load & prepare leads ----
leads <- fread(leads_file, stringsAsFactors = FALSE)
leads <- leads[!is.na(leads$snp), ]
leads$pair <- paste0(leads$snp, "-", leads$Gene)
# parse SNP into chr,pos,ref,alt if needed
tmp <- reshape::colsplit(leads$snp, ":", names = c("chr", "pos_hg19", "ref", "alt"))
leads[, c("chr","pos_hg19","ref","alt")] <- tmp

# split into chunks
f <- ceiling(seq_len(nrow(leads)) / nrow(leads) * no_of_chunks)
chunks <- split(leads, f)
if (is.na(chunk_number) || chunk_number < 1 || chunk_number > length(chunks)) stop("chunk_number out of range")
selected_chunk <- as.data.frame(chunks[[chunk_number]])
message("Processing chunk ", chunk_number, " with ", nrow(selected_chunk), " pairs")

# ---- helper: safe coefficient extraction for interaction ----
extract_interaction <- function(coef_table, interaction_prefix = "snp:TB_Status") {
  if (is.null(coef_table) || nrow(coef_table) == 0) return(list(beta=NA, se=NA, p=NA))
  rn <- rownames(coef_table)
  # find a row that includes the interaction term; be permissive about naming patterns
  idx <- grep(paste0("^", interaction_prefix), rn)               # exact start match
  if (length(idx) == 0) idx <- grep("snp:.*TB_Status", rn)       # broader patterns
  if (length(idx) == 0) idx <- grep("snp\\:.*TB_Status", rn)    # alternate escapes
  if (length(idx) == 0) idx <- grep("snp.*TB_Status", rn)        # fallback
  if (length(idx) == 0) return(list(beta=NA, se=NA, p=NA))
  row <- coef_table[idx[1], , drop = FALSE]
  # coefficients: Estimate, Std. Error, t value, Pr(>|t|)
  beta <- as.numeric(row[1, 1])
  se   <- as.numeric(row[1, 2])
  pval <- as.numeric(row[1, 4])  # lmerTest summary places p-value in 4th column
  return(list(beta=beta, se=se, p=pval))
}

# ---- loop through selected pairs ----
start.time <- Sys.time()
results <- vector("list", nrow(selected_chunk))

for (i in seq_len(nrow(selected_chunk))) {
  pair_row <- selected_chunk[i, , drop = FALSE]
  snp_id <- as.character(pair_row$snp)
  gene_id <- as.character(pair_row$Gene)
  chr_id <- as.character(pair_row$chr)
  message(sprintf("[%d/%d] Pair: %s | Gene: %s | CHR: %s", i, nrow(selected_chunk), pair_row$pair, gene_id, chr_id))
  
  # load genotype dosage for the chromosome (rows=variants, cols=samples expected)
  geno_file <- file.path(dosage_dir, paste0("CHR.", chr_id, ".SNP_rsid_ALL.267.TB.txt"))
  if (!file.exists(geno_file)) {
    warning("Genotype file not found: ", geno_file, " — skipping pair.")
    results[[i]] <- data.frame(gene = gene_id, snp = snp_id, Estimate.int.day7 = NA, se.int.day7 = NA, P.int.day7 = NA, stringsAsFactors = FALSE)
    next
  }
  dosagefiltered_DR <- as.data.frame(fread(geno_file, stringsAsFactors = FALSE))
  rownames(dosagefiltered_DR) <- dosagefiltered_DR$ID
  dosagefiltered_DR$ID <- NULL
  
  genotype <- as.data.frame(t(dosagefiltered_DR))
  genotype$genotype.ID <- rownames(genotype)
  
  # merge with MetaTable to attach RNA sequencing IDs
  merged <- merge(MetaTable, genotype, by.x = "genotype.ID_corrected", by.y = "genotype.ID", all = TRUE)
  # set rownames as RNA sequencing ID (must exist)
  if (!"RNA.Sequencing_ID" %in% colnames(merged)) {
    warning("MetaTable lacks RNA.Sequencing_ID column; skipping pair.")
    results[[i]] <- data.frame(gene = gene_id, snp = snp_id, Estimate.int.day7 = NA, se.int.day7 = NA, P.int.day7 = NA, stringsAsFactors = FALSE)
    next
  }
  rownames(merged) <- merged$RNA.Sequencing_ID
  genotype2 <- as.data.frame(t(merged[ , -seq_len( which(colnames(merged) == "RNA.Sequencing_ID") ) , drop = FALSE])) # keep consistent behavior
  # Reorder genotype columns to match RNAmatrixLM samples
  reorder_idx <- match(colnames(RNAmatrixLM), colnames(genotype2))
  if (any(is.na(reorder_idx))) {
    warning("Sample ID mismatch between expression and genotype files; skipping pair.")
    results[[i]] <- data.frame(gene = gene_id, snp = snp_id, Estimate.int.day7 = NA, se.int.day7 = NA, P.int.day7 = NA, stringsAsFactors = FALSE)
    next
  }
  genotype2_reorder <- genotype2[, reorder_idx, drop = FALSE]
  
  # sanity: ensure covariates and expression sample order match
  if (!all(colnames(RNAmatrixLM) == covariatesPCsLM$RNA.Sequencing_ID)) {
    stop("Covariate sample order does not match expression samples. Fix inputs.")
  }
  
  # prepare vectors
  if (!(gene_id %in% rownames(RNAmatrixLM))) {
    warning("Gene not found in expression matrix: ", gene_id, "; skipping")
    results[[i]] <- data.frame(gene = gene_id, snp = snp_id, Estimate.int.day7 = NA, se.int.day7 = NA, P.int.day7 = NA, stringsAsFactors = FALSE)
    next
  }
  gene_vec <- as.numeric(RNAmatrixLM[gene_id, ])
  snp_vec  <- as.numeric(genotype2_reorder[snp_id, ])
  
  # Build model data.frame containing gene, snp, donor, TB_Status, covariates
  model_df <- data.frame(
    gene = gene_vec,
    snp  = snp_vec,
    donor = MetaTable$Sample_ID,
    TB_Status = MetaTable$TB_Status,
    stringsAsFactors = FALSE
  )
  # Attach covariates (PCs) from covariatesPCsLM aligned by RNA.Sequencing_ID
  cov_mat <- covariatesPCsLM[, ..pc_cols]
  # covariatesPCsLM rows correspond to samples in the same order as RNAmatrixLM (asserted above)
  cov_df <- as.data.frame(cov_mat)
  # ensure same column names for covariates in model_df
  colnames(cov_df) <- pc_cols
  model_df <- cbind(model_df, cov_df)
  
  # factorise TB_Status and relevel
  model_df$TB_Status <- as.factor(model_df$TB_Status)
  if (!("Latent TB" %in% levels(model_df$TB_Status))) {
    # if "Latent TB" not present, keep first level as reference
    model_df$TB_Status <- relevel(model_df$TB_Status, ref = levels(model_df$TB_Status)[1])
  } else {
    model_df$TB_Status <- relevel(model_df$TB_Status, ref = "Latent TB")
  }
  
  # build formula: gene ~ snp * TB_Status + covariates + (1|donor)
  fixed_terms <- paste(c("snp * TB_Status", cov_formula_terms), collapse = " + ")
  form_str <- paste("gene ~", fixed_terms, "+ (1 | donor)")
  form <- as.formula(form_str)
  
  # fit model with tryCatch to avoid crashing whole run
  fit <- tryCatch({
    suppressMessages(lmer(form, data = model_df))
  }, error = function(e) {
    warning("lmer failed for pair ", snp_id, "-", gene_id, ": ", e$message)
    return(NULL)
  })
  
  if (is.null(fit)) {
    results[[i]] <- data.frame(gene = gene_id, snp = snp_id, Estimate.int.day7 = NA, se.int.day7 = NA, P.int.day7 = NA, stringsAsFactors = FALSE)
    next
  }
  
  ssum <- summary(fit)
  coef_table <- as.matrix(ssum$coefficients)  # columns: Estimate, Std. Error, t value, Pr(>|t|)
  
  # Determine interaction term name dynamically (non-reference TB_Status level)
  tb_levels <- levels(model_df$TB_Status)
  ref_level <- tb_levels[1]
  nonref_levels <- setdiff(tb_levels, ref_level)
  # choose the first non-reference level (Day7 in your dataset)
  interaction_term <- if (length(nonref_levels) >= 1) paste0("snp:TB_Status", nonref_levels[1]) else "snp:TB_Status"
  # extract robustly
  out <- extract_interaction(coef_table, interaction_prefix = "snp:TB_Status")
  # If missing but there are multiple possible naming patterns, try explicit pattern
  if (is.na(out$beta)) {
    out <- extract_interaction(coef_table, interaction_prefix = paste0("snp:TB_Status", nonref_levels[1]))
  }
  
  results[[i]] <- data.frame(
    gene = gene_id,
    snp  = snp_id,
    Estimate.int.day7 = out$beta,
    se.int.day7       = out$se,
    P.int.day7        = out$p,
    stringsAsFactors = FALSE
  )
} # end for loop

# ---- bind and write ----
res_df <- do.call(rbind, results)
out_file <- gzfile(paste0(out_prefix, ".txt.gz"), "w")
write.table(res_df, file = out_file, quote = FALSE, row.names = FALSE, sep = "\t")
close(out_file)

end.time <- Sys.time()
message("Done. Time elapsed: ", round(as.numeric(difftime(end.time, start.time, units = "mins")), 2), " minutes")
