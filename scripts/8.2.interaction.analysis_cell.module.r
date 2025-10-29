#!/usr/bin/env Rscript
# interaction.analysis_cell.module.r - test SNP x cell module interaction analysis with lmerTest
# Usage: Rscript interaction.analysis_cell.module.r <chunk_number>
#!/usr/bin/env Rscript
# interaction_chunk.R — quiet LMM interaction scan using plain NA

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(lmerTest)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript interaction_chunk.R <chunk_number>")
chunk_number <- as.integer(args[1])

# ---- user paths ----
rna_file   <- "../input/GE_415.txt.gz"
meta_file  <- "./meta.data_415_module.included.txt"
pcs_file   <- "../input/Covariates.PCs_415.txt"
leads_file <- "../conditional_analysis/backwards_threshod.pairs_lmerTest_ALL.txt"
dosage_dir <- "../input"
out_file   <- paste0("interaction_chunk", chunk_number, ".txt")
no_of_chunks <- 100

# ---- load data ----
RNAmatrixLM <- as.data.frame(fread(rna_file, stringsAsFactors = FALSE))
rownames(RNAmatrixLM) <- RNAmatrixLM$gene_name
RNAmatrixLM$gene_name <- NULL

MetaTable <- as.data.frame(fread(meta_file, stringsAsFactors = FALSE))
covariatesPCsLM <- fread(pcs_file, stringsAsFactors = FALSE)

leads <- fread(leads_file, stringsAsFactors = FALSE)
leads <- leads[!is.na(leads$snp), ]
leads$pair <- paste0(leads$snp, "-", leads$Gene)
parsed <- reshape::colsplit(leads$snp, ":", names = c("chr", "pos_hg19", "ref", "alt"))
leads[, c("chr","pos_hg19","ref","alt")] <- parsed

# chunking
f <- ceiling(seq_len(nrow(leads)) / nrow(leads) * no_of_chunks)
chunked <- split(leads, f)
if (is.na(chunk_number) || chunk_number < 1 || chunk_number > length(chunked)) stop("chunk_number out of range")
selected_chunk <- as.data.frame(chunked[[chunk_number]])

interaction_terms <- c("CD4_T", "CD8_T", "NK", "antimicrobial", "Ag_presenting",
                       "pan_myeloid", "pan_T", "antimicrobial.vs.CD4_T", "antimicrobial.vs.NK")

pc_cols <- colnames(covariatesPCsLM)[2:47]
cov_df_template <- as.data.frame(covariatesPCsLM[, ..pc_cols])
colnames(cov_df_template) <- pc_cols

results_list <- vector("list", nrow(selected_chunk))

for (i in seq_len(nrow(selected_chunk))) {
  pair_row <- selected_chunk[i, , drop = FALSE]
  snp_id <- as.character(pair_row$snp)
  gene_id <- as.character(pair_row$Gene)
  chr_id <- as.character(pair_row$chr)
  
  geno_file <- file.path(dosage_dir, paste0("CHR.", chr_id, ".SNP_rsid_ALL.267.TB.txt"))
  if (!file.exists(geno_file)) {
    results_list[[i]] <- data.frame(Gene = gene_id, snp = snp_id,
                                    matrix(NA, nrow = 1, ncol = length(interaction_terms) * 3,
                                           dimnames = list(NULL, unlist(lapply(interaction_terms, function(t) paste0(c("Estimate.", "SE.", "P."), t))))))
    next
  }
  dosage <- as.data.frame(fread(geno_file, stringsAsFactors = FALSE))
  rownames(dosage) <- dosage$ID
  dosage$ID <- NULL
  genotype <- as.data.frame(t(dosage))
  genotype$genotype.ID <- rownames(genotype)
  
  merged <- merge(MetaTable, genotype, by.x = "genotype.ID_corrected", by.y = "genotype.ID", all = TRUE)
  if (!("RNA.Sequencing_ID" %in% colnames(merged))) {
    results_list[[i]] <- data.frame(Gene = gene_id, snp = snp_id,
                                    matrix(NA, nrow = 1, ncol = length(interaction_terms) * 3,
                                           dimnames = list(NULL, unlist(lapply(interaction_terms, function(t) paste0(c("Estimate.", "SE.", "P."), t))))))
    next
  }
  rownames(merged) <- merged$RNA.Sequencing_ID
  genotype2 <- as.data.frame(t(merged[, -c(1:16), drop = FALSE]))
  
  reorder_idx <- match(colnames(RNAmatrixLM), colnames(genotype2))
  if (any(is.na(reorder_idx))) {
    results_list[[i]] <- data.frame(Gene = gene_id, snp = snp_id,
                                    matrix(NA, nrow = 1, ncol = length(interaction_terms) * 3,
                                           dimnames = list(NULL, unlist(lapply(interaction_terms, function(t) paste0(c("Estimate.", "SE.", "P."), t))))))
    next
  }
  genotype2_reorder <- genotype2[, reorder_idx, drop = FALSE]
  
  if (!(gene_id %in% rownames(RNAmatrixLM)) || !(snp_id %in% rownames(genotype2_reorder))) {
    results_list[[i]] <- data.frame(Gene = gene_id, snp = snp_id,
                                    matrix(NA, nrow = 1, ncol = length(interaction_terms) * 3,
                                           dimnames = list(NULL, unlist(lapply(interaction_terms, function(t) paste0(c("Estimate.", "SE.", "P."), t))))))
    next
  }
  
  gene_vec <- as.numeric(RNAmatrixLM[gene_id, ])
  snp_vec  <- as.numeric(genotype2_reorder[snp_id, ])
  
  model_df <- data.frame(gene = gene_vec, snp = snp_vec, donor = MetaTable$Sample_ID, stringsAsFactors = FALSE)
  model_df <- cbind(model_df, cov_df_template)
  
  term_results <- lapply(interaction_terms, function(term) {
    if (!term %in% colnames(MetaTable)) {
      return(c(Estimate = NA, SE = NA, P = NA))
    }
    model_df[[term]] <- MetaTable[[term]]
    fixed_part <- paste0("snp * ", term, " + ", paste0(pc_cols, collapse = " + "))
    form <- as.formula(paste("gene ~", fixed_part, "+ (1|donor)"))
    fit <- tryCatch(suppressMessages(lmer(form, data = model_df)), error = function(e) NULL)
    if (is.null(fit)) return(c(Estimate = NA, SE = NA, P = NA))
    coef_tbl <- summary(fit)$coefficients
    term_pat <- paste0("^snp[:*]?", term, "$")
    rn <- rownames(coef_tbl)
    idx <- grep(term_pat, rn, ignore.case = TRUE)
    if (length(idx) == 0) idx <- grep(paste0("snp.*", term), rn, ignore.case = TRUE)
    if (length(idx) == 0) return(c(Estimate = NA, SE = NA, P = NA))
    row <- coef_tbl[idx[1], ]
    c(Estimate = as.numeric(row[1]), SE = as.numeric(row[2]), P = as.numeric(row[4]))
  })
  
  flat <- unlist(lapply(term_results, function(x) x))
  names(flat) <- unlist(lapply(interaction_terms, function(t) paste0(c("Estimate.", "SE.", "P."), t)))
  results_list[[i]] <- c(Gene = gene_id, snp = snp_id, as.list(flat))
}

res_df <- rbindlist(lapply(results_list, function(x) as.data.table(as.list(x))), fill = TRUE)
fwrite(res_df, out_file, sep = "\t", quote = FALSE, na = "NA")
