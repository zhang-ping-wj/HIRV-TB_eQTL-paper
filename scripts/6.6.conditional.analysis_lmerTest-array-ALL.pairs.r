#!/usr/bin/env Rscript
# conditional.analysis_lmerTest-array-ALL.pairs.r - Stepwise forward/backward conditional mapping (lmerTest)
# Usage: Rscript conditional_mapping_stepwise.R <chunk_number>

suppressPackageStartupMessages({
  library(data.table)
  library(lmerTest)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript conditional_mapping_stepwise.R <chunk_number>")

chunk_number <- as.integer(args[1])

# ---- user-editable paths and parameters ----
input_dir    <- "../input"                        # base input directory
rna_file     <- file.path(input_dir, "GE_415.txt.gz")
meta_file    <- file.path(input_dir, "meta.data_415.txt")
pcs_file     <- file.path(input_dir, "Covariates.PCs_415.txt")
dosage_dir   <- file.path(input_dir)             # where CHR.<chr>.SNP... files live
pairs_prefix <- file.path(input_dir, "snp.gene.pairs.CHR.") # add CHR and .txt.gz
eigen_file   <- "eigenMT.result_final_thresholds.txt.gz"
nominal_thr  <- "nominal_pval_thresholds.txt"    # table with columns 'gene' and 'threshold'
chunks_total <- 200

# ---- load shared resources ----
RNAmatrixLM <- as.data.frame(fread(rna_file, stringsAsFactors = FALSE))
rownames(RNAmatrixLM) <- RNAmatrixLM$gene_name; RNAmatrixLM$gene_name <- NULL

MetaTable <- as.data.frame(fread(meta_file, stringsAsFactors = FALSE))

covariatesPCsLM <- fread(pcs_file, stringsAsFactors = FALSE)
# build fixed-effect string using PCs columns 2..47 (as in your original)
pc_cols <- colnames(covariatesPCsLM)[2:47]
coFixedLM <- paste(paste0("covariatesPCsLM$", pc_cols), collapse = " + ")

# ---- thresholds & chunking ----
thresholds.eQTLs <- fread(eigen_file, stringsAsFactors = FALSE)
# expect SNP column like "chr:pos:ref:alt" in thresholds.eQTLs$SNP
thresholds.eQTLs[, c("chr","pos_hg19","ref","alt") := tstrsplit(SNP, ":", fixed = TRUE)]
thresholds.eQTLs <- thresholds.eQTLs[qval <= 0.05]

# split into chunks (preserve order)
nrows <- nrow(thresholds.eQTLs)
fvec <- ceiling(seq_len(nrows) / nrows * chunks_total)
chunk_list <- split(thresholds.eQTLs, fvec)

if (chunk_number < 1 || chunk_number > length(chunk_list)) stop("chunk_number out of range")
selected_chunk <- chunk_list[[chunk_number]]
message("Processing chunk ", chunk_number, " containing ", length(unique(selected_chunk$gene)), " genes")

# load nominal thresholds file
thresholds_tbl <- fread(nominal_thr, stringsAsFactors = FALSE)

# ---- helper to load genotype matrix for a chromosome ----
load_genotype_matrix <- function(chr) {
  geno_file <- file.path(dosage_dir, paste0("CHR.", chr, ".SNP_rsid_ALL.267.TB.txt"))
  if (!file.exists(geno_file)) stop("Genotype file missing: ", geno_file)
  dosagefiltered_DR <- as.data.frame(fread(geno_file, stringsAsFactors = FALSE))
  rownames(dosagefiltered_DR) <- dosagefiltered_DR$ID
  dosagefiltered_DR$ID <- NULL
  genotype <- as.data.frame(t(dosagefiltered_DR))
  genotype$genotype.ID <- rownames(genotype)
  return(genotype)
}

# ---- iterate per gene in the selected chunk ----
for (target_gene in unique(selected_chunk$gene)) {
  message("Processing gene: ", target_gene)
  this.gene <- selected_chunk[selected_chunk$gene == target_gene, ]
  CHR <- unique(this.gene$chr)
  if (length(CHR) != 1) {
    warning("Multiple chromosomes for gene ", target_gene, "; skipping")
    next
  }
  # gene-level threshold
  thr_val <- thresholds_tbl$threshold[match(target_gene, thresholds_tbl$gene)]
  if (is.na(thr_val)) {
    warning("No nominal threshold for gene ", target_gene, "; skipping")
    next
  }
  threshold <- thr_val
  
  # initial significant SNP(s) from selected_chunk (lead hits)
  sig.snps <- this.gene[, .(SNP, gene, beta, `t-stat`, `p-value`)]
  setnames(sig.snps, c("SNP","gene","beta","t-stat","p-value"), c("snp","Gene","eQTL_beta","eQTL_t","pvalue"))
  sig.snps$pvalue <- as.numeric(sig.snps$pvalue)
  
  # load pairs for the chromosome and restrict to this gene
  pairs_file_chr <- paste0(pairs_prefix, CHR, ".txt.gz")
  if (!file.exists(pairs_file_chr)) {
    warning("Pairs file not found: ", pairs_file_chr); next
  }
  pairs <- as.data.frame(fread(pairs_file_chr))
  pairs.eqtl <- pairs[pairs$gene %in% target_gene, , drop = FALSE]
  if (nrow(pairs.eqtl) == 0) {
    warning("No pairs for gene ", target_gene); next
  }
  pairs.eqtl$snp <- pairs.eqtl$snps
  
  # load genotype matrix for this chromosome and align to samples
  genotype <- load_genotype_matrix(CHR)
  merged <- merge(MetaTable, genotype, by.x = "genotype.ID_corrected", by.y = "genotype.ID", all = TRUE)
  rownames(merged) <- merged$RNA.Sequencing_ID
  genotype2 <- as.data.frame(t(merged[ , -seq_len(16), drop = FALSE]))
  reorder_idx <- match(colnames(RNAmatrixLM), colnames(genotype2))
  genotype2_reorder <- genotype2[, reorder_idx, drop = FALSE]
  
  # basic checks
  if (!all(colnames(RNAmatrixLM) == colnames(genotype2_reorder))) {
    stop("Sample ordering mismatch for gene ", target_gene)
  }
  if (!all(colnames(RNAmatrixLM) == covariatesPCsLM$RNA.Sequencing_ID)) {
    stop("Covariate sample ordering mismatch")
  }
  
  # prepare candidate SNP list excluding lead(s)
  pairs.eqtl.it <- pairs.eqtl[!(pairs.eqtl$snp %in% sig.snps$snp), , drop = FALSE]
  if (nrow(pairs.eqtl.it) == 0) {
    message("No candidate SNPs left for forward selection for gene ", target_gene)
    next
  }
  
  # forward selection loop
  last.was.sig <- min(sig.snps$pvalue, na.rm = TRUE) <= threshold
  while (isTRUE(last.was.sig)) {
    pairs.eqtl.it <- pairs.eqtl[!(pairs.eqtl$snp %in% sig.snps$snp), , drop = FALSE]
    if (nrow(pairs.eqtl.it) == 0) break
    results_list <- vector("list", nrow(pairs.eqtl.it))
    for (ii in seq_len(nrow(pairs.eqtl.it))) {
      snp_id <- pairs.eqtl.it$snp[ii]
      gene_vec <- as.numeric(RNAmatrixLM[target_gene, ])
      snp_vec  <- as.numeric(genotype2_reorder[snp_id, ])
      # convert existing sig.snps genotypes into matrix for inclusion
      if (nrow(sig.snps) > 0) {
        siggeno <- genotype2_reorder[rownames(genotype2_reorder) %in% sig.snps$snp, , drop = FALSE]
        # when there are multiple sig.snps, create matrix form; else a vector
        if (nrow(siggeno) == 0) siggeno_mat <- NULL else siggeno_mat <- apply(siggeno, 1, as.numeric)
      } else {
        siggeno_mat <- NULL
      }
      # check for NA
      if (any(is.na(snp_vec))) {
        results_list[[ii]] <- data.frame(Gene = target_gene, snp = snp_id, eQTL_beta = NA, se = NA, eQTL_t = NA, pvalue = 1)
        next
      }
      # build formula dynamically: include siggeno as "sig.snp" placeholder
      if (!is.null(siggeno_mat)) {
        # create a combined data.frame for lmer
        # bind columns: gene, snp, sig.snp.* , covariates, donor
        # easier to use formula with a combined data.frame in environment
        df <- data.frame(gene = gene_vec, snp = snp_vec, donor = MetaTable$Sample_ID, stringsAsFactors = FALSE)
        # add siggeno columns
        if (is.matrix(siggeno_mat)) {
          siggenonames <- paste0("sigSNP", seq_len(ncol(siggeno_mat)))
          df[siggenonames] <- as.data.frame(t(siggeno_mat))
          sig_term <- paste(siggenonames, collapse = " + ")
        } else {
          # single existing signal
          df[["sigSNP1"]] <- as.numeric(siggeno_mat)
          sig_term <- "sigSNP1"
        }
        # bind covariates (PCs) columns from covariatesPCsLM matching sample order
        cov_df <- as.data.frame(t(covariatesPCsLM[, -1, with = FALSE]))
        colnames(cov_df) <- covariatesPCsLM$RNA.Sequencing_ID
        cov_df <- as.data.frame(t(covariatesPCsLM[, -1, with = FALSE])) # fallback — kept minimal
        form_str <- paste0("gene ~ snp + ", sig_term, " + ", coFixedLM, " + (1|donor)")
      } else {
        form_str <- paste0("gene ~ snp + ", coFixedLM, " + (1|donor)")
        df <- data.frame(gene = gene_vec, snp = snp_vec, donor = MetaTable$Sample_ID, stringsAsFactors = FALSE)
      }
      
      form <- as.formula(form_str)
      # fit model, handle rank-deficiency by returning p=1
      fit_ok <- TRUE
      res_row <- tryCatch({
        fit <- suppressMessages(lmer(form, data = df))
        coef_row <- summary(fit)$coefficients
        if (nrow(coef_row) < 2 || any(is.na(coef_row[2, ]))) {
          data.frame(Gene = target_gene, snp = snp_id, eQTL_beta = NA, se = NA, eQTL_t = NA, pvalue = 1)
        } else {
          data.frame(Gene = target_gene,
                     snp = snp_id,
                     eQTL_beta = as.numeric(coef_row[2, 1]),
                     se = as.numeric(coef_row[2, 2]),
                     eQTL_t = as.numeric(coef_row[2, 3]),
                     pvalue = as.numeric(coef_row[2, 5]),
                     stringsAsFactors = FALSE)
        }
      }, error = function(e) {
        data.frame(Gene = target_gene, snp = snp_id, eQTL_beta = NA, se = NA, eQTL_t = NA, pvalue = 1)
      })
      results_list[[ii]] <- res_row
    } # end per-SNP loop
    
    last.result <- rbindlist(results_list, use.names = TRUE, fill = TRUE)
    last.result[, pvalue := as.numeric(pvalue)]
    min_p <- min(last.result$pvalue, na.rm = TRUE)
    last.was.sig <- (min_p <= threshold)
    if (isTRUE(last.was.sig)) {
      best_row <- last.result[which.min(pvalue), .(snp, Gene, eQTL_beta, eQTL_t, pvalue)]
      sig.snps <- rbind(sig.snps, best_row[, .(snp, Gene, eQTL_beta, eQTL_t, pvalue)])
    } else break
  } # end while forward
  
  # ---- backward selection if more than one signal ----
  if (exists("sig.snps") && nrow(sig.snps) > 0) {
    output_dir <- paste0("lmerTest_results_CHR", CHR)
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    write.table(sig.snps, file = file.path(output_dir, paste0("forward_", target_gene, ".txt")), 
                quote = FALSE, row.names = FALSE, sep = "\t")
    
    n.snps <- nrow(sig.snps)
    if (n.snps == 1) {
      sig.snps.backward <- sig.snps
    } else {
      # perform backward selection: test each candidate by conditioning on the rest
      backward_list <- vector("list", n.snps)
      for (j in seq_len(n.snps)) {
        test.snp <- sig.snps$snp[j]
        other.snps <- sig.snps$snp[-j]
        other_geno <- genotype2_reorder[rownames(genotype2_reorder) %in% other.snps, , drop = FALSE]
        if (nrow(other_geno) == 0) {
          backward_list[[j]] <- data.frame(snp = NA, Gene = target_gene, eQTL_beta = NA, eQTL_t = NA, pvalue = NA, stringsAsFactors = FALSE)
          next
        }
        # build other.snps.geno term similarly to forward (kept concise here)
        # test each candidate SNP against genome-wide pairs excluding other.snps
        pairs.eqtl.it <- pairs.eqtl[!(pairs.eqtl$snp %in% other.snps), , drop = FALSE]
        if (nrow(pairs.eqtl.it) == 0) {
          backward_list[[j]] <- data.frame(snp = NA, Gene = target_gene, eQTL_beta = NA, eQTL_t = NA, pvalue = NA, stringsAsFactors = FALSE)
          next
        }
        # Evaluate best SNP in this backward pass (simple implementation)
        br_results <- vector("list", nrow(pairs.eqtl.it))
        for (ii in seq_len(nrow(pairs.eqtl.it))) {
          snp_id <- pairs.eqtl.it$snp[ii]
          gene_vec <- as.numeric(RNAmatrixLM[target_gene, ])
          snp_vec  <- as.numeric(genotype2_reorder[snp_id, ])
          # other snps as covariates
          other_mat <- apply(other_geno, 1, as.numeric)
          df <- data.frame(gene = gene_vec, snp = snp_vec, donor = MetaTable$Sample_ID, stringsAsFactors = FALSE)
          if (is.matrix(other_mat)) {
            om_names <- paste0("other", seq_len(ncol(other_mat)))
            df[om_names] <- as.data.frame(t(other_mat))
            other_term <- paste(om_names, collapse = " + ")
          } else {
            df[["other1"]] <- as.numeric(other_mat)
            other_term <- "other1"
          }
          form_str <- paste0("gene ~ snp + ", other_term, " + ", coFixedLM, " + (1|donor)")
          form <- as.formula(form_str)
          res_row <- tryCatch({
            fit <- suppressMessages(lmer(form, data = df))
            coef_row <- summary(fit)$coefficients
            if (nrow(coef_row) < 2 || any(is.na(coef_row[2, ]))) {
              data.frame(Gene = target_gene, snp = snp_id, eQTL_beta = NA, se = NA, eQTL_t = NA, pvalue = 1)
            } else {
              data.frame(Gene = target_gene,
                         snp = snp_id,
                         eQTL_beta = as.numeric(coef_row[2, 1]),
                         se = as.numeric(coef_row[2, 2]),
                         eQTL_t = as.numeric(coef_row[2, 3]),
                         pvalue = as.numeric(coef_row[2, 5]),
                         stringsAsFactors = FALSE)
            }
          }, error = function(e) {
            data.frame(Gene = target_gene, snp = snp_id, eQTL_beta = NA, se = NA, eQTL_t = NA, pvalue = 1)
          })
          br_results[[ii]] <- res_row
        } # end loop over candidate SNPs
        br_df <- rbindlist(br_results, use.names = TRUE, fill = TRUE)
        br_df[, pvalue := as.numeric(pvalue)]
        if (any(br_df$pvalue <= threshold, na.rm = TRUE)) {
          backward_list[[j]] <- br_df[which.min(pvalue), .(snp, Gene, eQTL_beta, eQTL_t, pvalue)]
        } else {
          backward_list[[j]] <- data.frame(snp = NA, Gene = target_gene, eQTL_beta = NA, eQTL_t = NA, pvalue = NA, stringsAsFactors = FALSE)
        }
      } # end j loop
      sig.snps.backward <- rbindlist(backward_list, use.names = TRUE, fill = TRUE)
    } # end else multiple snps
    
    write.table(sig.snps.backward, file = file.path(output_dir, paste0("backward_", target_gene, ".txt")), 
                quote = FALSE, row.names = FALSE, sep = "\t")
  } else {
    message("Not significant for gene: ", target_gene)
  }
} # end per-gene loop
