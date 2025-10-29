#!/usr/bin/env Rscript
# 7.2.Coloc_SGE_array_LMM.resul_with.GWAS_opti.r - coloc between TB eQTLs and GWASs 
# Usage: Rscript run_coloc_chunks.R <chunk_id> <gwas_tabix_file>

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(seqminer)
  library(coloc)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript chunked_coloc_withN267.R <chunk_id> <gwas_tabix_file>")
chunk_id <- as.integer(args[1])
gwas_file <- args[2]

message("Chunk: ", chunk_id, " | GWAS file: ", gwas_file)

# ---------- parameters ----------
tb_tabix_file   <- "../../LMM.resul_with.se.N.maf_rsid_sorted_pairs.txt.gz" # TB LMM sumstats (hg19) tabix-indexed
gene_file       <- "eigenMT.result_final.txt.gz"   # eigenMT results (contains SNP, pos_hg19, qval, gene, etc.)
window_bp       <- 200000
chunks_total    <- 20
output_prefix   <- "chunk.coloc.result"

# ---------- helpers ----------
safe_tabix_read <- function(tabixFile, region) {
  res <- tryCatch(seqminer::tabix.read.table(tabixFile = tabixFile, tabixRange = region, stringsAsFactors = FALSE),
                  error = function(e) NULL)
  if (is.null(res) || length(res) == 0) return(NULL)
  as_tibble(res)
}

# ---------- load eGenes and build regions ----------
eigenMT <- fread(gene_file)
eigenMT.sig <- eigenMT[eigenMT$qval <= 0.05, ]
# parse SNP column if present
if ("SNP" %in% names(eigenMT.sig)) {
  parsed <- reshape::colsplit(eigenMT.sig$SNP, ":", names = c("chr", "pos_hg19", "ref", "alt"))
  eigenMT.sig[, c("chr", "pos_hg19", "ref", "alt")] <- parsed
}
# ensure pos_hg19 numeric
eigenMT.sig$pos_hg19 <- as.integer(eigenMT.sig$pos_hg19)
eigenMT.sig$coloc.start <- pmax(0, eigenMT.sig$pos_hg19 - window_bp)
eigenMT.sig$coloc.end   <- eigenMT.sig$pos_hg19 + window_bp
eigenMT.sig$position    <- paste0(eigenMT.sig$chr, ":", eigenMT.sig$coloc.start, "-", eigenMT.sig$coloc.end)

if (nrow(eigenMT.sig) == 0) stop("No significant eigenMT entries found.")

# chunking
f <- ceiling(seq_len(nrow(eigenMT.sig)) / nrow(eigenMT.sig) * chunks_total)
chunk_list <- split(eigenMT.sig, f)
if (is.na(chunk_id) || chunk_id < 1 || chunk_id > length(chunk_list)) stop("chunk_id out of range")
selected_chunk <- as.data.frame(chunk_list[[chunk_id]])
message("Processing ", nrow(selected_chunk), " rows (", length(unique(selected_chunk$gene)), " genes) in chunk ", chunk_id)

# ---------- import functions ----------
import_TB_eQTL <- function(region, gene_name) {
  df <- safe_tabix_read(tb_tabix_file, region)
  if (is.null(df) || nrow(df) == 0) return(NULL)
  # expected TB columns (adjust if your file differs)
  colnames(df) <- c('chr', 'position_hg19', 'rsid', 'ref', 'alt', 'variant_id', 'gene_name',
                    'Estimate', 'se', 't.value', 'p', 'adj_P', 'sig', 'maf', 'pair')[1:ncol(df)]
  dt <- as.data.table(df)
  dt <- dt[gene_name == gene_name]
  if (nrow(dt) == 0) return(NULL)
  dt[, `:=`(beta = as.numeric(Estimate),
            se = as.numeric(se),
            varbeta = as.numeric(se)^2,
            snp = rsid,
            N = 267,
            MAF = as.numeric(maf))]
  dt <- dt[!duplicated(snp)]
  return(dt)
}

import_GWAS_region <- function(region, gwas_tabix) {
  df <- safe_tabix_read(gwas_tabix, region)
  if (is.null(df) || nrow(df) == 0) return(NULL)
  df <- as_tibble(df)
  # Normalize common column names (attempt a few common variants)
  names(df) <- make.names(names(df), unique = TRUE)
  # Try to detect required fields: chr, position (or position_hg19), ref/alt, beta, standard_error, p_value
  if (!("position_hg19" %in% names(df)) && ("position" %in% names(df))) df <- rename(df, position_hg19 = position)
  # Attempt to coerce standard names
  if (!("standard_error" %in% names(df)) && ("se" %in% names(df))) df <- rename(df, standard_error = se)
  if (!("p_value" %in% names(df)) && ("pval" %in% names(df))) df <- rename(df, p_value = pval)
  # convert important columns to numeric where present
  if ("standard_error" %in% names(df)) df$standard_error <- as.numeric(as.character(df$standard_error))
  if ("p_value" %in% names(df)) df$p_value <- as.numeric(as.character(df$p_value))
  if ("beta" %in% names(df)) df$beta <- as.numeric(as.character(df$beta))
  return(as.data.table(df))
}

run_coloc_pair <- function(tbe, gwas) {
  if (is.null(tbe) || is.null(gwas) || nrow(tbe) == 0 || nrow(gwas) == 0) return(NULL)
  ds1 <- list(beta = tbe$beta, varbeta = tbe$varbeta, type = "quant", snp = tbe$snp, MAF = tbe$MAF, N = 267)
  # require beta and standard_error in gwas
  if (!("beta" %in% names(gwas)) || !("standard_error" %in% names(gwas))) return(NULL)
  ds2 <- list(beta = as.numeric(gwas$beta.corrected), varbeta = (as.numeric(gwas$standard_error))^2, type = "cc", snp = gwas$id)
  res <- tryCatch(coloc::coloc.abf(dataset1 = ds1, dataset2 = ds2, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5),
                  error = function(e) NULL)
  if (is.null(res)) return(NULL)
  return(as.data.table(as_tibble(t(as.data.frame(res$summary)))))
}

# ---------- main loop ----------
results_list <- list()
for (g in unique(selected_chunk$gene)) {
  gene_rows <- selected_chunk[selected_chunk$gene == g, ]
  region <- gene_rows$position[1]
  if (is.na(region) || region == "") next
  message("Processing gene ", g, " region ", region)
  # TB eQTLs
  tbeqtl <- import_TB_eQTL(region, g)
  if (is.null(tbeqtl) || nrow(tbeqtl) == 0) { message("  no TB eQTLs -> skip"); next }
  tbeqtl[, id := paste0(chr, ":", position_hg19, ":", ref, ":", alt)]
  # GWAS region
  gwas_dt <- import_GWAS_region(region, gwas_file)
  if (is.null(gwas_dt) || nrow(gwas_dt) == 0) { message("  no GWAS rows -> skip"); next }
  # try to create id for GWAS (chr:pos:ref:alt) if possible
  if (!("ref" %in% names(gwas_dt) && "alt" %in% names(gwas_dt))) {
    # try common allele column pairs
    possible_ref <- intersect(c("Allele1","A1","REF","ref"), names(gwas_dt))
    possible_alt <- intersect(c("Allele2","A2","ALT","alt"), names(gwas_dt))
    if (length(possible_ref)) setnames(gwas_dt, possible_ref[1], "ref")
    if (length(possible_alt)) setnames(gwas_dt, possible_alt[1], "alt")
  }
  if ("chr" %in% names(gwas_dt) && "position_hg19" %in% names(gwas_dt) && "ref" %in% names(gwas_dt) && "alt" %in% names(gwas_dt)) {
    gwas_dt[, id := paste0(chr, ":", position_hg19, ":", toupper(ref), ":", toupper(alt))]
  } else if ("variant_id" %in% names(gwas_dt)) {
    gwas_dt[, id := variant_id]
  } else {
    gwas_dt[, id := paste0(chr, ":", position_hg19)]
  }
  # harmonise rsids/ids
  tbeqtl$id <- toupper(tbeqtl$id)
  gwas_dt$id <- toupper(gwas_dt$id)
  # join to keep shared ids
  shared_ids <- intersect(tbeqtl$id, gwas_dt$id)
  if (length(shared_ids) == 0) { message("  no overlapping SNPs -> skip"); next }
  tbe_sub <- tbeqtl[id %in% shared_ids]
  gwas_sub <- gwas_dt[id %in% shared_ids]
  # ensure numeric standard_error
  if ("standard_error" %in% names(gwas_sub)) gwas_sub$standard_error <- as.numeric(as.character(gwas_sub$standard_error))
  # drop zero or NA SEs
  if ("standard_error" %in% names(gwas_sub)) gwas_sub <- gwas_sub[!is.na(standard_error) & standard_error > 0]
  if (nrow(gwas_sub) == 0) { message("  no valid GWAS SEs -> skip"); next }
  # align allele direction: if GWAS alt allele equals TB alt, keep beta, else flip sign
  if ("alt" %in% names(gwas_sub) && "alt" %in% names(tbe_sub)) {
    gwas_sub <- merge(gwas_sub, tbe_sub[, .(id, alt_tb = alt)], by = "id", all.x = TRUE)
    gwas_sub$beta.corrected <- ifelse(toupper(gwas_sub$alt) == toupper(gwas_sub$alt_tb), as.numeric(gwas_sub$beta), -as.numeric(gwas_sub$beta))
    gwas_sub$beta.corrected[is.na(gwas_sub$beta.corrected)] <- as.numeric(gwas_sub$beta[is.na(gwas_sub$beta.corrected)])
  } else {
    gwas_sub$beta.corrected <- as.numeric(gwas_sub$beta)
  }
  # run coloc
  coloc_res <- run_coloc_pair(tbe_sub, gwas_sub)
  if (is.null(coloc_res)) { message("  coloc failed or returned NULL for ", g); next }
  coloc_res[, `:=`(
    gwas.nsnps = nrow(gwas_sub),
    eQTL.nsnps = nrow(tbe_sub),
    gene_name  = g,
    dataset    = basename(gwas_file)
  )]
  results_list[[g]] <- coloc_res
} # end loop over genes

# combine and write output
if (length(results_list) == 0) {
  message("No coloc results for chunk ", chunk_id)
} else {
  out_dt <- rbindlist(results_list, use.names = TRUE, fill = TRUE)
  out_file <- paste0(output_prefix, "_", chunk_id, "_", tools::file_path_sans_ext(basename(gwas_file)), ".txt")
  fwrite(out_dt, out_file, sep = "\t", quote = FALSE)
  message("Wrote coloc results: ", out_file)
}
