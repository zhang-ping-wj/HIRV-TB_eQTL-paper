#!/usr/bin/env Rscript
# Coloc_run_arrayR-eQTL.r — coloc analysis for eQTL catalogue vs TB summary stats
# Usage:
#   Rscript Coloc_run_arrayR-eQTL.r <dataset_id>
# Requirements:
#   - R packages: data.table, dplyr, readr, seqminer, coloc, GenomicRanges, ggplot2
#   - tabix index files accessible at paths listed in tabix_paths.csv


suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(readr)
  library(seqminer)
  library(coloc)
  library(GenomicRanges)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript coloc_loop.R <dataset_id>\nExample: Rscript coloc_loop.R study123_groupA")
}
dataset_id <- args[1]

# ----------------------- user-editable paths -----------------------
tabix_paths_file <- "tabix_paths.csv"   # CSV with columns: study, qtl_group, ftp_path, quant_method
tb_summary_file  <- "TB_LMM_results_rsid_sorted_pairs.txt.gz" # local tabix-indexed TB sumstats (hg19)
gene_regions_file<- "eigenMT.sig.hg38.gene.id_lmm_qvalue0.05.csv" # gene regions with hg38/hg19 coords & ids
# output directory
out_dir <- "coloc_results"
dir.create(out_dir, showWarnings = FALSE)

# ----------------------- helper functions -------------------------
import_eQTLCatalogue <- function(ftp_path, region, selected_gene_id, column_names, verbose = FALSE) {
  # Fetch summary statistics (hg38) for selected_gene_id from eQTL Catalogue tabix file
  if (verbose) message("Fetching: ", ftp_path, " region: ", region)
  fetch_table <- tryCatch(
    seqminer::tabix.read.table(tabixFile = ftp_path, tabixRange = region, stringsAsFactors = FALSE),
    error = function(e) return(NULL)
  )
  if (is.null(fetch_table) || nrow(fetch_table) == 0) return(NULL)
  fetch_table <- as_tibble(fetch_table)
  if (!missing(column_names) && length(column_names) == ncol(fetch_table)) {
    colnames(fetch_table) <- column_names
  }
  summary_stats <- fetch_table %>%
    filter(gene_id == selected_gene_id) %>%
    mutate(id = paste(chromosome, position, sep = ":")) %>%
    group_by(id) %>%
    mutate(row_count = n()) %>%
    ungroup() %>%
    filter(row_count == 1)  # remove multi-allelic entries (duplicated id)
  return(as.data.table(summary_stats))
}

import_TBregion_hg19 <- function(region, selected_gene_name, tb_file) {
  # Fetch TB sumstats (hg19) for selected_gene_name
  fetch_table <- tryCatch(
    seqminer::tabix.read.table(tabixFile = tb_file, tabixRange = region, stringsAsFactors = FALSE),
    error = function(e) return(NULL)
  )
  if (is.null(fetch_table) || nrow(fetch_table) == 0) return(NULL)
  fetch_table <- as_tibble(fetch_table)
  # adjust these column names to match your TB file format
  colnames(fetch_table) <- c('chr', 'position_hg19', 'rsid', 'ref', 'alt',
                             'variant_id', 'gene_name', 'Estimate', 'se', 't.value', 'p', 'adj_P','sig',
                             'maf','pair')
  summary_stats <- fetch_table %>% filter(gene_name == selected_gene_name)
  return(as.data.table(summary_stats))
}

run_coloc <- function(TBeqtl, eqtl_sumstats) {
  # TBeqtl: data.table with columns beta, se, varbeta, snp/rsid, MAF, N
  # eqtl_sumstats: data.table with columns beta, se, maf, an (allele number) or N, rsid
  if (nrow(TBeqtl) == 0 || nrow(eqtl_sumstats) == 0) return(NULL)
  # prepare dataset1 (TB)
  dataset1 <- list(
    beta = TBeqtl$beta,
    varbeta = TBeqtl$varbeta,
    type = "quant",
    snp = TBeqtl$snp,
    MAF = TBeqtl$MAF,
    N = TBeqtl$N
  )
  # prepare dataset2 (eQTL catalogue)
  dataset2 <- list(
    beta = eqtl_sumstats$beta,
    varbeta = (eqtl_sumstats$se)^2,
    type = "quant",
    snp = eqtl_sumstats$rsid,
    MAF = eqtl_sumstats$maf,
    N = if ("an" %in% colnames(eqtl_sumstats)) eqtl_sumstats$an[1] / 2 else eqtl_sumstats$N[1]
  )
  coloc_res <- tryCatch(
    coloc::coloc.abf(dataset1 = dataset1, dataset2 = dataset2, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5),
    error = function(e) return(NULL) # default coloc ABF priors
  )
  if (is.null(coloc_res)) return(NULL)
  res_formatted <- as_tibble(t(as.data.frame(coloc_res$summary)))
  return(as.data.table(res_formatted))
}

# ----------------------- main loop setup --------------------------
tabix_paths <- fread(tabix_paths_file, stringsAsFactors = FALSE)
tabix_paths <- tabix_paths %>%
  mutate(datasets = paste0(study, "_", qtl_group)) %>%
  filter(quant_method %in% c("ge", "microarray")) %>%
  as_tibble()

# find selected dataset row
platelet_df <- tabix_paths %>% filter(datasets == dataset_id)
if (nrow(platelet_df) == 0) stop("Dataset id not found in tabix_paths.csv: ", dataset_id)
platelet_df <- platelet_df[1, ]  # take first match

# get column names from the remote file header (safe-guarded)
colnames_remote <- tryCatch({
  read_tsv(platelet_df$ftp_path, n_max = 1) %>% colnames()
}, error = function(e) {
  warning("Could not read header from remote file; proceeding without setting column names.")
  NULL
})

# load gene regions to iterate
gene_regions <- fread(gene_regions_file, stringsAsFactors = FALSE)
genes <- unique(gene_regions$gene_name)

results_list <- list()

for (g in genes) {
  region_hg38 <- gene_regions$hg38_ID[gene_regions$gene_name == g][1]
  region_hg19 <- gene_regions$hg19_ID[gene_regions$gene_name == g][1]
  gene_id     <- gene_regions$gene_id[gene_regions$gene_name == g][1]
  if (is.na(region_hg38) || is.na(region_hg19) || is.na(gene_id)) next
  
  # quick check if remote region has any data (seqminer returns NULL or empty)
  try_fetch <- try(seqminer::tabix.read.table(tabixFile = platelet_df$ftp_path, tabixRange = region_hg38, stringsAsFactors = FALSE), silent = TRUE)
  if (inherits(try_fetch, "try-error") || length(try_fetch) == 0) next
  
  summary_stats <- import_eQTLCatalogue(platelet_df$ftp_path, region_hg38, gene_id, colnames_remote, verbose = FALSE)
  if (is.null(summary_stats) || nrow(summary_stats) == 0) next
  
  # TB sumstats (hg19)
  TBeqtl <- import_TBregion_hg19(region_hg19, selected_gene_name = g, tb_file = tb_summary_file)
  if (is.null(TBeqtl) || nrow(TBeqtl) == 0) next
  
  # harmonise rsid strings and subset overlapping variants
  TBeqtl$rsid <- gsub("\\s+|\\r", "", TBeqtl$rsid)
  summary_stats$rsid <- gsub("\\s+|\\r", "", summary_stats$rsid)
  summary_stats2 <- summary_stats[summary_stats$rsid %in% TBeqtl$rsid, ]
  if (nrow(summary_stats2) == 0) next
  
  # prepare TBeqtl fields expected by run_coloc
  TBeqtl <- copy(TBeqtl)
  TBeqtl[, `:=`(
    beta = Estimate,
    se = se,
    varbeta = se^2,
    snp = rsid,
    position = position_hg19,
    type = "quant",
    N = 267,
    MAF = maf
  )]
  TBeqtl <- TBeqtl[!duplicated(rsid)]
  
  # run coloc
  coloc_res <- run_coloc(TBeqtl, summary_stats2)
  if (is.null(coloc_res) || nrow(coloc_res) == 0) next
  coloc_res[, datasets := platelet_df$datasets]
  
  # add TB top SNP if available
  top_variant_id <- gene_regions$variant_id[gene_regions$gene_name == g][1]
  if (!is.na(top_variant_id) && top_variant_id %in% TBeqtl$variant_id) {
    coloc_res[, tb_top_rsid := TBeqtl[variant_id == top_variant_id, rsid][1]]
  } else {
    coloc_res[, tb_top_rsid := "Not.in.the.eQTL.catalog"]
  }
  
  coloc_res[, c("skin.nsnps", "TB.nsnps", "gene_name") := .(length(unique(summary_stats$rsid)), length(unique(TBeqtl$rsid)), g)]
  results_list[[g]] <- coloc_res
}

# combine and write results
if (length(results_list) > 0) {
  results_dt <- rbindlist(results_list, use.names = TRUE, fill = TRUE)
  out_file <- file.path(out_dir, paste0("coloc.res.datasets_", platelet_df$datasets, "_LMM.txt"))
  fwrite(results_dt, out_file, sep = "\t", quote = FALSE)
  message("Wrote results to: ", out_file)
} else {
  message("No coloc results generated for dataset: ", dataset_id)
}
