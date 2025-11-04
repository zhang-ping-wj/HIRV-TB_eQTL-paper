# HIRV-TB_eQTL-paper

## Overview

This repository contains scripts for the paper:

Author list xxxxx

“**Common genetic determinants of longitudinal human *in vivo* immune response to Mycobacterium tuberculosis.**”

xxxxxxx


## Pipeline

| Step    | Script                                            | Description                                                                         |
| ------- | ------------------------------------------------- | ----------------------------------------------------------------------------------- |
| **1.1** | `vcf_for_imputation.sh`                           | Prepare VCF for imputation by filtering and splitting per chromosome                |
| **1.2** | `for.PCA.sh`                                      | Filter imputed VCFs and convert to PLINK format for PCA                             |
| **1.3** | `pca.sh`                                          | Perform PCA on merged genotype data with 1000 Genomes reference data                |
| **2.1** | `hisat2_mapping.sh`                               | RNA-seq read mapping using HISAT2                                                   |
| **2.2** | `featurecount.sh`                                 | Generate read counts using featureCounts                                            |
| **3.1** | `dedup_rg_reorder.sh`                             | Remove duplicates, add read groups, and reorder BAM files                           |
| **3.2** | `crosscheck.sh`                                   | Perform cross-checking of DNA and RNA samples using Picard’s CrosscheckFingerprints |
| **4.1** | `MatrixeQTL_cis.R`                                | Run MatrixEQTL using different expression PCs                                       |
| **4.2** | `trans-day7.PC25_subset.SNPs.R`                   | trans-eQTL mapping using day7 lead cis-eQTLs                                        |
| **5.1**   | `run_fastqtl.sh`                                  | Run FastQTL for permutations and nominal passes                                     |
| **6.1** | `lmm.R`                                           | LMM model for cis-eQTL mapping                                                      |
| **6.2** | `eigenMT-input.r`                                 | Prepare eigenMT input                                                               |
| **6.3** | `eigenMT_python3.py`                              | Python eigenMT script                                                               |
| **6.4** | `run_eigenMT_array.sh`                            | SLURM array script to run eigenMT per chromosome                                    |
| **6.5** | `eigenMT-qvalue-significance.threshold.r`         | Compute q-values for eigenMT results and define significance thresholds             |
| **6.6** | `conditional.analysis_lmerTest-array-ALL.pairs.r` | Stepwise forward/backward conditional mapping using lmerTest                        |
| **7.1** | `Coloc_run_arrayR-eQTL.r`                         | Colocalisation analysis for eQTL Catalogue vs TB summary statistics                 |
| **7.2** | `Coloc_SGE_array_LMM.resul_with.GWAS_opti.r`      | Colocalisation between TB eQTLs and GWAS summary statistics                         |
| **8.1** | `LD Score (ldsc) regression analysis`             | Partitioned heritability estimation                                                 |
| **9.1** | `interaction.analysis-Day7vsDay2.r`               | SNP × TB_Status (Day7 vs Day2) interaction analysis                                 |
| **9.2** | `interaction.analysis_cell.module.r`              | SNP x cell module interaction analysis                                              |
---

** directory layout:**

```
project/
├── input/
│   ├── GE_415.txt.gz
│   ├── Covariates.PCs_415.txt
│   ├── meta.data_415.txt
│   └── CHR.*.SNP_rsid_ALL.267.TB.txt
├── scripts/
│   ├── run_fastqtl.sh
│   ├── run_matrixeqtl_day7.R
│   ├── conditional_mapping_stepwise.R
│   ├── chunked_coloc_withN267.R
│   ├── interaction_chunk_quiet.R
│   └── ...
├── nontebooks/
│   ├── run_fastqtl.sh
│   ├── run_matrixeqtl_day7.R
│   ├── conditional_mapping_stepwise.R
│   ├── chunked_coloc_withN267.R
│   ├── interaction_chunk_quiet.R
│   └── ...
│   ├── run_fastqtl.sh
│   ├── run_matrixeqtl_day7.R
│   ├── conditional_mapping_stepwise.R
│   ├── chunked_coloc_withN267.R
│   ├── interaction_chunk_quiet.R
│   └── ...
└── logs/
```


---

## Environment & dependencies

### R (≥ 4.2.1)

| Package | Purpose |
|----------|----------|
| **argparse** | Command-line argument parsing for R scripts |
| **data.table** | Fast file I/O and manipulation |
| **tidyverse** | Core data wrangling (`dplyr`, `tidyr`, `ggplot2`, `readr`) |
| **reshape** | Parse variant IDs (`colsplit`) |
| **MatrixEQTL** | cis/trans-eQTL mapping |
| **lme4**, **lmerTest** | Linear mixed models and p-values |
| **foreach**, **doMC** | Parallel loops for SNP/gene scans |
| **coloc** | Bayesian colocalisation |
| **seqminer** | Tabix-based read of summary statistics |
| **GenomicRanges** | Interval operations for eQTLs |
| **ggplot2**, **circlize**, **scales** | Visualisation |

---

### Python & command-line tools

| Tool | Version (example) | Purpose |
|------|-------------------|---------|
| **FastQTL** | ≥ 2.184 | cis/trans eQTL discovery |
| **HISAT2** | ≥ 2.1.0 | RNA-seq alignment |
| **featureCounts** | ≥ 1.6.4 | Gene quantification |
| **samtools** | ≥ 1.9 | BAM processing |
| **bcftools** | ≥ 1.10 | VCF filtering / annotation |
| **htslib** (bgzip, tabix) | ≥ 1.8 | Compression / indexing |
| **Picard Tools** | ≥ 2.21.1 | BAM deduplication, read-group assignment |
| **eigenMT** | Python 3 | Conditional cis-eQTL fine-mapping |
| **Anaconda** | ≥ 2023.07 | Environment management |

---

### HPC environment

| Tool | Purpose |
|------|----------|
| **SBATCH (SLURM)** | Array & batch job submission |
| **module load / purge** | HPC environment setup |

---

## ✉️ Contact
For questions, please contact:  
**Ping Zhang**  (<ping.zhang@well.ox.ac.uk>)
