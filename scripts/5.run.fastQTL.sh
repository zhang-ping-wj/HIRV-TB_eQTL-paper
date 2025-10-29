#!/usr/bin/env bash
# run_fastqtl.sh — run FastQTL for permutations and nominal passes
# Example (SLURM): sbatch --array=1-100 --cpus-per-task=4 run_fastqtl.sh

set -euo pipefail
IFS=$'\n\t'

# ---- user settings ----
THREADS=4
CHUNKS=100

# Input files
VCF_D2="temp.d2.vcf.gz"
VCF_D7="temp.d7.vcf.gz"
VCF_ACTIVE="temp.active.vcf.gz"

BED_D2="PhenotypesFASTqtl_sorted.d2.bed.gz"
BED_D7="PhenotypesFASTqtl_sorted.d7.bed.gz"
BED_ACTIVE="PhenotypesFASTqtl_sorted.active.bed.gz"

COV_D2="Covariates.d2.31.PCs.txt"
COV_D7="Covariates.d7.25.PCs.txt"
COV_ACTIVE="Covariates.active.11.PCs.txt"

# Output prefixes
OUT_DIR="fastqtl_results"
mkdir -p "$OUT_DIR"

# FastQTL executable (update if not in PATH)
FASTQTL=${FASTQTL:-$(command -v fastQTL || true)}

echo "------------------------------------------------"
echo "$(date): Running task ${SLURM_ARRAY_TASK_ID:-1} on $(hostname)"
echo "------------------------------------------------"

# ----------------------------
# 1. Nominal pass
# ----------------------------
"$FASTQTL" -V "$VCF_D2" -B "$BED_D2" -C "$COV_D2" \
  -O "${OUT_DIR}/res_d2_${SLURM_ARRAY_TASK_ID:-1}_${CHUNKS}.txt" \
  --chunk "${SLURM_ARRAY_TASK_ID:-1}" "$CHUNKS"

"$FASTQTL" -V "$VCF_D7" -B "$BED_D7" -C "$COV_D7" \
  -O "${OUT_DIR}/res_d7_${SLURM_ARRAY_TASK_ID:-1}_${CHUNKS}.txt" \
  --chunk "${SLURM_ARRAY_TASK_ID:-1}" "$CHUNKS"

"$FASTQTL" -V "$VCF_ACTIVE" -B "$BED_ACTIVE" -C "$COV_ACTIVE" \
  -O "${OUT_DIR}/res_active_${SLURM_ARRAY_TASK_ID:-1}_${CHUNKS}.txt" \
  --chunk "${SLURM_ARRAY_TASK_ID:-1}" "$CHUNKS"

# ----------------------------
# 2. Permutation pass
# ----------------------------
"$FASTQTL" -V "$VCF_D2" -B "$BED_D2" -C "$COV_D2" \
  -O "${OUT_DIR}/permute_d2_${SLURM_ARRAY_TASK_ID:-1}_${CHUNKS}.txt" \
  --chunk "${SLURM_ARRAY_TASK_ID:-1}" "$CHUNKS" --permute 1000 --seed 123456

"$FASTQTL" -V "$VCF_D7" -B "$BED_D7" -C "$COV_D7" \
  -O "${OUT_DIR}/permute_d7_${SLURM_ARRAY_TASK_ID:-1}_${CHUNKS}.txt" \
  --chunk "${SLURM_ARRAY_TASK_ID:-1}" "$CHUNKS" --permute 1000 --seed 123456

"$FASTQTL" -V "$VCF_ACTIVE" -B "$BED_ACTIVE" -C "$COV_ACTIVE" \
  -O "${OUT_DIR}/permute_active_${SLURM_ARRAY_TASK_ID:-1}_${CHUNKS}.txt" \
  --chunk "${SLURM_ARRAY_TASK_ID:-1}" "$CHUNKS" --permute 1000 --seed 123456

echo "Completed at: $(date)"
