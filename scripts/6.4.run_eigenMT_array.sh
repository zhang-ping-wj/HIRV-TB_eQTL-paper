#!/usr/bin/env bash
# run_eigenMT_array.sh — SLURM array script to run eigenMT per chromosome
# Usage: sbatch --array=1-22 --cpus-per-task=4 run_eigenMT_array.sh

set -euo pipefail
IFS=$'\n\t'

# ---- user settings ----
THREADS=4
CHROM=${SLURM_ARRAY_TASK_ID:-1}
OUT_PREFIX="QTL.chr${CHROM}"
INPUT_DIR="../input"

# ---- 1. R pre-processing step ----
module load R/4.2.1-foss-2022a

Rscript 2.1.eigenMT-input.r \
  --CHR "chr${CHROM}" \
  --output "${OUT_PREFIX}.txt" \
  --verbose

# ---- 2. Run eigenMT ----
module purge
module load Anaconda3/2023.07-2
eval "$(conda shell.bash hook)"
conda activate eigenMT

python eigenMT_python3.py \
  --CHROM "chr${CHROM}" \
  --QTL "${OUT_PREFIX}.txt" \
  --GEN "${INPUT_DIR}/CHR.${CHROM}.SNP_rsid_ALL.267.TB.txt" \
  --GENPOS "${INPUT_DIR}/CHR.${CHROM}.snpsloc_rsid_ALL.267.TB.txt" \
  --PHEPOS "${INPUT_DIR}/geneloc_267.txt" \
  --window 200 \
  --OUT "Results.${OUT_PREFIX}.txt"

echo "Finished chromosome ${CHROM} at $(date)"
