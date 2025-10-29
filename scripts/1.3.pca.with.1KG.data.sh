#!/usr/bin/env bash
# pca.sh — Perform PCA on merged genotype data with 1000 Genomes data.
# Example:
#   ./pca.sh -m Merge.1KG -i input.genotypes -r /path/to/common.snp.R -t 12

set -euo pipefail
IFS=$'\n\t'

THREADS=12
MERGE_PREFIX="Merge.1KG"
INPUT_PREFIX="input.genotypes"
RSCRIPT_PATH=""
SNPS_FILE="list.snps"
PLINK=${PLINK:-$(command -v plink || true)}
RSCRIPT=${RSCRIPT:-$(command -v Rscript || true)}

while getopts "m:i:r:s:t:" opt; do
  case $opt in
    m) MERGE_PREFIX=$OPTARG ;;
    i) INPUT_PREFIX=$OPTARG ;;
    r) RSCRIPT_PATH=$OPTARG ;;
    s) SNPS_FILE=$OPTARG ;;
    t) THREADS=$OPTARG ;;
  esac
done

echo "Started at: $(date)"

"$RSCRIPT" "$RSCRIPT_PATH" "$MERGE_PREFIX.bim" "${INPUT_PREFIX}.bim" "$SNPS_FILE"

"$PLINK" --threads "$THREADS" --bfile "$MERGE_PREFIX" --extract "$SNPS_FILE" --make-bed --out "${MERGE_PREFIX}.shared"
"$PLINK" --threads "$THREADS" --bfile "$INPUT_PREFIX" --extract "$SNPS_FILE" --make-bed --out "${INPUT_PREFIX}.shared"

"$PLINK" --threads "$THREADS" --bfile "${MERGE_PREFIX}.shared" \
  --bmerge "${INPUT_PREFIX}.shared.bed" "${INPUT_PREFIX}.shared.bim" "${INPUT_PREFIX}.shared.fam" \
  --make-bed --out "${MERGE_PREFIX}.input"

"$PLINK" --threads "$THREADS" --bfile "${MERGE_PREFIX}.input" --maf 0.001 --indep 50 5 1.5 --out "${MERGE_PREFIX}.input"
"$PLINK" --threads "$THREADS" --bfile "${MERGE_PREFIX}.input" --extract "${MERGE_PREFIX}.input.prune.in" \
  --make-bed --out "${MERGE_PREFIX}.input.pruned"

"$PLINK" --threads "$THREADS" --bfile "${MERGE_PREFIX}.input.pruned" --pca

echo "Ended at: $(date)"
