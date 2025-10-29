#!/usr/bin/env bash
# for.PCA.sh - Filter imputed VCFs and convert to PLINK format for PCA
# Usage: ./for.PCA.sh -i "chr*.dose.vcf.gz" -r /path/to/human_g1k_v37.fasta -k /path/to/1KG_dir -o input.genotypes -t 12
#
# Steps:
# 1) index per-chromosome VCFs
# 2) concat into Concat.MIS_HRC.vcf.gz
# 3) count SNPs with INFO/R2>=0.1
# 4) filter INFO/R2<0.8
# 5) normalise / annotate and convert to PLINK bed/bim/fam

set -euo pipefail
IFS=$'\n\t'

# defaults
THREADS=12
PATTERN="chr*.dose.vcf.gz"
REF=""
ONEKG_DIR=""
OUT_PREFIX="input.genotypes"
PLINK_PATH=${PLINK_PATH:-$(command -v plink || true)}
BCFTOOLS=${BCFTOOLS:-$(command -v bcftools || true)}

usage() {
  cat <<EOF
Usage: $0 -i <vcf_glob> -r <reference_fasta> [-k <1KG_dir>] [-o <out_prefix>] [-t <threads>]
  -i VCF glob pattern (default: ${PATTERN})
  -r reference FASTA (required for bcftools norm)
  -k optional directory with Merge.1KG.{bed,bim,fam} and common.snp.R
  -o output prefix for PLINK files (default: ${OUT_PREFIX})
  -t threads (default: ${THREADS})
EOF
  exit 1
}

while getopts "i:r:k:o:t:h" opt; do
  case $opt in
    i) PATTERN=$OPTARG ;;
    r) REF=$OPTARG ;;
    k) ONEKG_DIR=${OPTARG%/} ;;
    o) OUT_PREFIX=$OPTARG ;;
    t) THREADS=$OPTARG ;;
    h) usage ;;
    *) usage ;;
  esac
done

[[ -n "$REF" ]] || { echo "ERROR: reference FASTA (-r) is required"; usage; }

# tool checks
for t in "$BCFTOOLS" "$PLINK_PATH"; do
  [[ -x "$t" ]] || { echo "ERROR: required tool missing or not executable: $t"; exit 1; }
done

echo "Started at: $(date)"
echo "Pattern: $PATTERN"
echo "Threads: $THREADS"
echo "Reference: $REF"

# 1) index each input VCF if needed
echo "Indexing input VCFs (if needed)..."
shopt -s nullglob
FILES=( $PATTERN )
if [[ ${#FILES[@]} -eq 0 ]]; then
  echo "ERROR: no files matched pattern: $PATTERN"
  exit 1
fi

for z in "${FILES[@]}"; do
  echo "Indexing $z"
  $BCFTOOLS index --threads "$THREADS" "$z"
done

# 2) concat
OUT_CONCAT="Concat.MIS_HRC.vcf.gz"
echo "Concatenating to ${OUT_CONCAT}..."
$BCFTOOLS concat --threads "$THREADS" "${FILES[@]}" -Oz -o "$OUT_CONCAT"
$BCFTOOLS index --threads "$THREADS" "$OUT_CONCAT"

# 3) count SNPs with INFO/R2>=0.1 (safer than zcat | grep)
echo "### Check SNP number with INFO/R2>=0.1"
$BCFTOOLS view -i 'INFO/R2>=0.1' -H "$OUT_CONCAT" | wc -l

# 4) filter out INFO/R2 < 0.8
OUT_FILT="${OUT_CONCAT%.vcf.gz}.filtered.vcf.gz"
echo "Filtering out INFO/R2<0.8 -> ${OUT_FILT}"
$BCFTOOLS filter --threads "$THREADS" -e 'INFO/R2<0.8' "$OUT_CONCAT" -Oz -o "$OUT_FILT"
$BCFTOOLS index --threads "$THREADS" "$OUT_FILT"

echo "### Check SNP number with R2>=0.8"
$BCFTOOLS view -H "$OUT_FILT" | wc -l

# 5) normalize / annotate ID -> CHROM:POS:REF:ALT and convert to PLINK
echo "Normalising, annotating and converting to PLINK ($OUT_PREFIX)..."
# The pipeline writes to stdout BCF which plink reads from /dev/stdin
$BCFTOOLS norm --threads "$THREADS" -m -any "$OUT_FILT" \
  | $BCFTOOLS norm --threads "$THREADS" -f "$REF" -Ou \
  | $BCFTOOLS annotate --threads "$THREADS" -Ob -x ID -I +'%CHROM:%POS:%REF:%ALT' \
  | "$PLINK_PATH" --bcf /dev/stdin \
       --keep-allele-order \
       --vcf-idspace-to _ \
       --const-fid \
       --allow-extra-chr 0 \
       --split-x b37 no-fail \
       --make-bed \
       --out "$OUT_PREFIX"

# 6) collect outputs into PCA.with.1KG.data and  copy 1KG files
mkdir -p PCA.with.1KG.data
mv ${OUT_PREFIX}".{bed,bim,fam} PCA.with.1KG.data/ 

if [[ -n "$ONEKG_DIR" ]]; then
  echo "Copying 1KG reference files from ${ONEKG_DIR} ..."
  for f in Merge.1KG.bed Merge.1KG.bim Merge.1KG.fam common.snp.R; do
    src="${ONEKG_DIR}/${f}"
    if [[ -f "$src" ]]; then
      cp "$src" PCA.with.1KG.data/
    else
      echo "WARN: missing ${f} in ${ONEKG_DIR}"
    fi
  done
fi

echo "Ended at: $(date)"
