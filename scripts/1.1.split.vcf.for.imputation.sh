#!/usr/bin/env bash
# vcf_for_imputation.sh - Prepare VCF for imputation by filtering and splitting per chromosome
# Usage: ./vcf_for_imputation.sh -i input.bcf -o output_prefix -t THREADS
# Example: ./vcf_for_imputation.sh -i TB.bcf -o TB_for_impute -t 12

set -euo pipefail
IFS=$'\n\t'

# defaults
THREADS=4
INPUT_BCF=""
OUT_PREFIX="output"

usage() {
  cat <<EOF
Usage: $0 -i <input.bcf> -o <out_prefix> [-t <threads>]
  -i input bcf (required)
  -o output prefix (default: ${OUT_PREFIX})
  -t threads (default: ${THREADS})
EOF
  exit 1
}

while getopts "i:o:t:h" opt; do
  case ${opt} in
    i) INPUT_BCF=${OPTARG} ;;
    o) OUT_PREFIX=${OPTARG} ;;
    t) THREADS=${OPTARG} ;;
    h) usage ;;
    *) usage ;;
  esac
done

if [[ -z "${INPUT_BCF}" ]]; then
  echo "ERROR: input file required."
  usage
fi

# Tools (allow overriding via env vars)
BCFTOOLS=${BCFTOOLS:-$(command -v bcftools || true)}
BGZIP=${BGZIP:-$(command -v bgzip || true)}
TABIX=${TABIX:-$(command -v tabix || true)}

# Optional module load (if on a cluster)
if command -v module >/dev/null 2>&1; then
  # Example; remove or comment for public repo
  module load BCFtools/1.10.2-GCC-8.3.0 || true
fi

# Validate tools
for t in "$BCFTOOLS" "$BGZIP" "$TABIX"; do
  if [[ -z "$t" || ! -x "$t" ]]; then
    echo "ERROR: required tool not found or not executable: $t"
    echo "Please install or set BCFTOOLS/BGZIP/TABIX environment variables."
    exit 1
  fi
done

echo "Running pipeline on: ${INPUT_BCF} (threads=${THREADS})"

# 1) remove XY/MT chromosomes and samples failing QC if necessary 
$BCFTOOLS view --threads "$THREADS" "${INPUT_BCF}" -t ^X,Y,MT -Ov -o "${OUT_PREFIX}.forPCA.vcf"

# 2) counts
echo "Sample count:"
$BCFTOOLS query -l "${OUT_PREFIX}.forPCA.vcf" | wc -l
echo "Variant count (lines excluding header):"
$BCFTOOLS view -H "${OUT_PREFIX}.forPCA.vcf" | wc -l

# 3) MAF filter >= 0.05
$BCFTOOLS view --threads "$THREADS" -e 'MAF<0.05' "${OUT_PREFIX}.forPCA.vcf" -Ov -o "${OUT_PREFIX}.maf0.05.vcf"
$BCFTOOLS view -H "${OUT_PREFIX}.maf0.05.vcf" | wc -l

# 4) Missingness filter <= 0.05
$BCFTOOLS view --threads "$THREADS" -i 'F_MISSING<0.05' "${OUT_PREFIX}.maf0.05.vcf" -Ov -o "${OUT_PREFIX}.maf0.05_missing0.05.vcf"
$BCFTOOLS view -H "${OUT_PREFIX}.maf0.05_missing0.05.vcf" | wc -l

# 5) SNPs only
$BCFTOOLS view --threads "$THREADS" --types snps "${OUT_PREFIX}.maf0.05_missing0.05.vcf" -Ov -o "${OUT_PREFIX}.maf0.05_missing0.05_snps.vcf"
$BCFTOOLS view -H "${OUT_PREFIX}.maf0.05_missing0.05_snps.vcf" | wc -l

# 6) Prepare for imputation: compress, index, list chromosomes, split
VCF="${OUT_PREFIX}.maf0.05_missing0.05_snps.vcf"
VCFGZ="${VCF}.gz"

echo "Compressing and indexing ${VCF}..."
"$BGZIP" -c "$VCF" > "$VCFGZ"
"$TABIX" -p vcf "$VCFGZ"

echo "Listing chromosomes..."
"$TABIX" --list-chroms "$VCFGZ" > "${OUT_PREFIX}.chromosomes.txt"

echo "Splitting per chromosome (will produce CHR.<name>.vcf.gz)..."
while IFS= read -r chr; do
  echo "Processing chromosome: $chr"
  # extract header+chr region, sort, compress
  "$TABIX" -h "$VCFGZ" "$chr" | $BCFTOOLS sort --threads "$THREADS" -Oz -o "CHR.${chr}.vcf.gz"
done < "${OUT_PREFIX}.chromosomes.txt"

echo "Done. Output prefix: ${OUT_PREFIX}"
