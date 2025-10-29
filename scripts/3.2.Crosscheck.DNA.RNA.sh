#!/usr/bin/env bash
# crosscheck.sh — perform cross-checking of DNA and RNA samples using Picard's CrosscheckFingerprints.
# Usage:
#   sbatch --cpus-per-task=24 crosscheck.sh

set -euo pipefail
IFS=$'\n\t'

# ---- user-editable settings ----
THREADS=24
BAM_DIR="RG.bam"                       # directory containing *.RG_sorted.bam
MERGED_BAM="RG_sorted.merged.bam"
VCF_INPUT="TB_352.vcf.gz"             # VCF to cross-check
HAPLOTYPE_MAP="hg19_nochr.map"        # haplotype map for CrosscheckFingerprints
OUT_TXT="DNA.RNA.crosscheck_352.txt"
OUT_MATRIX="DNA.RNA.crosscheck_Matrix_352.txt"

# tools (override via env if needed)
SAMTOOLS=${SAMTOOLS:-$(command -v samtools || true)}
JAVA=${JAVA:-$(command -v java || true)}
PICARD_JAR=${PICARD_JAR:-picard.jar}   # set to full path of picard.jar if required

# merge BAMs
mkdir -p "$(dirname "$MERGED_BAM")"
"$SAMTOOLS" merge -@ "$THREADS" "$MERGED_BAM" "${BAM_DIR}"/*.RG_sorted.bam

# optional: index merged BAM (uncomment if desired)
# "$SAMTOOLS" index -@ "$THREADS" "$MERGED_BAM"

# run CrosscheckFingerprints
"$JAVA" -Xmx80g -jar "$PICARD_JAR" CrosscheckFingerprints \
  INPUT="$MERGED_BAM" \
  INPUT="$VCF_INPUT" \
  HAPLOTYPE_MAP="$HAPLOTYPE_MAP" \
  NUM_THREADS="$THREADS" \
  OUTPUT="$OUT_TXT" \
  MATRIX_OUTPUT="$OUT_MATRIX"
