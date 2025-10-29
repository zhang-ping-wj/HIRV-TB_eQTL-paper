#!/usr/bin/env bash
# featurecount.sh - generate read counts using featureCounts
# Runs featureCounts on all BAM files from HISAT2 mapping

set -euo pipefail
IFS=$'\n\t'

# ---- basic settings ----
THREADS=8
GTF="/path/to/gencode.v31lift37.annotation.gtf.gz"   # ← edit this path
BAM_DIR="Mapping/mapped.hisat2"
OUT_DIR="featureCounts"
OUT_FILE="${OUT_DIR}/featureCounts_all.txt"
FEATURECOUNTS=${FEATURECOUNTS:-$(command -v featureCounts || true)}

mkdir -p "$OUT_DIR"

echo "Started at: $(date)"

"$FEATURECOUNTS" -T "$THREADS" \
  -p -s 0 \
  -a "$GTF" -g gene_name \
  --extraAttributes gene_id \
  -o "$OUT_FILE" \
  "${BAM_DIR}"/*.bam

echo "Ended at: $(date)"
