#!/usr/bin/env bash
# hisat2_mapping.sh - RNA-seq read mapping with HISAT2   
# Runs HISAT2 + samtools sort/index for each FASTQ pair listed in sample.key.txt

set -euo pipefail
IFS=$'\n\t'

# basic settings
THREADS=8
KEY_FILE="sample.key.txt"
INDEX="/path/to/genome_index"        # ← change this to your HISAT2 index base name
ADAPTER_DIR="Mapping/AdapterTrimmed"
OUT_DIR="Mapping/mapped.hisat2"

# tools 
HISAT2=${HISAT2:-$(command -v hisat2 || true)}
SAMTOOLS=${SAMTOOLS:-$(command -v samtools || true)}

mkdir -p "$OUT_DIR"

# read sample name from array index (or default to first)
LINE=${SLURM_ARRAY_TASK_ID:-1}
SAMPLE=$(awk "NR==$LINE {print \$1}" "$KEY_FILE")

FASTQ1="${ADAPTER_DIR}/${SAMPLE}_R1_001_val_1.fq.gz"
FASTQ2="${ADAPTER_DIR}/${SAMPLE}_R2_001_val_2.fq.gz"

"$HISAT2" --threads "$THREADS" -x "$INDEX" -1 "$FASTQ1" -2 "$FASTQ2" \
  | "$SAMTOOLS" sort -@ "$THREADS" -o "${OUT_DIR}/${SAMPLE}.hisat2.bam"

"$SAMTOOLS" index -@ "$THREADS" "${OUT_DIR}/${SAMPLE}.hisat2.bam"
