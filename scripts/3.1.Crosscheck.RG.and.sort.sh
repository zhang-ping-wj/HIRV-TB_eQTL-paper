#!/usr/bin/env bash
# dedup_rg_reorder.sh — remove duplicates, add read groups snf and reorder BAM files
# Example (SLURM): sbatch --array=1-424 --cpus-per-task=12 dedup_rg_reorder.sh
# Edit the variables below as needed.

set -euo pipefail
IFS=$'\n\t'

THREADS=12
FOLDER="."                                # base folder containing Mapping/mapped.hisat2 and sample.key_424.txt
KEY_FILE="${FOLDER}/sample.key_424.txt"
OUT_DEDUP_DIR="Deduplicated.bam"
OUT_RG_DIR="RG.bam"
TMPDIR="${PWD}/tmp"
PICARD_JAR="${PICARD_JAR:-picard.jar}"   # set PICARD_JAR to full path if needed
JAVA=${JAVA:-java}
SAMTOOLS=${SAMTOOLS:-samtools}

mkdir -p "$OUT_DEDUP_DIR" "$OUT_RG_DIR" "$TMPDIR"

LINE=${SLURM_ARRAY_TASK_ID:-1}
SAMPLE=$(awk "NR==$LINE {print \$1}" "$KEY_FILE")

# Mark duplicates
"$JAVA" -Xmx8g -jar "$PICARD_JAR" MarkDuplicates \
  INPUT="${FOLDER}/Mapping/mapped.hisat2/${SAMPLE}.hisat2.bam" \
  OUTPUT="${OUT_DEDUP_DIR}/${SAMPLE}.dedup.bam" \
  REMOVE_DUPLICATES=true \
  METRICS_FILE="${OUT_DEDUP_DIR}/${SAMPLE}.metrics.txt"

"$SAMTOOLS" index -@ "$THREADS" "${OUT_DEDUP_DIR}/${SAMPLE}.dedup.bam"

# Keep uniquely-mapped reads (mapq >= 60)
"$SAMTOOLS" view -@ "$THREADS" -q 60 -b "${OUT_DEDUP_DIR}/${SAMPLE}.dedup.bam" > "${OUT_DEDUP_DIR}/${SAMPLE}.dedup.uniqmap.bam"
"$SAMTOOLS" index -@ "$THREADS" "${OUT_DEDUP_DIR}/${SAMPLE}.dedup.uniqmap.bam"

# Add / replace read groups
"$JAVA" -Xmx8g -jar "$PICARD_JAR" AddOrReplaceReadGroups \
  INPUT="${OUT_DEDUP_DIR}/${SAMPLE}.dedup.uniqmap.bam" \
  OUTPUT="${OUT_RG_DIR}/${SAMPLE}.RG.bam" \
  RGID="$SAMPLE" RGLB=RNAseq RGPL=illumina RGPU="$SAMPLE" RGSM="$SAMPLE"

"$SAMTOOLS" index -@ "$THREADS" "${OUT_RG_DIR}/${SAMPLE}.RG.bam"

# Reorder to match reference sequence dictionary
# Set SEQ_DICT to the fasta sequence dictionary 
SEQ_DICT="${SEQ_DICT:-}"   # e.g. "GRCh37.dict"
if [[ -n "$SEQ_DICT" ]]; then
  "$JAVA" -Xmx80g -Djava.io.tmpdir="$TMPDIR" -jar "$PICARD_JAR" ReorderSam \
    INPUT="${OUT_RG_DIR}/${SAMPLE}.RG.bam" \
    OUTPUT="${OUT_RG_DIR}/${SAMPLE}.RG_sorted.bam" \
    SEQUENCE_DICTIONARY="$SEQ_DICT" \
    TMP_DIR="$TMPDIR"
fi
