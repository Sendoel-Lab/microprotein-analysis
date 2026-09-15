#!/usr/bin/env bash
# Adapter/quality trimming of paired-end RNA-seq reads with fastp.
#
# Usage: rnaseq/01_trim.sh <samples.tsv>
#   samples.tsv columns (tab-separated, header): sample  fastq_dir
#   fastq_dir holds one or more lanes named <lane>_1.fq.gz / <lane>_2.fq.gz.
# Output: ${MP_DATA}/rnaseq/trimmed/<sample>/<lane>_{1,2}.trim.fq.gz
set -euo pipefail
: "${MP_DATA:?source config.env first}"

SAMPLES="$1"
OUT="$MP_DATA/rnaseq/trimmed"
THREADS="${MP_THREADS:-8}"
MIN_LEN=36
mkdir -p "$OUT/qc"

tail -n +2 "$SAMPLES" | while IFS=$'\t' read -r sample fastq_dir; do
    [[ -z "$sample" ]] && continue
    mkdir -p "$OUT/$sample"
    for r1 in "$fastq_dir"/*_1.fq.gz; do
        r2="${r1%_1.fq.gz}_2.fq.gz"
        lane=$(basename "$r1" _1.fq.gz)
        fastp \
            -i "$r1" -I "$r2" \
            -o "$OUT/$sample/${lane}_1.trim.fq.gz" \
            -O "$OUT/$sample/${lane}_2.trim.fq.gz" \
            --detect_adapter_for_pe \
            -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
            -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
            --trim_poly_g \
            -l "$MIN_LEN" \
            -h "$OUT/qc/${lane}.fastp.html" -j "$OUT/qc/${lane}.fastp.json" \
            --thread "$THREADS"
    done
done
