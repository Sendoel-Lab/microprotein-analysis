#!/usr/bin/env bash
# riboWaltz QC of one UMI-deduplicated transcriptome BAM.
#
# Usage: riboseq/03_ribowaltz_qc.sh <sample> <bam_dir> <bam_prefix> <out_dir>
#   <bam_dir>/<bam_prefix>.bam is the transcriptome BAM (e.g. prefix <sample>.tx.dedup).
# Annotation: ${MP_REFS}/ribowaltz_annotation.csv (columns transcript, l_tr, l_utr5, l_cds, l_utr3)
set -euo pipefail
: "${MP_REFS:?source config.env first}"

SAMPLE="$1"; BAM_DIR="$2"; PREFIX="$3"; OUT="$4"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

Rscript "$HERE/ribowaltz_qc.R" \
    --bam_dir "$BAM_DIR" \
    --bam_prefix "$PREFIX" \
    --sample "$SAMPLE" \
    --annotation "$MP_REFS/ribowaltz_annotation.csv" \
    --output_dir "$OUT" \
    --cds_only
