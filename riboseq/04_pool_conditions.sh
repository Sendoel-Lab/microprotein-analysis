#!/usr/bin/env bash
# Merge per-sample transcriptome BAMs by condition and run riboWaltz QC on each pool.
#
# Usage: riboseq/04_pool_conditions.sh <pools.tsv> <out_dir>
#   pools.tsv columns (tab-separated, header): pool  bam
#   One row per member BAM; ${VAR} in bam paths is expanded from the environment.
# Output: <out_dir>/<pool>/bam/<pool>.bam and riboWaltz outputs in <out_dir>/<pool>/
set -euo pipefail

POOLS="$1"; OUT="$2"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
THREADS="${MP_THREADS:-8}"
mkdir -p "$OUT"

for pool in $(tail -n +2 "$POOLS" | cut -f1 | awk '!seen[$0]++'); do
    bams=()
    while IFS=$'\t' read -r p bam; do
        [[ "$p" == "$pool" ]] && bams+=("$(eval echo "$bam")")
    done < <(tail -n +2 "$POOLS")
    mkdir -p "$OUT/$pool/bam"
    samtools merge -f -@ "$THREADS" "$OUT/$pool/bam/$pool.bam" "${bams[@]}"
    samtools index "$OUT/$pool/bam/$pool.bam"
    "$HERE/03_ribowaltz_qc.sh" "$pool" "$OUT/$pool/bam" "$pool" "$OUT/$pool"
done

Rscript "$HERE/05_pooled_inframe_plot.R" "$OUT" $(tail -n +2 "$POOLS" | cut -f1 | awk '!seen[$0]++')
