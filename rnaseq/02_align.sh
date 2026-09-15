#!/usr/bin/env bash
# STAR alignment of trimmed paired-end RNA-seq reads; lanes of a sample are
# passed to STAR as comma-separated lists.
#
# Usage: rnaseq/02_align.sh <sample> [<sample> ...]
# Input:  ${MP_DATA}/rnaseq/trimmed/<sample>/*_{1,2}.trim.fq.gz
# Output: ${MP_DATA}/rnaseq/star/<sample>/<sample>.Aligned.sortedByCoord.out.bam (+ .bai)
#         ${MP_DATA}/rnaseq/star/<sample>/<sample>.Aligned.toTranscriptome.out.bam
set -euo pipefail
: "${MP_DATA:?source config.env first}"
: "${MP_REFS:?source config.env first}"

GENOME_DIR="$MP_REFS/star_index"
THREADS="${MP_THREADS:-16}"

for sample in "$@"; do
    in_dir="$MP_DATA/rnaseq/trimmed/$sample"
    out_dir="$MP_DATA/rnaseq/star/$sample"
    mkdir -p "$out_dir"
    r1_list=$(ls "$in_dir"/*_1.trim.fq.gz | sort | paste -sd, -)
    r2_list=$(ls "$in_dir"/*_2.trim.fq.gz | sort | paste -sd, -)
    (
        cd "$out_dir"
        STAR \
            --runThreadN "$THREADS" \
            --genomeDir "$GENOME_DIR" \
            --readFilesIn "$r1_list" "$r2_list" \
            --readFilesCommand zcat \
            --outFileNamePrefix "${sample}." \
            --twopassMode Basic \
            --twopass1readsN -1 \
            --alignIntronMin 20 \
            --outFilterMultimapNmax 20 \
            --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
            --outFilterMismatchNmax 1 \
            --outFilterMismatchNoverReadLmax 0.04 \
            --outFilterType BySJout \
            --outSAMattributes AS NH HI nM MD \
            --outSAMtype BAM SortedByCoordinate \
            --quantMode TranscriptomeSAM \
            --limitBAMsortRAM 60000000000
    )
    samtools index -@ 4 "$out_dir/${sample}.Aligned.sortedByCoord.out.bam"
done
