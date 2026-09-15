#!/usr/bin/env bash
# Ribo-seq processing of one sequencing batch: in-read demultiplexing and UMI
# extraction, rRNA/tRNA depletion, STAR alignment, UMI deduplication and CDS
# read counting.
#
# Usage: riboseq/01_process_batch.sh <manifest.tsv> <raw_dir> <out_dir>
#   manifest.tsv columns (tab-separated, header): sample  barcode  fastq
#   fastq is the pooled read-1 file in <raw_dir> that contains the sample.
# Read layout: 2 nt UMI | insert (>= 20 nt) | 5 nt UMI | 5 nt sample barcode | adapter.
# Output per sample in <out_dir>:
#   <sample>.fp.fastq.gz                     extracted reads (UMI moved to read name)
#   <sample>.bowtie.fastq.gz                 reads not mapping to rRNA/tRNA
#   <sample>.Aligned.sortedByCoord.out.bam   genome alignment
#   <sample>.dedup.bam                       genome alignment, UMI-deduplicated
#   <sample>.tx.dedup.bam                    transcriptome alignment, UMI-deduplicated
#   featureCounts_cds_<sample>.dedup.tsv     CDS counts per gene
set -euo pipefail
: "${MP_REFS:?source config.env first}"

MANIFEST="$(readlink -f "$1")"; RAW="$(readlink -f "$2")"; OUT="$3"
GENOME_DIR="$MP_REFS/star_index"
RRNA_INDEX="$MP_REFS/rrna_trna/rrna_trna"
GTF="$MP_REFS/gencode.vM25.annotation.gtf"
THREADS="${MP_THREADS:-16}"
mkdir -p "$OUT"
cd "$OUT"

tail -n +2 "$MANIFEST" | while IFS=$'\t' read -r sample barcode fastq; do
    [[ -z "$sample" ]] && continue
    regex="^(?P<umi_1>.{2}).{20,}(?P<umi_2>.{5})(?P<cell_1>${barcode}){s<2}(?P<discard_1>AGATCGGAAG){s<3}(?P<discard_2>.*)"

    umi_tools extract --extract-method=regex --bc-pattern="$regex" \
        --stdin "$RAW/$fastq" --stdout "$sample.fp.fastq.gz" -L "$sample.umi.log"

    bowtie2 -p "$THREADS" --very-sensitive -U "$sample.fp.fastq.gz" -x "$RRNA_INDEX" \
        --un-gz "$sample.bowtie.fastq.gz" 2> "$sample.bowtie.log" > /dev/null

    STAR \
        --runThreadN "$THREADS" \
        --genomeDir "$GENOME_DIR" \
        --readFilesIn "$sample.bowtie.fastq.gz" \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix "$sample." \
        --outFilterMultimapNmax 20 \
        --outFilterMismatchNmax 1 \
        --quantMode TranscriptomeSAM

    samtools index "$sample.Aligned.sortedByCoord.out.bam"
    umi_tools dedup -I "$sample.Aligned.sortedByCoord.out.bam" \
        --output-stats="$sample.dedup.stats" -S "$sample.dedup.bam" -L "$sample.dedup.log"
    samtools index "$sample.dedup.bam"

    samtools sort -@ "$THREADS" -o "$sample.tx.bam" "$sample.Aligned.toTranscriptome.out.bam"
    samtools index "$sample.tx.bam"
    umi_tools dedup -I "$sample.tx.bam" \
        --output-stats="$sample.tx.dedup.stats" -S "$sample.tx.dedup.bam" -L "$sample.tx.dedup.log"

    featureCounts -t CDS -s 1 -g gene_name -a "$GTF" \
        -o "featureCounts_cds_$sample.dedup.tsv" "$sample.dedup.bam"
done
