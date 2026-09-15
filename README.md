# microprotein-analysis

Analysis code for the sequencing experiments of the Gm10076 / Rpl41 microprotein study:
RNA-seq and ribosome-profiling processing, differential expression and translational
efficiency, gene-set enrichment, and single-cell CRISPR screen (perturb-seq) summaries.

## Layout

| Directory | Content |
|---|---|
| `rnaseq/` | fastp trimming, STAR alignment, featureCounts gene counts |
| `riboseq/` | UMI extraction and demultiplexing, rRNA depletion, STAR, UMI deduplication, CDS counts, riboWaltz QC, pooled in-frame summary |
| `de_te/` | DESeq2 on matched RNA/Ribo libraries (`~ batch + genotype`), RNA-vs-Ribo scatter and volcano plots, Enrichr GO enrichment, overlap between contrasts |
| `perturbseq/` | volcano plots and DEG burden per ORF class from glmGamPoi differential expression |
| `envs/` | conda environments (`seq.yml` for `rnaseq/`, `riboseq/`, `de_te/`; `perturbseq.yml` for `perturbseq/`) |

Sample sheets and manifests are in `riboseq/manifests/` and `de_te/samples/`.

## Setup

```bash
cp config.env.example config.env   # edit paths
source config.env
mamba env create -f envs/seq.yml
mamba env create -f envs/perturbseq.yml
```

Expected under `MP_REFS`: `star_index/` (GENCODE vM25, STAR 2.7.7a), `gencode.vM25.annotation.gtf`,
`rrna_trna/rrna_trna` (bowtie2 index), `ribowaltz_annotation.csv`.
Expected under `MP_EXTERNAL`: `riboseq_initial/` (transcriptome BAMs and CDS count tables of the
initial ribosome-profiling cohort), `perturbseq/glm/` (glmGamPoi tables per cluster),
`perturbseq/sgRNAs_ORFtype.xlsx`.

Ribosome-profiling libraries were processed with umi_tools 1.1.1 and featureCounts 2.0.1;
RNA-seq counting used featureCounts 2.1.1, which `envs/seq.yml` installs.

## Workflow

```bash
# RNA-seq
rnaseq/01_trim.sh samples.tsv
rnaseq/02_align.sh <sample> ...
rnaseq/03_count.sh $MP_DATA/rnaseq/counts/counts.tsv <sample> ...

# Ribo-seq
riboseq/01_process_batch.sh riboseq/manifests/run1.tsv $MP_DATA/riboseq/raw/run1 $MP_DATA/riboseq/run1
riboseq/01_process_batch.sh riboseq/manifests/run2.tsv $MP_DATA/riboseq/raw/run2 $MP_DATA/riboseq/run2
riboseq/02_build_count_matrix.py --batch run1 --indir $MP_DATA/riboseq/run1 --outdir $MP_DATA/riboseq/counts
riboseq/02_build_count_matrix.py --batch run2 --indir $MP_DATA/riboseq/run2 --outdir $MP_DATA/riboseq/counts
riboseq/04_pool_conditions.sh riboseq/manifests/pools.tsv $MP_DATA/riboseq/ribowaltz_pooled

# Differential expression / translational efficiency
Rscript de_te/deseq2_te.R --samples de_te/samples/gm10076_6v6.tsv --rna_counts $MP_DATA/rnaseq/counts/counts.tsv \
    --reference Ctrl --exclude Rpl41 --title "Gm10076 KO" --out $MP_RESULTS/te_gm10076_6v6
Rscript de_te/deseq2_te.R --samples de_te/samples/rpl41_2v2.tsv --rna_counts $MP_DATA/rnaseq/counts/counts.tsv \
    --reference mCherry --exclude Rpl41 --title "Rpl41 KO" --out $MP_RESULTS/te_rpl41_2v2
Rscript de_te/enrichr_go.R --merged $MP_RESULTS/te_gm10076_6v6/merged_RNA_vs_Ribo_C7KO_vs_Ctrl.tsv \
    --label gm10076 --out $MP_RESULTS/enrichr
Rscript de_te/enrichr_go.R --merged $MP_RESULTS/te_rpl41_2v2/merged_RNA_vs_Ribo_Rpl41F10_vs_mCherry.tsv \
    --label rpl41 --out $MP_RESULTS/enrichr
Rscript de_te/venn_overlap.R --a $MP_RESULTS/te_gm10076_6v6/merged_RNA_vs_Ribo_C7KO_vs_Ctrl.tsv \
    --b $MP_RESULTS/te_rpl41_2v2/merged_RNA_vs_Ribo_Rpl41F10_vs_mCherry.tsv \
    --label_a "Gm10076 KO" --label_b "Rpl41 KO" --out $MP_RESULTS/venn

# Perturb-seq
perturbseq/volcano_glm.py --cluster 0 --perturbations Emg1 Gm13935 --lfc 0.5
perturbseq/deg_burden_by_orf_class.py
```

## Conventions

- DESeq2: Wald test, `independentFiltering = FALSE`, no fold-change shrinkage; genes with
  RPKM >= 2 in all RNA libraries of a contrast and present in the Ribo count matrix.
- Significance: `padj < 0.05` and `|log2FC| > 0.25` per assay. Homodirectional = significant in
  both assays with the same sign; transcriptional = RNA only; translational = Ribo only.
- Enrichr: one-sided Fisher's exact test against the library background, Benjamini-Hochberg
  adjusted, top 5 terms by raw p-value.
- Vector figures are written as PDF with editable text.
