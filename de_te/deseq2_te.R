#!/usr/bin/env Rscript
# DESeq2 differential expression of matched RNA-seq and ribosome-profiling
# libraries (design ~ batch + genotype), with an RNA-vs-Ribo log2FC scatter
# and per-assay volcano plots for every genotype vs the reference.
#
# Usage:
#   Rscript de_te/deseq2_te.R --samples <sheet.tsv> --rna_counts <counts.tsv> \
#       --reference <genotype> --out <dir> [--exclude Rpl41] [--rpkm 2] [--title "Gm10076 KO"]
# Sample sheet columns: rna_sample  ribo_sample  ribo_counts  batch  genotype
#   ribo_counts: a count matrix with a <ribo_sample> column, or a single-sample
#   featureCounts table (last column used); ${VAR} is expanded from the environment.
# Outputs per contrast <c> = <genotype>_vs_<reference>:
#   deseq2_RNA_<c>.tsv, deseq2_Ribo_<c>.tsv, merged_RNA_vs_Ribo_<c>.tsv,
#   scatter_RNA_vs_Ribo_<c>.{pdf,png}, volcano_{RNA,Ribo}_<c>.{pdf,png}
suppressPackageStartupMessages({
  library(DESeq2); library(ggplot2); library(ggrepel); library(cowplot); library(optparse)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--samples",    type = "character"),
  make_option("--rna_counts", type = "character"),
  make_option("--reference",  type = "character"),
  make_option("--out",        type = "character"),
  make_option("--exclude",    type = "character", default = "", help = "Comma-separated genes dropped from scatter"),
  make_option("--rpkm",       type = "double",    default = 2),
  make_option("--title",      type = "character", default = ""),
  make_option("--padj",       type = "double",    default = 0.05),
  make_option("--lfc",        type = "double",    default = 0.25)
)))
for (a in c("samples", "rna_counts", "reference", "out"))
  if (is.null(opt[[a]])) stop(paste0("--", a, " is required"))
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)

expand_env <- function(x) {
  m <- gregexpr("\\$\\{[A-Za-z_][A-Za-z0-9_]*\\}", x)
  regmatches(x, m) <- lapply(regmatches(x, m), function(v) Sys.getenv(gsub("[${}]", "", v)))
  x
}

## ---- sample sheet ----
S <- read.delim(opt$samples, stringsAsFactors = FALSE)
S$ribo_counts <- expand_env(S$ribo_counts)
genotypes <- c(opt$reference, setdiff(unique(S$genotype), opt$reference))
batches   <- unique(S$batch)

## ---- RNA counts, RPKM filter ----
rna_raw  <- read.table(opt$rna_counts, header = TRUE, sep = "\t", comment.char = "#",
                       check.names = FALSE, row.names = 1)
stopifnot(all(S$rna_sample %in% colnames(rna_raw)))
rna  <- as.matrix(rna_raw[, S$rna_sample])
rpkm <- sweep(sweep(rna, 1, rna_raw$Length, "/") * 1e3, 2, colSums(rna), "/") * 1e6
keep <- rowSums(rpkm >= opt$rpkm) == ncol(rpkm)
cat(sprintf("RNA: %d genes, %d with RPKM >= %g in all %d samples\n",
            nrow(rna), sum(keep), opt$rpkm, ncol(rna)))

## ---- Ribo counts ----
read_ribo <- function(path, sample) {
  d <- read.table(path, header = TRUE, sep = "\t", comment.char = "#",
                  check.names = FALSE, row.names = 1)
  if (sample %in% colnames(d)) d[, sample, drop = FALSE] else d[, ncol(d), drop = FALSE]
}
ribo_cols <- Map(read_ribo, S$ribo_counts, S$ribo_sample)
ribo_genes <- Reduce(intersect, lapply(ribo_cols, rownames))
ribo <- do.call(cbind, lapply(ribo_cols, function(d) d[ribo_genes, , drop = FALSE]))
colnames(ribo) <- S$ribo_sample
ribo <- as.matrix(ribo)

genes <- intersect(rownames(rna)[keep], rownames(ribo))
cat(sprintf("Genes in both assays after filtering: %d\n", length(genes)))
rna  <- rna[genes, ]
ribo <- ribo[genes, ]

## ---- DESeq2 per assay ----
run_assay <- function(cts, samples, label) {
  cd <- data.frame(batch    = factor(S$batch, levels = batches),
                   genotype = factor(S$genotype, levels = genotypes),
                   row.names = samples)
  dds <- DESeq(DESeqDataSetFromMatrix(countData = cts, colData = cd, design = ~ batch + genotype))
  res <- list()
  for (g in genotypes[-1]) {
    contrast <- sprintf("%s_vs_%s", g, opt$reference)
    df <- as.data.frame(results(dds, name = paste0("genotype_", contrast),
                                alpha = opt$padj, independentFiltering = FALSE))
    df$gene_name <- rownames(df)
    df <- df[order(df$padj, df$pvalue), c("gene_name", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj")]
    cat(sprintf("[%s %s] %d genes, %d with padj < %g\n", label, contrast, nrow(df),
                sum(!is.na(df$padj) & df$padj < opt$padj), opt$padj))
    write.table(df, file.path(opt$out, sprintf("deseq2_%s_%s.tsv", label, contrast)),
                sep = "\t", quote = FALSE, row.names = FALSE)
    res[[contrast]] <- df
  }
  res
}
rna_res  <- run_assay(rna,  S$rna_sample,  "RNA")
ribo_res <- run_assay(ribo, S$ribo_sample, "Ribo")

## ---- RNA vs Ribo scatter with marginal histograms ----
PAL <- c("Not significant" = "grey75", "Transcriptional" = "#3F6FB0",
         "Homodirectional" = "#E68A1E", "Translational" = "#9B1B1B", "Opposite" = "#5A3D8C")
AX  <- c(-6.5, 6.5)
exclude <- strsplit(opt$exclude, ",")[[1]]

make_scatter <- function(rna_df, ribo_df, contrast, n_case, n_ref) {
  m <- merge(rna_df [, c("gene_name", "log2FoldChange", "pvalue", "padj")],
             ribo_df[, c("gene_name", "log2FoldChange", "pvalue", "padj")],
             by = "gene_name", suffixes = c("_RNA", "_Ribo"))
  m <- m[is.finite(m$log2FoldChange_RNA) & is.finite(m$log2FoldChange_Ribo) & !(m$gene_name %in% exclude), ]
  pr <- cor(m$log2FoldChange_RNA, m$log2FoldChange_Ribo, method = "pearson")
  sp <- cor(m$log2FoldChange_RNA, m$log2FoldChange_Ribo, method = "spearman")

  sig_rna  <- !is.na(m$padj_RNA)  & m$padj_RNA  < opt$padj & abs(m$log2FoldChange_RNA)  > opt$lfc
  sig_ribo <- !is.na(m$padj_Ribo) & m$padj_Ribo < opt$padj & abs(m$log2FoldChange_Ribo) > opt$lfc
  same_dir <- sign(m$log2FoldChange_RNA) == sign(m$log2FoldChange_Ribo)
  m$category <- factor(
    ifelse(sig_rna & sig_ribo &  same_dir, "Homodirectional",
    ifelse(sig_rna & sig_ribo & !same_dir, "Opposite",
    ifelse(sig_rna,                        "Transcriptional",
    ifelse(sig_ribo,                       "Translational", "Not significant")))),
    levels = names(PAL))

  m$rank_score <- pmin(replace(m$padj_RNA, is.na(m$padj_RNA), 1),
                       replace(m$padj_Ribo, is.na(m$padj_Ribo), 1))
  m$dist2 <- m$log2FoldChange_RNA^2 + m$log2FoldChange_Ribo^2
  top_by <- function(cat, n) { d <- m[m$category == cat, ]; head(d[order(d$rank_score), ], n) }
  sig <- m[m$category != "Not significant", ]
  ns  <- m[m$category == "Not significant", ]
  d_lab <- unique(rbind(top_by("Homodirectional", 25), top_by("Transcriptional", 15),
                        top_by("Translational", 15),
                        head(sig[order(-sig$dist2), ], 30), head(ns[order(-ns$dist2), ], 5)))
  m <- m[order(m$category), ]

  scatter <- ggplot(m, aes(x = log2FoldChange_RNA, y = log2FoldChange_Ribo, color = category)) +
    geom_vline(xintercept = c(-opt$lfc, opt$lfc), linetype = "dotted", color = "grey50") +
    geom_hline(yintercept = c(-opt$lfc, opt$lfc), linetype = "dotted", color = "grey50") +
    geom_point(size = 1.0, alpha = 0.75) +
    geom_text_repel(data = d_lab, aes(label = gene_name), color = "black", size = 3.2,
                    max.overlaps = Inf, min.segment.length = 0, box.padding = 0.3,
                    fontface = "italic", show.legend = FALSE) +
    scale_color_manual(values = PAL, name = NULL,
                       breaks = c("Transcriptional", "Homodirectional", "Translational", "Opposite", "Not significant")) +
    scale_x_continuous(limits = AX, expand = c(0, 0)) +
    scale_y_continuous(limits = AX, expand = c(0, 0)) +
    labs(x = sprintf("Transcriptome log2FC (%s)", contrast),
         y = sprintf("Translatome log2FC (%s)", contrast),
         subtitle = sprintf("%d genes (%d vs %d)  |  Pearson r = %.3f  |  Spearman rho = %.3f",
                            nrow(m), n_case, n_ref, pr, sp)) +
    guides(color = guide_legend(override.aes = list(size = 4, alpha = 1), ncol = 1)) +
    theme_bw(base_size = 12) +
    theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1),
          legend.background = element_rect(fill = alpha("white", 0.9), color = NA),
          legend.key = element_blank(), legend.text = element_text(size = 12),
          legend.spacing.y = unit(2, "pt"),
          plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA))

  hist_theme <- theme_bw(base_size = 9) +
    theme(plot.margin = margin(2, 2, 2, 2), legend.position = "none",
          panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA),
          panel.grid.major = element_line(color = "grey88", linewidth = 0.3),
          panel.grid.minor = element_line(color = "grey94", linewidth = 0.2))
  top_hist <- ggplot(m, aes(x = log2FoldChange_RNA, fill = category)) +
    geom_histogram(bins = 80, position = "stack", color = "black", linewidth = 0.15) +
    scale_fill_manual(values = PAL, guide = "none") +
    scale_x_continuous(limits = AX, expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0), n.breaks = 3) + labs(x = NULL, y = NULL) + hist_theme
  right_hist <- ggplot(m, aes(y = log2FoldChange_Ribo, fill = category)) +
    geom_histogram(bins = 80, position = "stack", color = "black", linewidth = 0.15, orientation = "y") +
    scale_fill_manual(values = PAL, guide = "none") +
    scale_y_continuous(limits = AX, expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0), n.breaks = 3) + labs(x = NULL, y = NULL) + hist_theme
  corner <- ggplot() + theme_void() + annotate("text", x = 0.5, y = 0.5, label = "Genes", size = 4) +
    theme(plot.background = element_rect(fill = "white", color = NA))
  title <- ggdraw() + draw_label(
    sprintf("Transcriptional vs translational changes (%s)", if (nzchar(opt$title)) opt$title else contrast),
    size = 14, fontface = "bold")
  body <- plot_grid(top_hist, corner, scatter, right_hist, ncol = 2, nrow = 2,
                    rel_widths = c(10, 2), rel_heights = c(2.15, 10), align = "hv", axis = "tblr")
  fig <- plot_grid(title, body, ncol = 1, rel_heights = c(0.06, 1)) +
    theme(plot.background = element_rect(fill = "white", color = NA))

  stem <- file.path(opt$out, sprintf("scatter_RNA_vs_Ribo_%s", contrast))
  ggsave(paste0(stem, ".png"), fig, width = 11, height = 11, dpi = 150, bg = "white")
  ggsave(paste0(stem, ".pdf"), fig, width = 11, height = 11, bg = "white",
         device = pdf, useDingbats = FALSE, useKerning = FALSE)
  write.table(m[, setdiff(colnames(m), "dist2")], file.path(opt$out, sprintf("merged_RNA_vs_Ribo_%s.tsv", contrast)),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("[scatter %s] r = %.3f, rho = %.3f\n", contrast, pr, sp))
  print(table(m$category))
}

## ---- volcano per assay ----
volcano <- function(df, assay, contrast, top_n = 25) {
  d <- df[!is.na(df$padj), ]
  d$sig <- ifelse(d$padj < opt$padj & abs(d$log2FoldChange) > opt$lfc,
                  ifelse(d$log2FoldChange > 0, "up", "down"), "ns")
  d_lab <- head(d[d$sig != "ns", ][order(d[d$sig != "ns", ]$padj), ], top_n)
  p <- ggplot(d, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = sig), size = 1.1, alpha = 0.7) +
    geom_vline(xintercept = c(-opt$lfc, opt$lfc), linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = -log10(opt$padj), linetype = "dashed", color = "grey50") +
    scale_color_manual(values = c(ns = "grey70", up = "firebrick", down = "steelblue3")) +
    geom_text_repel(data = d_lab, aes(label = gene_name), size = 3, max.overlaps = Inf, min.segment.length = 0) +
    labs(title = sprintf("%s: %s", assay, contrast),
         subtitle = sprintf("dashed = |LFC| > %g & padj < %g", opt$lfc, opt$padj),
         x = "log2 fold change", y = "-log10(padj)") +
    theme_bw(base_size = 12)
  stem <- file.path(opt$out, sprintf("volcano_%s_%s", assay, contrast))
  ggsave(paste0(stem, ".png"), p, width = 8, height = 6, dpi = 150, bg = "white")
  ggsave(paste0(stem, ".pdf"), p, width = 8, height = 6, bg = "white",
         device = pdf, useDingbats = FALSE, useKerning = FALSE)
}

for (contrast in names(rna_res)) {
  g <- sub(sprintf("_vs_%s$", opt$reference), "", contrast)
  make_scatter(rna_res[[contrast]], ribo_res[[contrast]], contrast,
               sum(S$genotype == g), sum(S$genotype == opt$reference))
  volcano(rna_res[[contrast]],  "RNA",  contrast)
  volcano(ribo_res[[contrast]], "Ribo", contrast)
}
