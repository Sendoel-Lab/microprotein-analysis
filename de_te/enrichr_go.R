#!/usr/bin/env Rscript
# Enrichr over-representation analysis (default Enrichr background) of the
# transcriptional, homodirectional and translational gene sets of one contrast,
# separately for up-, down-regulated and combined genes; top terms are plotted.
#
# Usage:
#   Rscript de_te/enrichr_go.R --merged <merged_RNA_vs_Ribo_<c>.tsv> --label <name> --out <dir> \
#       [--library GO_Biological_Process_2026] [--top 5]
# Output per set: <dir>/<label>_<category>_<direction>_<library>.{tsv,pdf,png}
suppressPackageStartupMessages({ library(enrichR); library(ggplot2); library(optparse) })

opt <- parse_args(OptionParser(option_list = list(
  make_option("--merged",  type = "character"),
  make_option("--label",   type = "character"),
  make_option("--out",     type = "character"),
  make_option("--library", type = "character", default = "GO_Biological_Process_2026"),
  make_option("--top",     type = "integer",   default = 5)
)))
for (a in c("merged", "label", "out")) if (is.null(opt[[a]])) stop(paste0("--", a, " is required"))
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
setEnrichrSite("Enrichr")

CATEGORIES <- c("Homodirectional", "Transcriptional", "Translational")
lib_tag <- tolower(gsub("[^A-Za-z0-9]+", "_", opt$library))

run_set <- function(cat_name, direction, genes) {
  if (length(genes) < 5) return(invisible(NULL))
  df <- enrichr(unique(genes), databases = opt$library)[[opt$library]]
  if (is.null(df) || nrow(df) == 0) return(invisible(NULL))
  df <- df[order(df$P.value), ]
  df$Rank <- seq_len(nrow(df))
  n_fdr <- sum(df$Adjusted.P.value < 0.05, na.rm = TRUE)
  cat(sprintf("[%s %s] n = %d genes, %d terms, %d with FDR < 0.05\n",
              cat_name, direction, length(genes), nrow(df), n_fdr))

  base <- sprintf("%s_%s_%s_%s", opt$label, tolower(cat_name), tolower(direction), lib_tag)
  write.table(df, file.path(opt$out, paste0(base, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)

  p_df <- head(df, opt$top)
  p_df$Count      <- as.integer(sapply(strsplit(p_df$Overlap, "/"), `[`, 1))
  p_df$nlog10p    <- -log10(p_df$P.value)
  p_df$nlog10padj <- -log10(p_df$Adjusted.P.value)
  p_df$fdr        <- ifelse(p_df$Adjusted.P.value < 0.05, "FDR<0.05", "FDR>=0.05")
  p_df$Term       <- factor(p_df$Term, levels = rev(p_df$Term))
  p <- ggplot(p_df, aes(x = nlog10p, y = Term, size = Count, color = nlog10padj, shape = fdr)) +
    geom_point() +
    scale_color_gradient(low = "steelblue3", high = "firebrick", name = "-log10(adj.p)") +
    scale_size_continuous(range = c(2, 8), name = "Count") +
    scale_shape_manual(values = c("FDR<0.05" = 16, "FDR>=0.05" = 1), name = NULL) +
    labs(title = sprintf("%s: %s %s genes (n = %d)", opt$label, cat_name, direction, length(genes)),
         subtitle = sprintf("%s, ranked by raw p-value, %d/%d with FDR < 0.05", opt$library, n_fdr, nrow(df)),
         x = "-log10(raw p-value)", y = NULL) +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(size = 11, hjust = 0.5),
          plot.subtitle = element_text(size = 9, hjust = 0.5), plot.title.position = "plot")
  ggsave(file.path(opt$out, paste0(base, ".png")), p, width = 11, height = 7, dpi = 150)
  ggsave(file.path(opt$out, paste0(base, ".pdf")), p, width = 11, height = 7, bg = "white",
         device = pdf, useDingbats = FALSE, useKerning = FALSE)
}

m <- read.table(opt$merged, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
for (cat_name in CATEGORIES) {
  sub <- m[m$category == cat_name, ]
  lfc <- if (cat_name == "Translational") sub$log2FoldChange_Ribo else sub$log2FoldChange_RNA
  run_set(cat_name, "Up",       sub$gene_name[lfc > 0])
  run_set(cat_name, "Down",     sub$gene_name[lfc < 0])
  run_set(cat_name, "Combined", sub$gene_name)
}
