#!/usr/bin/env Rscript
# Overlap of differentially expressed gene sets between two contrasts, per
# category (homodirectional / transcriptional / translational) and direction,
# with a one-sided Fisher's exact test on the shared testable universe.
#
# Usage:
#   Rscript de_te/venn_overlap.R --a <merged_a.tsv> --b <merged_b.tsv> \
#       --label_a <name> --label_b <name> --out <dir>
# Output: <dir>/venn_<category>_<direction>.{pdf,png}, venn_overview.{pdf,png}, venn_stats.tsv
suppressPackageStartupMessages({ library(ggvenn); library(ggplot2); library(cowplot); library(optparse) })

opt <- parse_args(OptionParser(option_list = list(
  make_option("--a",       type = "character"),
  make_option("--b",       type = "character"),
  make_option("--label_a", type = "character"),
  make_option("--label_b", type = "character"),
  make_option("--out",     type = "character")
)))
for (x in c("a", "b", "label_a", "label_b", "out")) if (is.null(opt[[x]])) stop(paste0("--", x, " is required"))
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)

CATEGORIES <- c(Homodirectional = "HD", Transcriptional = "Tx", Translational = "Tl")

load_contrast <- function(path) {
  m <- read.table(path, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
  list(table = m, universe = m$gene_name[!is.na(m$padj_RNA) | !is.na(m$padj_Ribo)])
}
get_set <- function(m, cat_name, direction) {
  sub <- m[m$category == cat_name, ]
  lfc <- if (cat_name == "Translational") sub$log2FoldChange_Ribo else sub$log2FoldChange_RNA
  if (direction == "up") sub$gene_name[lfc > 0] else sub$gene_name[lfc < 0]
}
fisher_overlap <- function(set_a, set_b, universe) {
  a <- intersect(set_a, universe); b <- intersect(set_b, universe)
  n11 <- length(intersect(a, b)); n10 <- length(a) - n11; n01 <- length(b) - n11
  n00 <- length(universe) - n11 - n10 - n01
  ft <- fisher.test(matrix(c(n11, n10, n01, n00), nrow = 2), alternative = "greater")
  list(n_a = length(a), n_b = length(b), overlap = n11, universe = length(universe),
       odds = unname(ft$estimate), p = ft$p.value)
}

A <- load_contrast(opt$a); B <- load_contrast(opt$b)
universe <- intersect(A$universe, B$universe)
cat(sprintf("Shared universe: %d genes\n", length(universe)))

buckets <- list()
for (cat_name in names(CATEGORIES)) for (direction in c("up", "down")) {
  tag <- sprintf("%s_%s", tolower(cat_name), direction)
  set_a <- get_set(A$table, cat_name, direction); set_b <- get_set(B$table, cat_name, direction)
  buckets[[tag]] <- list(cat_name = cat_name, direction = direction, set_a = set_a, set_b = set_b,
                         fish = fisher_overlap(set_a, set_b, universe))
}
adj <- p.adjust(sapply(buckets, function(b) b$fish$p), method = "BH")
for (tag in names(buckets)) buckets[[tag]]$fish$padj <- adj[[tag]]

stats <- data.frame(); plots <- list()
for (tag in names(buckets)) {
  b <- buckets[[tag]]; f <- b$fish
  d <- setNames(list(intersect(b$set_a, universe), intersect(b$set_b, universe)), c(opt$label_a, opt$label_b))
  p <- ggvenn(d, show_percentage = FALSE, fill_color = c("#3F6FB0", "#9B1B1B"), fill_alpha = 0.45,
              stroke_color = "black", stroke_size = 0.5, set_name_size = 5, text_size = 5) +
    labs(title = sprintf("%s %s", CATEGORIES[[b$cat_name]], b$direction),
         subtitle = sprintf("overlap = %d  |  odds = %.2f  |  universe = %d\nFisher p = %.2e  |  BH p = %.2e",
                            f$overlap, f$odds, f$universe, f$p, f$padj)) +
    theme(plot.title = element_text(size = 13, hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(size = 10, hjust = 0.5),
          plot.background = element_rect(fill = "white", color = NA))
  ggsave(file.path(opt$out, sprintf("venn_%s.png", tag)), p, width = 6, height = 6, dpi = 150, bg = "white")
  ggsave(file.path(opt$out, sprintf("venn_%s.pdf", tag)), p, width = 6, height = 6, bg = "white",
         device = pdf, useDingbats = FALSE, useKerning = FALSE)
  plots[[tag]] <- p
  stats <- rbind(stats, data.frame(category = b$cat_name, direction = b$direction,
                                   n_a = f$n_a, n_b = f$n_b, overlap = f$overlap, universe = f$universe,
                                   odds = round(f$odds, 3), pvalue = signif(f$p, 4), padj_BH = signif(f$padj, 4)))
}
write.table(stats, file.path(opt$out, "venn_stats.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(stats)

overview <- plot_grid(plotlist = plots[c("homodirectional_up", "homodirectional_down",
                                         "transcriptional_up", "transcriptional_down",
                                         "translational_up", "translational_down")],
                      ncol = 2, nrow = 3) +
  theme(plot.background = element_rect(fill = "white", color = NA))
ggsave(file.path(opt$out, "venn_overview.png"), overview, width = 12, height = 18, dpi = 150, bg = "white")
ggsave(file.path(opt$out, "venn_overview.pdf"), overview, width = 12, height = 18, bg = "white",
       device = pdf, useDingbats = FALSE, useKerning = FALSE)
