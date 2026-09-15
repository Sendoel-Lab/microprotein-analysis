#!/usr/bin/env Rscript
# In-frame P-site percentage per pooled condition (all read lengths and 26-34 nt).
#
# Usage: Rscript riboseq/05_pooled_inframe_plot.R <pooled_dir> <pool> [<pool> ...]
# Input:  <pooled_dir>/<pool>/frame_distribution.csv
# Output: <pooled_dir>/pooled_inframe_summary.csv, <pooled_dir>/pooled_inframe_by_condition.pdf
suppressPackageStartupMessages({ library(ggplot2); library(data.table) })

args  <- commandArgs(trailingOnly = TRUE)
OUT   <- args[1]
conds <- args[-1]

rows <- rbindlist(lapply(conds, function(cc) {
  d   <- fread(file.path(OUT, cc, "frame_distribution.csv"))
  sub <- d[length >= 26 & length <= 34]
  data.table(condition    = cc,
             inframe_all  = d[is.na(length), inframe_pct],
             inframe_2634 = round(100 * sum(sub$frame0) / sum(sub$n_cds), 2),
             n_cds        = d[is.na(length), n_cds])
}))
fwrite(rows, file.path(OUT, "pooled_inframe_summary.csv"))
print(rows)

m <- melt(rows, id.vars = c("condition", "n_cds"),
          measure.vars = c("inframe_all", "inframe_2634"),
          variable.name = "set", value.name = "pct")
m[, set := factor(set, levels = c("inframe_all", "inframe_2634"),
                  labels = c("All lengths", "26-34 nt"))]
m[, condition := factor(condition, levels = conds)]

p <- ggplot(m, aes(condition, pct, fill = set)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.72) +
  geom_text(aes(label = sprintf("%.1f", pct)),
            position = position_dodge(width = 0.8), vjust = -0.4, size = 3) +
  geom_hline(yintercept = 33.3, linetype = "dashed", colour = "grey50") +
  annotate("text", x = 0.6, y = 35.5, label = "random (33%)",
           hjust = 0, size = 3, colour = "grey40") +
  scale_fill_manual(values = c("All lengths" = "#9ecae1", "26-34 nt" = "#08519c")) +
  labs(title = "In-frame P-sites by condition",
       x = NULL, y = "In-frame P-sites (%)", fill = "Read lengths") +
  ylim(0, max(m$pct) * 1.12) +
  theme_bw(base_size = 12) + theme(legend.position = "top")

ggsave(file.path(OUT, "pooled_inframe_by_condition.pdf"), p, width = 8, height = 5,
       device = pdf, useDingbats = FALSE, useKerning = FALSE)
