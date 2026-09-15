#!/usr/bin/env Rscript
# riboWaltz QC of a transcriptome-aligned ribosome-profiling BAM.
#
# Outputs in --output_dir:
#   psite_offsets.csv, length_distribution.csv, frame_distribution.csv
#   qc_plots/          length distribution, P-site region, frame, 5' end heatmap, metagene
#   periodicity_plots/ wide metagene profiles (all lengths and each length 20-34 nt)
suppressPackageStartupMessages({
  library(riboWaltz); library(data.table); library(ggplot2); library(viridis); library(optparse)
})

option_list <- list(
  make_option("--bam_dir",     type = "character", help = "Directory containing the BAM"),
  make_option("--bam_prefix",  type = "character", help = "BAM filename without .bam"),
  make_option("--sample",      type = "character", help = "Sample name used in plots"),
  make_option("--annotation",  type = "character", help = "riboWaltz annotation CSV"),
  make_option("--output_dir",  type = "character", help = "Output directory"),
  make_option("--cds_only",    type = "logical", default = FALSE, action = "store_true",
              help = "Keep only transcripts with a CDS")
)
opt <- parse_args(OptionParser(option_list = option_list))
for (arg in c("bam_dir", "bam_prefix", "sample", "annotation", "output_dir"))
  if (is.null(opt[[arg]])) stop(paste0("--", arg, " is required"))
dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)

## ---- annotation and reads ----
an <- fread(opt$annotation)
if (opt$cds_only) an <- subset(an, l_cds > 0)

sam <- c(opt$sample); names(sam) <- opt$bam_prefix
reads_list <- bamtolist(bamfolder = opt$bam_dir, annotation = an, name_samples = sam)

## ---- P-site offsets and positions ----
psite_offset <- psite(reads_list, flanking = 6, extremity = "auto")
write.csv(psite_offset, file.path(opt$output_dir, "psite_offsets.csv"), row.names = FALSE)
reads_psite_list <- psite_info(reads_list, psite_offset)
plot_name <- paste0("plot_", opt$sample)

## ---- frame distribution of CDS P-sites, per read length and overall ----
dt   <- reads_psite_list[[opt$sample]]
cdsd <- dt[psite_region == "cds" & !is.na(psite_from_start)]
cdsd[, frame := psite_from_start %% 3]
frame_tab <- cdsd[, .(n_cds = .N, frame0 = sum(frame == 0),
                      frame1 = sum(frame == 1), frame2 = sum(frame == 2)),
                  by = length][order(length)]
overall <- data.table(length = NA_integer_, n_cds = sum(frame_tab$n_cds),
                      frame0 = sum(frame_tab$frame0), frame1 = sum(frame_tab$frame1),
                      frame2 = sum(frame_tab$frame2))
frame_out <- rbind(frame_tab, overall)
frame_out[, inframe_pct := round(100 * frame0 / n_cds, 2)]
frame_out[, sample := opt$sample]
setcolorder(frame_out, c("sample", "length", "n_cds", "frame0", "frame1", "frame2", "inframe_pct"))
write.csv(frame_out, file.path(opt$output_dir, "frame_distribution.csv"), row.names = FALSE)

## ---- plot helpers (vector PDF with editable text, plus PNG) ----
save_plot <- function(dir, filename, p, w, h) {
  pdf(file.path(dir, filename), width = w, height = h, useDingbats = FALSE, useKerning = FALSE)
  print(p); dev.off()
  png(file.path(dir, sub("\\.pdf$", ".png", filename)), width = w, height = h, units = "in", res = 150)
  print(p); dev.off()
}
add_line <- function(p) p + geom_line(color = "dodgerblue4", linewidth = 1)
as_log10 <- function(mp) {
  mp[[plot_name]][["data"]][["mean_scaled_count"]] <- log10(mp[[plot_name]][["data"]][["mean_scaled_count"]])
  mp[[plot_name]][["labels"]][["y"]] <- "log10 (# P-sites)"
  mp
}

## ---- QC plots ----
qc_dir <- file.path(opt$output_dir, "qc_plots")
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
qc <- function(f, p) save_plot(qc_dir, f, p, 10, 6)

length_dist <- rlength_distr(reads_list, sample = opt$sample, colour = "dodgerblue4")
qc("01_length_dist.pdf", length_dist[[plot_name]])
qc("02_length_dist_zoom.pdf",
   rlength_distr(reads_list, sample = opt$sample, colour = "dodgerblue4", length_range = 20:60)[[plot_name]])
qc("03_psite_region.pdf",
   region_psite(reads_psite_list, an, sample = opt$sample,
                colour = c("darkorange", "darkgreen", "dodgerblue4"))[["plot"]])
qc("04_frame.pdf", frame_psite(reads_psite_list, an, sample = opt$sample, colour = "dodgerblue4")[[plot_name]])
qc("05_frame_by_length.pdf",
   frame_psite_length(reads_psite_list, an, sample = opt$sample, cl = 100, colour = "dodgerblue4")[[plot_name]])
qc("06_ends_heatmap.pdf",
   rends_heat(reads_list, an, sample = opt$sample, cl = 100, utr5l = 25, cdsl = 50, utr3l = 25,
              colour = viridis(20))[[plot_name]])
mp <- metaprofile_psite(reads_psite_list, an, sample = opt$sample,
                        utr5l = 20, cdsl = 50, utr3l = 20, colour = "dodgerblue4")
qc("07_metagene_scaled.pdf", add_line(mp[[plot_name]]))
mp_uns <- metaprofile_psite(reads_psite_list, an, sample = opt$sample,
                            utr5l = 20, cdsl = 50, utr3l = 20, colour = "dodgerblue4", scale_factors = "none")
qc("08_metagene_raw.pdf", add_line(mp_uns[[plot_name]]))
qc("09_metagene_log10.pdf", add_line(as_log10(mp_uns)[[plot_name]]))

## ---- periodicity plots (50 nt UTR / 200 nt CDS windows) ----
per_dir <- file.path(opt$output_dir, "periodicity_plots")
dir.create(per_dir, recursive = TRUE, showWarnings = FALSE)
per <- function(f, p) save_plot(per_dir, f, p, 25, 6)

mp_all <- metaprofile_psite(reads_psite_list, an, sample = opt$sample,
                            utr5l = 50, cdsl = 200, utr3l = 50, colour = "dodgerblue4")
per("metagene_all_scaled.pdf", add_line(mp_all[[plot_name]]))
mp_all_uns <- metaprofile_psite(reads_psite_list, an, sample = opt$sample,
                                utr5l = 50, cdsl = 200, utr3l = 50, colour = "dodgerblue4", scale_factors = "none")
per("metagene_all_raw.pdf", add_line(mp_all_uns[[plot_name]]))
per("metagene_all_log10.pdf", add_line(as_log10(mp_all_uns)[[plot_name]]))

for (i in 20:34) {
  f <- sprintf("metagene_%dnt.pdf", i)
  tryCatch({
    mp_i <- metaprofile_psite(reads_psite_list, an, sample = opt$sample,
                              utr5l = 50, cdsl = 200, utr3l = 50, colour = "dodgerblue4",
                              scale_factors = "none", length_range = i)
    mp_i[[plot_name]][["labels"]][["title"]] <- i
    per(f, add_line(mp_i[[plot_name]]))
  }, error = function(e) {
    pdf(file.path(per_dir, f), width = 25, height = 6, useDingbats = FALSE, useKerning = FALSE)
    plot.new(); text(0.5, 0.5, sprintf("No reads at %d nt", i), cex = 3); dev.off()
  })
}

write.csv(length_dist[["count_dt"]], file.path(opt$output_dir, "length_distribution.csv"), row.names = FALSE)
