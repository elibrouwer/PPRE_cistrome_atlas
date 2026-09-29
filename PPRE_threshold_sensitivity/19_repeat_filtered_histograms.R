# Repeat the fig7/fig8/fig9 hits-per-gene histograms after dropping every
# FIMO-derived deduplicated interval (p<=1e-4, combined view -- same intervals
# used throughout 05/13/18) that overlaps an annotated UCSC RepeatMasker
# element (any class: SINE/Alu, LINE, LTR, DNA transposon, simple repeat,
# satellite, low-complexity). Motivation: nuclear-receptor half-sites are
# documented to be spread genome-wide by retrotransposon (esp. Alu) expansion
# for related motifs (DR2/RAR, HNF4a), so this checks whether the multi-PPRE
# signal seen in fig7-fig9 survives once repeat-embedded hits are removed.
#
# rmsk_hg38.txt.gz is the full UCSC hg38 RepeatMasker table
# (https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz),
# downloaded once into data/ alongside this pipeline's other cached resources.

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(ggplot2)
  library(scales)
})
source("./PPRE_threshold_sensitivity/00_config.R")

marker_genes <- c("CPT1A", "CPT1B", "HADHA", "HADHB", "ACADVL")
std_chroms <- paste0("chr", c(1:22, "X", "Y", "M"))

make_marked_histogram <- function(d, value_col, binwidth, x_lab, y_lab = "number of genes", title) {
  markers <- d[gene %in% marker_genes, c("gene", value_col), with = FALSE]
  setnames(markers, value_col, "value")
  missing_markers <- setdiff(marker_genes, markers$gene)
  if (length(missing_markers) > 0)
    warning("Marker gene(s) not found: ", paste(missing_markers, collapse = ", "))
  setorder(markers, value)

  max_count <- max(hist(d[[value_col]], breaks = seq(min(d[[value_col]]), max(d[[value_col]]) + binwidth, by = binwidth), plot = FALSE)$counts)
  markers[, y_label := max_count * seq(1.22, 0.9, length.out = .N)[frank(value, ties.method = "first")]]
  label_fmt <- if (is.integer(d[[value_col]])) "%s (%d)" else "%s (%.2f)"

  ggplot(d, aes(x = .data[[value_col]])) +
    geom_histogram(binwidth = binwidth, fill = "#3B7FB6", color = "white", linewidth = 0.2) +
    geom_vline(data = markers, aes(xintercept = value), color = "#D6604D", linetype = "dashed", linewidth = 0.5) +
    geom_point(data = markers, aes(x = value, y = y_label), color = "#B2182B", size = 1.2) +
    geom_text(data = markers, aes(x = value, y = y_label, label = sprintf(label_fmt, gene, value)),
              color = "#B2182B", size = 3, hjust = -0.08, vjust = 0.35) +
    scale_x_continuous(breaks = pretty_breaks(10)) +
    scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.22))) +
    coord_cartesian(clip = "off") +
    labs(x = x_lab, y = y_lab, title = title) +
    theme_minimal(base_size = 11)
}

## ---- load / cache RepeatMasker as a single merged mask (class info not needed for presence/absence) ----
rmsk_gz  <- file.path(out_data, "rmsk_hg38.txt.gz")
rmsk_rds <- file.path(out_data, "rmsk_hg38_merged_gr.rds")
if (file.exists(rmsk_rds)) {
  repeat_gr <- readRDS(rmsk_rds)
} else {
  stopifnot(file.exists(rmsk_gz))
  rmsk <- fread(rmsk_gz, header = FALSE, select = c(6, 7, 8),
                col.names = c("chrom", "start0", "end"))
  rmsk <- rmsk[chrom %in% std_chroms]
  repeat_gr <- reduce(GRanges(rmsk$chrom, IRanges(rmsk$start0 + 1L, rmsk$end)))
  saveRDS(repeat_gr, rmsk_rds)
}
cat(sprintf("RepeatMasker mask: %s merged intervals, %s bp (%.1f%% of the genome)\n",
            format(length(repeat_gr), big.mark = ","), format(sum(width(repeat_gr)), big.mark = ","),
            100 * sum(as.numeric(width(repeat_gr))) / 3.05e9))

## ---- deduplicated FIMO intervals (combined view, p<=1e-4 baseline), repeats removed ----
interval_list <- readRDS(file.path(out_data, "deduplicated_intervals_by_threshold.rds"))
tag <- paste0("p", format(baseline_threshold, scientific = TRUE))
itv <- interval_list[[tag]]
gr_itv <- GRanges(itv$chrom, IRanges(itv$start, itv$end))

in_repeat <- overlapsAny(gr_itv, repeat_gr, ignore.strand = TRUE)
cat(sprintf("Deduplicated genome-wide intervals: %s total, %s (%.1f%%) overlap a RepeatMasker element -- dropped\n",
            format(length(gr_itv), big.mark = ","), format(sum(in_repeat), big.mark = ","), 100 * mean(in_repeat)))
gr_itv_nr <- gr_itv[!in_repeat]

## =========================== fig10: raw hits per gene, repeat-filtered (fig7 equivalent) ===========================
read_promoter_bed <- function(bed) {
  df <- fread(bed, header = FALSE, col.names = c("chrom", "start0", "end", "gene", "score", "strand"))
  GRanges(df$chrom, IRanges(df$start0 + 1L, df$end), gene = df$gene)
}
promo_var <- read_promoter_bed(promoter_windows_nonzero_bed)   # variable-length, union-across-isoforms windows (same as 02/05)

count_hits <- function(gr_itv_set, promo_gr) {
  ov <- findOverlaps(gr_itv_set, promo_gr, ignore.strand = TRUE)
  dt <- data.table(gene = promo_gr$gene[subjectHits(ov)], interval_idx = queryHits(ov))[
    , .(n_hits = uniqueN(interval_idx)), by = gene]
  all_genes <- data.table(gene = unique(promo_gr$gene))
  dt <- merge(all_genes, dt, by = "gene", all.x = TRUE)
  dt[is.na(n_hits), n_hits := 0L]
  dt
}

d10 <- count_hits(gr_itv_nr, promo_var)
cat(sprintf("fig10: repeat-filtered raw hits per gene (variable window): n=%s, min=%d, median=%.0f, mean=%.1f, max=%d\n",
            format(nrow(d10), big.mark = ","), min(d10$n_hits), median(d10$n_hits), mean(d10$n_hits), max(d10$n_hits)))

p10 <- make_marked_histogram(d10, "n_hits", binwidth = 5,
                              x_lab = "repeat-filtered FIMO motif hits per gene (promoter-overlapping intervals, 30 motifs combined)",
                              title = "Distribution of FIMO motif hits per gene, repeat elements removed")
ggsave(file.path(out_figures, "fig10_hits_per_gene_repeat_filtered.png"), p10, width = 8, height = 5.5, dpi = 150, bg = "white")
fwrite(d10[order(-n_hits)], file.path(out_results, "hits_per_gene_repeat_filtered.csv"))

## =========================== fig11: promoter-length-corrected, repeat-filtered (fig8 equivalent) ===========================
bed <- fread(promoter_windows_nonzero_bed, header = FALSE,
             col.names = c("chrom", "start0", "end", "gene", "score", "strand"))
gene_width <- bed[, .(promoter_bp = sum(end - start0)), by = gene]
d11 <- merge(gene_width, d10, by = "gene")

genome_bp <- 3.05e9
rate_nr <- length(gr_itv_nr) / genome_bp   # genome-wide repeat-filtered interval rate, same denominator convention as 13
d11[, expected := promoter_bp * rate_nr]
d11[, poisson_fold := n_hits / pmax(expected, 1e-9)]
cat(sprintf("fig11: repeat-filtered, size-corrected fold: min=%.2f, median=%.2f, max=%.2f\n",
            min(d11$poisson_fold), median(d11$poisson_fold), max(d11$poisson_fold)))

p11 <- make_marked_histogram(d11, "poisson_fold", binwidth = 0.1,
                              x_lab = "repeat-filtered FIMO hits per gene, corrected for promoter length\n(observed / expected given that gene's total promoter bp; 1.0 = size alone predicts this)",
                              title = "Promoter-length-corrected FIMO motif hits per gene, repeat elements removed")
ggsave(file.path(out_figures, "fig11_hits_per_gene_size_corrected_repeat_filtered.png"), p11, width = 8, height = 5.5, dpi = 150, bg = "white")
fwrite(d11[order(-poisson_fold), .(gene, promoter_bp, n_hits, poisson_fold)],
       file.path(out_results, "hits_per_gene_size_corrected_repeat_filtered.csv"))

## =========================== fig12: fixed -3000/+1000 window, repeat-filtered (fig9 equivalent) ===========================
fixed_bed <- file.path(out_data, "fixed_window_promoters_nonzero.bed")
stopifnot(file.exists(fixed_bed))
promo_fixed <- read_promoter_bed(fixed_bed)

d12 <- count_hits(gr_itv_nr, promo_fixed)
cat(sprintf("fig12: repeat-filtered hits per gene (fixed window): n=%s, min=%d, median=%.0f, mean=%.1f, max=%d\n",
            format(nrow(d12), big.mark = ","), min(d12$n_hits), median(d12$n_hits), mean(d12$n_hits), max(d12$n_hits)))

p12 <- make_marked_histogram(d12, "n_hits", binwidth = 1,
                              x_lab = "repeat-filtered FIMO motif hits per gene (single -3000/+1000-of-TSS window, 30 motifs combined)",
                              title = "FIMO motif hits per gene, fixed promoter window, repeat elements removed")
ggsave(file.path(out_figures, "fig12_hits_per_gene_fixed_window_repeat_filtered.png"), p12, width = 8, height = 5.5, dpi = 150, bg = "white")
fwrite(d12[order(-n_hits)], file.path(out_results, "hits_per_gene_fixed_window_repeat_filtered.csv"))

## ---- before/after check on the previously extreme genes ----
extreme_genes <- c("ABR", "TBCD", "OBSCN", "RBFOX3", "EHMT1", marker_genes)
before_raw <- fread(file.path(out_results, "hits_per_gene_combined_baseline.csv"))
setnames(before_raw, "n_hits", "n_hits_before")
check <- merge(before_raw[gene %in% extreme_genes], d10[gene %in% extreme_genes, .(gene, n_hits_after = n_hits)], by = "gene")
cat("\n=== Before vs. after repeat-filtering, raw hit counts (variable window) ===\n")
print(check[order(-n_hits_before)])

message("\nWrote fig10-fig12 and their CSVs")
