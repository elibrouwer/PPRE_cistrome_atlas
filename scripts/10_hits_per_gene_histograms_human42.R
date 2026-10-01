# Distribution of FIMO motif hits per gene for the 42 human-source PPAR/RXR motifs (fig7 to fig12).
#
# Re-implementation of PPRE_threshold_sensitivity/17_hits_per_gene_histogram.R, 18_corrected_and_fixed_window_histograms.R
# and 19_repeat_filtered_histograms.R (30-motif run) for FIMO_results_human42_unmasked, with the same definitions:
#   * hits = deduplicated (merged, strand ignored) FIMO intervals at p <= 1e-4 on the standard chromosomes
#   * hits per gene = number of those intervals overlapping the gene's promoter windows (-3000/+1000 of the TSS, union over
#     isoforms), for the cardiomyocyte-expressed gene universe (GTEx v9 non-zero, 15,703 genes)
#   fig7  raw hits per gene
#   fig8  promoter-length corrected: observed / (promoter bp x genome-wide interval rate), genome = 3.05e9 bp
#   fig9  fixed single -3000/+1000 window (widest transcript per gene)
#   fig10-12  the same three after dropping every interval that overlaps a UCSC RepeatMasker element
# The promoter windows, the fixed-window file and the RepeatMasker mask are motif independent and are reused from
# PPRE_threshold_sensitivity/data.

suppressPackageStartupMessages({ library(data.table); library(GenomicRanges); library(ggplot2); library(scales) })

base_dir <- "C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Minor internship/PPAR_RXR_motifs_Hocomoco_Jaspar"
fimo_dir <- file.path(base_dir, "FIMO_results_human42_unmasked")
batches  <- c("1_4", "5_8", "9_12", "13_16", "17_20", "21_24", "25_28", "29_32", "33_36", "37_40", "41_42")
data_dir <- file.path(base_dir, "PPRE_threshold_sensitivity", "data")
out_dir  <- "C:/Users/brouw/PPRE_cistrome_atlas/results/human42_pipeline/hits_per_gene"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
std <- c(as.character(1:22), "X", "Y", "MT", "M")
marker_genes <- c("CPT1A", "CPT1B", "HADHA", "HADHB", "ACADVL")
genome_bp <- 3.05e9
n_motifs <- 42

## ---- deduplicated intervals, p <= 1e-4 ----
hits <- rbindlist(lapply(batches, function(b) {
  dt <- fread(file.path(fimo_dir, b, "fimo.tsv"), sep = "\t", header = TRUE, fill = TRUE)
  dt <- dt[!is.na(start) & !startsWith(motif_id, "#")]
  dt <- dt[as.character(sequence_name) %in% std, .(chr = as.character(sequence_name), start = as.integer(start), stop = as.integer(stop))]
  dt[, chr := ifelse(chr %in% c("MT", "M"), "chrM", paste0("chr", chr))]
  cat(b, "loaded\n"); flush.console(); dt
}))
merged <- reduce(GRanges(hits$chr, IRanges(hits$start, hits$stop)), ignore.strand = TRUE)
rm(hits); gc()
n_total_intervals <- length(merged)
cat(sprintf("deduplicated intervals at p<=1e-4: %s\n", format(n_total_intervals, big.mark = ",")))

## ---- shared inputs ----
read_bed <- function(bed) {
  df <- fread(bed, header = FALSE, col.names = c("chrom", "start0", "end", "gene", "score", "strand"))
  list(gr = GRanges(df$chrom, IRanges(df$start0 + 1L, df$end), gene = df$gene), df = df)
}
promo_var <- read_bed(file.path(data_dir, "promoter_windows_nonzero_expressed.bed"))
promo_fix <- read_bed(file.path(data_dir, "fixed_window_promoters_nonzero.bed"))
repeat_gr <- readRDS(file.path(data_dir, "rmsk_hg38_merged_gr.rds"))
gene_bp <- promo_var$df[, .(promoter_bp = sum(end - start0)), by = gene]

count_hits <- function(itv, promo_gr) {
  ov <- findOverlaps(itv, promo_gr, ignore.strand = TRUE)
  dt <- data.table(gene = promo_gr$gene[subjectHits(ov)], interval_idx = queryHits(ov))[, .(n_hits = uniqueN(interval_idx)), by = gene]
  dt <- merge(data.table(gene = unique(promo_gr$gene)), dt, by = "gene", all.x = TRUE)
  dt[is.na(n_hits), n_hits := 0L]
  dt[, n_hits := as.integer(n_hits)]
  dt
}

make_marked_histogram <- function(d, value_col, binwidth, x_lab, title, y_lab = "number of genes") {
  markers <- d[gene %in% marker_genes, c("gene", value_col), with = FALSE]
  setnames(markers, value_col, "value")
  if (length(setdiff(marker_genes, markers$gene)) > 0) warning("Marker gene(s) not found: ", paste(setdiff(marker_genes, markers$gene), collapse = ", "))
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
save_fig <- function(p, name) ggsave(file.path(out_dir, name), p, width = 8, height = 5.5, dpi = 150, bg = "white")
summ <- function(label, x) cat(sprintf("%s: n=%s min=%s median=%s mean=%.1f max=%s\n", label, format(length(x), big.mark = ","),
                                       format(min(x), digits = 4), format(median(x), digits = 4), mean(x), format(max(x), digits = 4)))

## ---- fig7: raw hits per gene ----
d7 <- count_hits(merged, promo_var$gr); d7 <- d7[n_hits > 0]
summ("fig7 hits per gene", d7$n_hits)
save_fig(make_marked_histogram(d7, "n_hits", 5, sprintf("FIMO motif hits per gene (promoter-overlapping intervals, %d motifs combined)", n_motifs),
                               "Distribution of FIMO motif hits per gene"), "fig7_hits_per_gene_histogram.png")
fwrite(d7[order(-n_hits)], file.path(out_dir, "hits_per_gene_combined_baseline.csv"))

## ---- fig8: promoter-length corrected ----
rate <- n_total_intervals / genome_bp
d8 <- merge(gene_bp, d7, by = "gene")
d8[, expected := promoter_bp * rate][, fold_corrected := n_hits / pmax(expected, 1e-9)]
summ("fig8 corrected fold", d8$fold_corrected)
save_fig(make_marked_histogram(d8, "fold_corrected", 0.1,
                               "FIMO hits per gene, corrected for promoter length\n(observed / expected given that gene's total promoter bp; 1.0 = size alone predicts this)",
                               "Promoter-length-corrected FIMO motif hits per gene"), "fig8_hits_per_gene_size_corrected.png")
fwrite(d8[order(-fold_corrected), .(gene, promoter_bp, n_hits_raw = n_hits, fold_corrected)], file.path(out_dir, "hits_per_gene_size_corrected.csv"))

## ---- fig9: fixed -3000/+1000 window ----
d9 <- count_hits(merged, promo_fix$gr)
summ("fig9 fixed window", d9$n_hits)
save_fig(make_marked_histogram(d9, "n_hits", 1, sprintf("FIMO motif hits per gene (single -3000/+1000-of-TSS window, %d motifs combined)", n_motifs),
                               "FIMO motif hits per gene, fixed promoter window"), "fig9_hits_per_gene_fixed_window.png")
fwrite(d9[order(-n_hits)], file.path(out_dir, "hits_per_gene_fixed_window.csv"))

## ---- fig10-12: repeat elements removed ----
in_repeat <- overlapsAny(merged, repeat_gr, ignore.strand = TRUE)
cat(sprintf("intervals overlapping a RepeatMasker element: %s of %s (%.1f%%) -- dropped\n",
            format(sum(in_repeat), big.mark = ","), format(length(merged), big.mark = ","), 100 * mean(in_repeat)))
merged_nr <- merged[!in_repeat]

d10 <- count_hits(merged_nr, promo_var$gr)
summ("fig10 repeat-filtered", d10$n_hits)
save_fig(make_marked_histogram(d10, "n_hits", 5, sprintf("repeat-filtered FIMO motif hits per gene (promoter-overlapping intervals, %d motifs combined)", n_motifs),
                               "Distribution of FIMO motif hits per gene, repeat elements removed"), "fig10_hits_per_gene_repeat_filtered.png")
fwrite(d10[order(-n_hits)], file.path(out_dir, "hits_per_gene_repeat_filtered.csv"))

d11 <- merge(gene_bp, d10, by = "gene")
rate_nr <- length(merged_nr) / genome_bp
d11[, expected := promoter_bp * rate_nr][, poisson_fold := n_hits / pmax(expected, 1e-9)]
summ("fig11 repeat-filtered corrected fold", d11$poisson_fold)
save_fig(make_marked_histogram(d11, "poisson_fold", 0.1,
                               "repeat-filtered FIMO hits per gene, corrected for promoter length\n(observed / expected given that gene's total promoter bp; 1.0 = size alone predicts this)",
                               "Promoter-length-corrected FIMO motif hits per gene, repeat elements removed"), "fig11_hits_per_gene_size_corrected_repeat_filtered.png")
fwrite(d11[order(-poisson_fold), .(gene, promoter_bp, n_hits, poisson_fold)], file.path(out_dir, "hits_per_gene_size_corrected_repeat_filtered.csv"))

d12 <- count_hits(merged_nr, promo_fix$gr)
summ("fig12 repeat-filtered fixed window", d12$n_hits)
save_fig(make_marked_histogram(d12, "n_hits", 1, sprintf("repeat-filtered FIMO motif hits per gene (single -3000/+1000-of-TSS window, %d motifs combined)", n_motifs),
                               "FIMO motif hits per gene, fixed promoter window, repeat elements removed"), "fig12_hits_per_gene_fixed_window_repeat_filtered.png")
fwrite(d12[order(-n_hits)], file.path(out_dir, "hits_per_gene_fixed_window_repeat_filtered.csv"))

cat("\nmarker genes (raw, variable window | repeat-filtered | fixed window):\n")
print(merge(merge(d7[gene %in% marker_genes, .(gene, raw = n_hits)], d10[gene %in% marker_genes, .(gene, repeat_filtered = n_hits)], by = "gene"),
            d9[gene %in% marker_genes, .(gene, fixed_window = n_hits)], by = "gene"))
cat("done\n")
