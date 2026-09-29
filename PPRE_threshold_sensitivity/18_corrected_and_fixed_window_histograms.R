# Two follow-ups to fig7 (raw FIMO motif hits per gene), both addressing the
# same size confound documented in 13_size_normalized_scoring.R (genes with
# many alternative-TSS promoter windows accumulate more raw hits by pure
# chance, nothing to do with PPAR biology):
#
#   fig8: promoter-length-CORRECTED hits per gene -- poisson_fold_combined
#         from 13 (observed hits / hits expected from that gene's own total
#         promoter bp at the genome-wide rate). 1.0 = exactly as many hits as
#         promoter size alone predicts.
#   fig9: hits per gene using a single FIXED -3000/+1000-of-TSS window per
#         gene (one canonical transcript, the widest one, instead of 13's
#         union across every annotated isoform) -- so every gene gets the
#         same ~4kb window and the size confound cannot arise in the first
#         place, at the cost of ignoring alternative promoters.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(GenomicRanges)
  library(readr)
})
source("./PPRE_threshold_sensitivity/00_config.R")

marker_genes <- c("CPT1A", "CPT1B", "HADHA", "HADHB", "ACADVL")

# shared plotting helper: histogram + staggered marker-gene labels, same look as fig7
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

## =========================== fig8: promoter-length-corrected ===========================
size_norm <- fread(file.path(out_results, "gene_scoring_size_normalized.csv"))
cat(sprintf("fig8 gene universe: %s genes\n", format(nrow(size_norm), big.mark = ",")))

p8 <- make_marked_histogram(size_norm, "poisson_fold_combined", binwidth = 0.1,
                             x_lab = "FIMO hits per gene, corrected for promoter length\n(observed / expected given that gene's total promoter bp; 1.0 = size alone predicts this)",
                             title = "Promoter-length-corrected FIMO motif hits per gene")
ggsave(file.path(out_figures, "fig8_hits_per_gene_size_corrected.png"), p8, width = 8, height = 5.5, dpi = 150, bg = "white")
fwrite(size_norm[, .(gene, promoter_bp, n_hits_raw = fisher_n_loci_combined, fold_corrected = poisson_fold_combined)][order(-fold_corrected)],
       file.path(out_results, "hits_per_gene_size_corrected.csv"))

## =========================== fig9: fixed -3000/+1000 single-TSS window ===========================
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
tx_by_gene <- transcriptsBy(txdb, by = "gene")

nonzero_genes <- unique(read_csv(expr_nonzero_csv, show_col_types = FALSE)$HGNC_gene_name)
entrez_map <- AnnotationDbi::select(org.Hs.eg.db, keys = nonzero_genes, keytype = "SYMBOL", columns = "ENTREZID")
entrez_map <- entrez_map[!is.na(entrez_map$ENTREZID), ]
entrez_map <- entrez_map[!duplicated(entrez_map$SYMBOL), ]

# one canonical transcript per gene (the widest annotated one), NOT the union
# across all isoforms -- this is what pins every gene's window to a fixed 4kb
promo_list <- vector("list", nrow(entrez_map))
n_no_tx <- 0L
for (i in seq_len(nrow(entrez_map))) {
  eid <- entrez_map$ENTREZID[i]
  if (!(eid %in% names(tx_by_gene))) { n_no_tx <- n_no_tx + 1L; next }
  tx_gr <- tx_by_gene[[eid]]
  canonical <- tx_gr[which.max(width(tx_gr))]
  promo <- promoters(canonical, upstream = 3000, downstream = 1000)
  mcols(promo) <- NULL
  mcols(promo)$gene <- entrez_map$SYMBOL[i]
  promo_list[[i]] <- promo
}
promo_list <- promo_list[!vapply(promo_list, is.null, logical(1))]
promo_gr <- sort(do.call(c, unname(promo_list)))
std_chroms <- paste0("chr", c(1:22, "X", "Y", "M"))
promo_gr <- promo_gr[as.character(seqnames(promo_gr)) %in% std_chroms]
cat(sprintf("fig9 gene universe: %d genes (%d had no txdb transcript, skipped), each a single %dbp window\n",
            length(promo_gr), n_no_tx, unique(width(promo_gr))[1]))

# saved for reuse (e.g. repeat-filtering in 19) without reloading the TxDb/org.Hs.eg.db packages
fixed_window_bed <- data.frame(chrom = as.character(seqnames(promo_gr)), start = start(promo_gr) - 1L,
                                end = end(promo_gr), gene = promo_gr$gene, score = 0,
                                strand = as.character(strand(promo_gr)))
write.table(fixed_window_bed, file.path(out_data, "fixed_window_promoters_nonzero.bed"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

interval_list <- readRDS(file.path(out_data, "deduplicated_intervals_by_threshold.rds"))
tag <- paste0("p", format(baseline_threshold, scientific = TRUE))
itv <- interval_list[[tag]]  # combined view = all deduplicated intervals, no family filter
gr_itv <- GRanges(itv$chrom, IRanges(itv$start, itv$end))

ov <- findOverlaps(gr_itv, promo_gr, ignore.strand = TRUE)
hit_counts <- data.table(gene = promo_gr$gene[subjectHits(ov)], interval_idx = queryHits(ov))[
  , .(n_hits = uniqueN(interval_idx)), by = gene]

fixed_window <- data.table(gene = promo_gr$gene)
fixed_window <- merge(fixed_window, hit_counts, by = "gene", all.x = TRUE)
fixed_window[is.na(n_hits), n_hits := 0L]
cat(sprintf("fig9: hits per gene, fixed window: min=%d, median=%.0f, mean=%.1f, max=%d\n",
            min(fixed_window$n_hits), median(fixed_window$n_hits), mean(fixed_window$n_hits), max(fixed_window$n_hits)))

p9 <- make_marked_histogram(fixed_window, "n_hits", binwidth = 1,
                             x_lab = "FIMO motif hits per gene (single -3000/+1000-of-TSS window, 30 motifs combined)",
                             title = "FIMO motif hits per gene, fixed promoter window")
ggsave(file.path(out_figures, "fig9_hits_per_gene_fixed_window.png"), p9, width = 8, height = 5.5, dpi = 150, bg = "white")
fwrite(fixed_window[order(-n_hits)], file.path(out_results, "hits_per_gene_fixed_window.csv"))

message("\nWrote fig8_hits_per_gene_size_corrected.png, fig9_hits_per_gene_fixed_window.png, and their CSVs")
