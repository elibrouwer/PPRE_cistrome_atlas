# Step 3: for each gene's promoter window, compute what fraction of FIMO motif hits
# fall inside a real PPARA ChIP-seq peak, using two independent datasets:
#   A) GSM4748812 - single HepG2 PPARA ChIP-seq experiment (MACS2 narrowPeak, hg38)
#   B) ChIP-Atlas "Oth.ALL.05.PPARA.AllCell" - PPARA peaks aggregated across every public
#      PPARA ChIP-seq experiment/cell type (hg38); reduced to non-redundant regions with a
#      per-region "supporting experiments" count, so we can ask for recurrent (not one-off)
#      binding.
# Also computes the reverse (recall) direction: of the real ChIP peaks in each promoter,
# how many actually contain a FIMO hit for any of the 30 motifs.

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
})

out_dir <- "PPARA_chipseq_overlap"

## ---- load promoter-restricted FIMO hits --------------------------------------------
fimo_cols <- c("motif_id","motif_alt_id","sequence_name","start","stop","strand",
               "score","p_value","q_value","matched_sequence","gene")
fimo <- fread(file.path(out_dir, "promoter_fimo_hits.tsv"), header = FALSE, col.names = fimo_cols)
fimo[, chrom := paste0("chr", sequence_name)]
fimo_gr <- GRanges(fimo$chrom, IRanges(fimo$start, fimo$stop), strand = fimo$strand,
                    motif_id = fimo$motif_id, gene = fimo$gene, score = fimo$score)

## ---- load dataset A: single-experiment HepG2 PPARA ChIP-seq ------------------------
gsm_cols <- c("chrom","start","end","name","score","strand","signalValue","pValue","qValue","summit","gene")
gsm_path <- file.path(out_dir, "promoter_GSM4748812_peaks.bed")
gsm <- if (file.size(gsm_path) > 0) {
  fread(gsm_path, header = FALSE, sep = "\t", quote = "", col.names = gsm_cols)
} else data.table(matrix(nrow = 0, ncol = length(gsm_cols), dimnames = list(NULL, gsm_cols)))
gsm_gr <- if (nrow(gsm) > 0) GRanges(gsm$chrom, IRanges(gsm$start + 1, gsm$end)) else GRanges()

## ---- load dataset B: ChIP-Atlas all-cell-type PPARA aggregate ----------------------
ca_cols <- c("chrom","start","end","info","score","strand","thickStart","thickEnd","itemRgb","gene")
ca_path <- file.path(out_dir, "promoter_chipatlas_peaks.bed")
ca <- if (file.size(ca_path) > 0) {
  fread(ca_path, header = FALSE, sep = "\t", quote = "", col.names = ca_cols, fill = TRUE)
} else data.table(matrix(nrow = 0, ncol = length(ca_cols), dimnames = list(NULL, ca_cols)))
ca_gr_raw <- if (nrow(ca) > 0) GRanges(ca$chrom, IRanges(ca$start + 1, ca$end)) else GRanges()
# collapse to non-redundant loci + count how many independent experiment-peaks support each
ca_reduced <- reduce(ca_gr_raw)
if (length(ca_reduced) > 0) {
  ca_reduced$n_experiments <- countOverlaps(ca_reduced, ca_gr_raw)
}
ca_reduced_recurrent <- ca_reduced[ca_reduced$n_experiments >= 3]

## ---- overlap FIMO hits against each ChIP dataset ------------------------------------
fimo$in_GSM4748812     <- overlapsAny(fimo_gr, gsm_gr)
fimo$in_ChIPAtlas_any  <- overlapsAny(fimo_gr, ca_reduced)
fimo$in_ChIPAtlas_ge3  <- overlapsAny(fimo_gr, ca_reduced_recurrent)

## ---- per-gene summary: precision (FIMO hit -> real ChIP peak?) ----------------------
genes <- sort(unique(fimo$gene))
summary_rows <- lapply(genes, function(g) {
  sub <- fimo[gene == g]
  data.table(
    gene = g,
    n_fimo_hits = nrow(sub),
    n_in_GSM4748812 = sum(sub$in_GSM4748812),
    pct_in_GSM4748812 = round(100 * mean(sub$in_GSM4748812), 1),
    n_in_ChIPAtlas_any = sum(sub$in_ChIPAtlas_any),
    pct_in_ChIPAtlas_any = round(100 * mean(sub$in_ChIPAtlas_any), 1),
    n_in_ChIPAtlas_ge3exp = sum(sub$in_ChIPAtlas_ge3),
    pct_in_ChIPAtlas_ge3exp = round(100 * mean(sub$in_ChIPAtlas_ge3), 1)
  )
})
precision_summary <- rbindlist(summary_rows)

## ---- per-gene reverse (recall): of real ChIP peaks in the promoter, how many contain
## a FIMO hit for ANY of the 30 motifs? -------------------------------------------------
recall_rows <- lapply(genes, function(g) {
  gsm_g <- gsm_gr[gsm$gene == g]
  gene_window <- range(fimo_gr[fimo_gr$gene == g], ignore.strand = TRUE)
  ca_g <- subsetByOverlaps(ca_reduced, gene_window)
  fimo_g_gr <- fimo_gr[fimo_gr$gene == g]
  data.table(
    gene = g,
    n_GSM4748812_peaks = length(gsm_g),
    n_GSM4748812_with_fimo = sum(overlapsAny(gsm_g, fimo_g_gr)),
    n_ChIPAtlas_loci = length(ca_g),
    n_ChIPAtlas_loci_with_fimo = sum(overlapsAny(ca_g, fimo_g_gr))
  )
})
recall_summary <- rbindlist(recall_rows)

cat("\n================= PRECISION: FIMO hit -> real ChIP-seq peak? =================\n")
print(precision_summary)
cat("\n================= RECALL: real ChIP-seq peak -> any FIMO hit? =================\n")
print(recall_summary)

fwrite(precision_summary, file.path(out_dir, "precision_summary_fimo_to_chip.csv"))
fwrite(recall_summary, file.path(out_dir, "recall_summary_chip_to_fimo.csv"))
fwrite(fimo, file.path(out_dir, "promoter_fimo_hits_annotated.csv"))

## ---- per-gene, per-motif breakdown (which PWMs are enriching for real binding) ------
motif_rows <- fimo[, .(
  n_hits = .N,
  n_in_GSM4748812 = sum(in_GSM4748812),
  pct_in_GSM4748812 = round(100 * mean(in_GSM4748812), 1),
  n_in_ChIPAtlas_any = sum(in_ChIPAtlas_any),
  pct_in_ChIPAtlas_any = round(100 * mean(in_ChIPAtlas_any), 1)
), by = .(motif_id)][order(-pct_in_ChIPAtlas_any)]
fwrite(motif_rows, file.path(out_dir, "per_motif_overlap_summary.csv"))
cat("\n================= PER-MOTIF (all genes combined) =================\n")
print(motif_rows)

cat("\nAll summaries written to", out_dir, "\n")
