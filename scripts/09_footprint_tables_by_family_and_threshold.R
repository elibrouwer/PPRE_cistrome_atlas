# PPRE motif genomic footprint for the 42 human-source PPAR/RXR motifs.
#
# Same definitions as PPRE_threshold_sensitivity/04_dedup_intervals_by_threshold.R and
# 05_promoter_overlap_gene_ranking.R (30-motif run), applied to FIMO_results_human42_unmasked:
#   * hits on the standard chromosomes (chr1-22, X, Y, M), merged with GenomicRanges::reduce (strand ignored)
#   * footprint % = merged length / 3,088,286,401 bp (chr1-22, X, Y, M of GRCh38)
#   * gene universe = cardiomyocyte-expressed genes (GTEx v9 snRNA-seq, non-zero) that have a promoter window
#     (-3000/+1000 of the TSS, union over isoforms): 15,703 genes, promoter_windows_nonzero_expressed.bed
#   * genes with >= 1 hit = genes whose promoter window overlaps at least one merged interval
#   * heterodimer / PPAR half-site / RXR half-site rows: merged intervals (merged over all motifs) that contain
#     at least one hit of that family; the three rows are not exclusive
#   * p <= 1e-5 and p <= 1e-6 are cumulative subsets of the p <= 1e-4 hits
# Outputs: pooled table (per tier and family) and per-motif table (per tier).

suppressPackageStartupMessages({ library(data.table); library(GenomicRanges) })

base_dir <- "C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Minor internship/PPAR_RXR_motifs_Hocomoco_Jaspar"
fimo_dir <- file.path(base_dir, "FIMO_results_human42_unmasked")
batches  <- c("1_4", "5_8", "9_12", "13_16", "17_20", "21_24", "25_28", "29_32", "33_36", "37_40", "41_42")
inv_csv  <- file.path(base_dir, "Motifs_update_2026-09-30", "motif_inventory_human_only.csv")
promoter_bed <- file.path(base_dir, "PPRE_threshold_sensitivity", "data", "promoter_windows_nonzero_expressed.bed")
out_dir  <- "C:/Users/brouw/PPRE_cistrome_atlas/results/human42_pipeline"
genome_bp <- 3088286401
tiers <- c(1e-4, 1e-5, 1e-6)
std <- c(as.character(1:22), "X", "Y", "MT", "M")

inv <- read.csv(inv_csv, check.names = FALSE, fileEncoding = "UTF-8-BOM", stringsAsFactors = FALSE)
family <- setNames(ifelse(inv$category == "PPAR:RXR heterodimer", "heterodimer",
                   ifelse(grepl("^PPAR", toupper(inv$name)), "PPAR_half", "RXR_half")), inv$motif_id)
stopifnot(length(family) == 42)
print(table(family))

hits <- rbindlist(lapply(batches, function(b) {
  dt <- fread(file.path(fimo_dir, b, "fimo.tsv"), sep = "\t", header = TRUE, fill = TRUE)
  dt <- dt[!is.na(start) & !startsWith(motif_id, "#")]
  dt <- dt[, .(motif_id, chr = as.character(sequence_name), start = as.integer(start), stop = as.integer(stop),
               pvalue = as.numeric(`p-value`))]
  dt <- dt[chr %in% std]
  dt[, chr := ifelse(chr %in% c("MT", "M"), "chrM", paste0("chr", chr))]
  cat(b, format(nrow(dt), big.mark = ","), "hits on standard chromosomes\n"); flush.console()
  dt
}))
hits[, family := family[motif_id]]
stopifnot(!anyNA(hits$family), uniqueN(hits$motif_id) == 42)
cat("total hits (standard chromosomes):", format(nrow(hits), big.mark = ","), "\n")

pdf <- fread(promoter_bed, header = FALSE, col.names = c("chrom", "start0", "end", "gene", "score", "strand"))
prom <- GRanges(pdf$chrom, IRanges(pdf$start0 + 1L, pdf$end), strand = pdf$strand, gene = pdf$gene)
n_universe <- uniqueN(pdf$gene)
cat("gene universe:", n_universe, "\n")

genes_hit <- function(merged) uniqueN(prom$gene[subjectHits(findOverlaps(merged, prom, ignore.strand = TRUE))])

pooled <- list()
for (pt in tiers) {
  sub <- hits[pvalue <= pt]
  gr <- GRanges(sub$chr, IRanges(sub$start, sub$stop))
  merged <- reduce(gr, ignore.strand = TRUE)
  iv <- subjectHits(findOverlaps(gr, merged, ignore.strand = TRUE))
  rows <- list(data.table(view = "all", bp = sum(width(merged)), n_intervals = length(merged), genes = genes_hit(merged)))
  for (f in c("heterodimer", "PPAR_half", "RXR_half")) {
    has <- logical(length(merged)); has[unique(iv[sub$family == f])] <- TRUE
    rows[[length(rows) + 1]] <- data.table(view = f, bp = sum(width(merged[has])), n_intervals = sum(has), genes = genes_hit(merged[has]))
  }
  res <- rbindlist(rows); res[, `:=`(threshold = pt, n_hits = nrow(sub))]
  pooled[[length(pooled) + 1]] <- res
  cat("tier", pt, ": merged intervals", format(length(merged), big.mark = ","), "\n"); flush.console()
  rm(gr, merged, iv); gc()
}
pooled <- rbindlist(pooled)
pooled[, `:=`(footprint_pct = 100 * bp / genome_bp, gene_coverage_pct = 100 * genes / n_universe, genome_bp = genome_bp, gene_universe = n_universe)]
fwrite(pooled, file.path(out_dir, "footprint_pooled_by_family_and_threshold.csv"))
print(pooled)

per_motif <- list()
for (m in inv$motif_id) for (pt in tiers) {
  s <- hits[motif_id == m & pvalue <= pt]
  if (nrow(s) == 0) { per_motif[[length(per_motif) + 1]] <- data.table(motif_id = m, threshold = pt, n_hits = 0L, bp = 0L, genes = 0L); next }
  merged <- reduce(GRanges(s$chr, IRanges(s$start, s$stop)), ignore.strand = TRUE)
  per_motif[[length(per_motif) + 1]] <- data.table(motif_id = m, threshold = pt, n_hits = nrow(s), bp = sum(width(merged)), genes = genes_hit(merged))
}
per_motif <- rbindlist(per_motif)
per_motif[, `:=`(family = family[motif_id], footprint_pct = 100 * bp / genome_bp, gene_coverage_pct = 100 * genes / n_universe)]
fwrite(per_motif, file.path(out_dir, "footprint_per_motif_by_threshold.csv"))
cat("done\n")
