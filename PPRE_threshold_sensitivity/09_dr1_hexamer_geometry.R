# Orthogonal, PWM-independent evidence: scan actual promoter DNA directly for
# the classical nuclear-receptor half-site hexamer (AGGTCA-like) and look for
# two such hexamers on the SAME strand separated by the canonical PPAR::RXR
# DR1 (direct repeat, 1bp spacer) architecture.
#
# Why not just measure the gap between our existing PPAR_half/RXR_half FIMO
# hits instead? Checked motif widths directly from the .meme files first:
# every one of the 27 half-site motifs in this project is 9-23bp wide (none
# are pure 6bp hexamers) -- HOCOMOCO/JASPAR/CIS-BP monomer PWMs are built from
# ChIP-seq data that already captures flanking context, not an isolated
# hexamer. So "gap between reported FIMO hit boundaries" would not correspond
# to a real spacer measurement. Scanning for the canonical hexamer directly in
# the genomic sequence sidesteps that ambiguity entirely and is the standard
# way direct-repeat response elements are identified independent of any
# particular PWM (matches max.mismatch=1 to allow common natural variants
# like AGGACA/AGTTCA, not just the textbook AGGTCA).
#
# This is intentionally a DIFFERENT, complementary line of evidence from the
# FIMO-threshold analysis (01-08) -- it doesn't use FIMO output at all, and
# is not affected by the p-value threshold question. It exists to test
# whether "half-site-only" genes from the family-contribution analysis have
# real DR1 geometry backing them, or whether that label is closer to noise.

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(Biostrings)
  library(BSgenome.Hsapiens.UCSC.hg38)
})
source("./PPRE_threshold_sensitivity/00_config.R")

genome <- BSgenome.Hsapiens.UCSC.hg38
HEXAMER <- DNAString("AGGTCA")
MAX_MISMATCH <- 1L

read_promoter_bed <- function(bed) {
  df <- fread(bed, header = FALSE,
              col.names = c("chrom", "start0", "end", "gene", "score", "strand"))
  GRanges(seqnames = df$chrom, ranges = IRanges(start = df$start0 + 1L, end = df$end),
          strand = "*", gene = df$gene)   # scan both strands independently regardless of gene strand
}
promo <- read_promoter_bed(promoter_windows_nonzero_bed)
promo <- promo[width(promo) >= 13]   # can't fit AGGTCA-N-AGGTCA (13bp) in anything shorter
seqlengths(promo) <- seqlengths(genome)[seqlevels(promo)]
promo <- trim(promo)

cat(sprintf("Scanning %d promoter windows (%d genes) for AGGTCA-like hexamers (<=%d mismatch)...\n",
            length(promo), length(unique(promo$gene)), MAX_MISMATCH))

seqs_fwd <- getSeq(genome, promo)
seqs_rev <- reverseComplement(seqs_fwd)

find_dr1_pairs <- function(seqs, strand_label) {
  hits <- vmatchPattern(HEXAMER, seqs, max.mismatch = MAX_MISMATCH)
  n_win <- length(seqs)
  out <- vector("list", n_win)
  for (i in seq_len(n_win)) {
    h <- hits[[i]]
    if (length(h) < 2) next
    st <- start(h); en <- end(h)
    ord <- order(st)
    st <- st[ord]; en <- en[ord]
    gaps <- st[-1] - en[-length(en)] - 1L
    if (length(gaps) == 0) next
    out[[i]] <- data.table(
      window_idx = i, strand = strand_label,
      n_hexamers = length(h),
      n_DR0 = sum(gaps == 0), n_DR1 = sum(gaps == 1), n_DR2 = sum(gaps == 2),
      min_gap = min(gaps)
    )
  }
  rbindlist(out)
}

res_fwd <- find_dr1_pairs(seqs_fwd, "+")
res_rev <- find_dr1_pairs(seqs_rev, "-")
res <- rbindlist(list(res_fwd, res_rev))
res[, gene := promo$gene[window_idx]]

## ---- gap-distance distribution across ALL adjacent hexamer pairs (sanity check: ----
## a real biological DR1 preference should show a local peak at gap==1 rather
## than a flat/monotonic distribution driven purely by hexamer combinatorics) ----
all_gaps <- rbindlist(list(
  { h <- vmatchPattern(HEXAMER, seqs_fwd, max.mismatch = MAX_MISMATCH)
    rbindlist(lapply(h, function(x) if (length(x) >= 2) {
      st <- sort(start(x)); en <- sort(end(x))
      data.table(gap = st[-1] - en[-length(en)] - 1L)
    } else NULL)) },
  { h <- vmatchPattern(HEXAMER, seqs_rev, max.mismatch = MAX_MISMATCH)
    rbindlist(lapply(h, function(x) if (length(x) >= 2) {
      st <- sort(start(x)); en <- sort(end(x))
      data.table(gap = st[-1] - en[-length(en)] - 1L)
    } else NULL)) }
))
all_gaps <- all_gaps[gap >= 0 & gap <= 20]
gap_dist <- all_gaps[, .N, by = gap][order(gap)]
fwrite(gap_dist, file.path(out_results, "dr1_hexamer_gap_distribution.csv"))
cat("\nGap-distance distribution (0-20bp), all adjacent hexamer pairs genome-wide-in-promoters:\n")
print(gap_dist)

## ---- per-gene summary: does this gene's promoter contain >=1 genuine DR1 pair? ----
gene_dr1 <- res[, .(
  n_DR1_pairs = sum(n_DR1), n_DR0_pairs = sum(n_DR0), n_DR2_pairs = sum(n_DR2),
  n_hexamer_windows = .N
), by = gene]
gene_dr1[, has_DR1 := n_DR1_pairs > 0]

# genes with a promoter window but zero hexamer pairs at all (not in `res`)
all_genes <- unique(promo$gene)
missing <- setdiff(all_genes, gene_dr1$gene)
if (length(missing) > 0) {
  gene_dr1 <- rbindlist(list(gene_dr1, data.table(
    gene = missing, n_DR1_pairs = 0L, n_DR0_pairs = 0L, n_DR2_pairs = 0L,
    n_hexamer_windows = 0L, has_DR1 = FALSE
  )))
}

fwrite(gene_dr1, file.path(out_results, "dr1_hexamer_geometry_per_gene.csv"))
cat(sprintf("\n%d / %d expressed genes (%.1f%%) have >=1 independent AGGTCA-DR1 hexamer pair in their promoter (PWM-independent evidence).\n",
            sum(gene_dr1$has_DR1), nrow(gene_dr1), 100 * mean(gene_dr1$has_DR1)))

## ---- cross-check: does DR1 hexamer evidence track with the FIMO heterodimer-PWM
## evidence, or with half-site-only status, from the threshold analysis? This
## validates (or undercuts) whether "half-site-only at strict threshold" is
## biologically hollow or still geometrically plausible. ----
fam_gene_scores <- fread(file.path(out_results, "gene_scores_by_threshold.csv"))
het_genes_p4 <- fam_gene_scores[threshold == 1e-4 & view == "heterodimer" & gene_set == "nonzero", gene]
het_genes_p6 <- fam_gene_scores[threshold == 1e-6 & view == "heterodimer" & gene_set == "nonzero", gene]
comb_genes_p6 <- fam_gene_scores[threshold == 1e-6 & view == "combined"    & gene_set == "nonzero", gene]
half_only_p6 <- setdiff(comb_genes_p6, het_genes_p6)

cat(sprintf("\nDR1 hexamer rate, heterodimer-PWM-supported genes (p<=1e-4): %.1f%% (n=%d)\n",
            100 * mean(gene_dr1$has_DR1[gene_dr1$gene %in% het_genes_p4]), length(het_genes_p4)))
cat(sprintf("DR1 hexamer rate, heterodimer-PWM-supported genes (p<=1e-6): %.1f%% (n=%d)\n",
            100 * mean(gene_dr1$has_DR1[gene_dr1$gene %in% het_genes_p6]), length(het_genes_p6)))
cat(sprintf("DR1 hexamer rate, half-site-ONLY genes at p<=1e-6 (no heterodimer PWM hit): %.1f%% (n=%d)\n",
            100 * mean(gene_dr1$has_DR1[gene_dr1$gene %in% half_only_p6]), length(half_only_p6)))
cat(sprintf("DR1 hexamer rate, ALL nonzero-expressed genes (background rate): %.1f%% (n=%d)\n",
            100 * mean(gene_dr1$has_DR1), nrow(gene_dr1)))

message("\nDone. Written: dr1_hexamer_geometry_per_gene.csv, dr1_hexamer_gap_distribution.csv")
