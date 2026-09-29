# Compute several gene-level "how strong is the PPRE evidence in this
# promoter" scores, each grounded in an established statistical or biological
# method rather than an ad hoc heuristic, so they can be compared and
# aggregated on equal footing:
#
#   1. sum(-log10 p) / max(-log10 p)  -- already computed in 05 (this
#      project's own pre-existing convention, see chipseeker_fimo_analysis.Rmd
#      "PPRE_score_sum/max/median/average"); kept as the baseline comparator.
#   2. Fisher's combined probability test (Fisher, 1925, "Statistical Methods
#      for Research Workers") -- combines the p-values of all DISTINCT
#      deduplicated intervals in a gene's promoter into one calibrated
#      combined p-value: -2*sum(ln p_i) ~ chi-sq(2k). Rewards genes with
#      several independent supporting loci, not just one strong hit, in a
#      statistically principled way that Fisher's test was specifically
#      designed for (independent evidence, same null hypothesis).
#   3. Stouffer's weighted Z-method (Stouffer et al. 1949; weighted form:
#      Liptak 1958; standard meta-analysis technique, e.g. as used for GWAS
#      meta-analysis) -- like Fisher's, but allows explicit differential
#      weighting of evidence sources. Used here to directly operationalize
#      the threshold-sensitivity analysis's own finding (01-08) that
#      heterodimer motif hits are the more specific evidence: heterodimer
#      intervals get weight 3, half-site intervals weight 1.
#   4. DR1 hexamer pair count (script 09) -- PWM-independent, grounded in the
#      classical nuclear-receptor direct-repeat spacing rule (Umesono et al.
#      1991, Cell; Naar et al. 1991, Cell): a real DR1 element is two
#      AGGTCA-like hexamers on the same strand with a 1bp spacer, independent
#      of any particular PWM's threshold.
#
# All four are computed at this analysis's recommended primary threshold
# (p<=1e-4, from 00_config.R baseline_threshold) since 01-08 already showed
# tightening the p-value threshold does not improve real ChIP-seq tracking --
# so evidence STRUCTURE (multiple loci, family composition, DR1 geometry) is
# used to differentiate genes here instead of a stricter cutoff.

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
})
source("./PPRE_threshold_sensitivity/00_config.R")

interval_list <- readRDS(file.path(out_data, "deduplicated_intervals_by_threshold.rds"))
tag <- paste0("p", format(baseline_threshold, scientific = TRUE))
itv <- interval_list[[tag]]
stopifnot(!is.null(itv))

read_promoter_bed <- function(bed) {
  df <- fread(bed, header = FALSE,
              col.names = c("chrom", "start0", "end", "gene", "score", "strand"))
  GRanges(seqnames = df$chrom, ranges = IRanges(start = df$start0 + 1L, end = df$end), gene = df$gene)
}
gs <- read_promoter_bed(promoter_windows_nonzero_bed)

view_filter <- function(families_present, view) {
  if (view == "combined") return(rep(TRUE, length(families_present)))
  grepl(view, families_present, fixed = TRUE)
}

## ---- gather, per gene, the list of DISTINCT interval p-values (and family) in its promoter ----
itv_all <- itv[view_filter(itv$families_present, "combined")]
gr_all <- GRanges(seqnames = itv_all$chrom, ranges = IRanges(start = itv_all$start, end = itv_all$end))
ov <- findOverlaps(gr_all, gs, ignore.strand = TRUE)

interval_gene_tbl <- data.table(
  interval_row = queryHits(ov),
  gene = gs$gene[subjectHits(ov)],
  pvalue = itv_all$best_pvalue[queryHits(ov)],
  families_present = itv_all$families_present[queryHits(ov)]
)
# one row per (gene, distinct interval) -- an interval overlapping multiple
# promoter windows of the SAME gene (isoform unions) should count once
interval_gene_tbl <- unique(interval_gene_tbl, by = c("gene", "interval_row"))
interval_gene_tbl[, is_heterodimer := grepl("heterodimer", families_present, fixed = TRUE)]

## ---- Fisher's combined probability test, per gene, per family view ----
# IMPORTANT: with a median of ~47 supporting loci per gene, the linear-scale
# combined p-value underflows to exactly 0 (double-precision floor) for the
# large majority of genes -- pchisq(..., lower.tail=FALSE) alone would tie
# almost everyone at p=0 and destroy all rank information exactly where it's
# needed. Computed on the log scale (log.p=TRUE) instead, which stays a
# genuine, rankable value however small the true p-value is.
fisher_combine <- function(p) {
  p <- pmax(p, 1e-300)  # avoid log(0) going into the statistic itself
  k <- length(p)
  if (k == 0) return(c(neglog10p_combined = 0, k = 0))
  stat <- -2 * sum(log(p))
  log_p <- pchisq(stat, df = 2 * k, lower.tail = FALSE, log.p = TRUE)  # natural log of the upper-tail p
  c(neglog10p_combined = -log_p / log(10), k = k)
}

fisher_by_view <- function(view) {
  sub <- if (view == "combined") interval_gene_tbl else interval_gene_tbl[grepl(view, families_present, fixed = TRUE)]
  sub[, {
    r <- fisher_combine(pvalue)
    list(neglog10p_combined = r["neglog10p_combined"], n_loci = r["k"])
  }, by = gene]
}
fisher_combined  <- fisher_by_view("combined");   setnames(fisher_combined,  c("neglog10p_combined","n_loci"), c("fisher_neglog10p_combined","fisher_n_loci_combined"))
fisher_het       <- fisher_by_view("heterodimer");setnames(fisher_het,       c("neglog10p_combined","n_loci"), c("fisher_neglog10p_heterodimer","fisher_n_loci_heterodimer"))
fisher_ppar_half <- fisher_by_view("PPAR_half");  setnames(fisher_ppar_half, c("neglog10p_combined","n_loci"), c("fisher_neglog10p_PPAR_half","fisher_n_loci_PPAR_half"))
fisher_rxr_half  <- fisher_by_view("RXR_half");   setnames(fisher_rxr_half,  c("neglog10p_combined","n_loci"), c("fisher_neglog10p_RXR_half","fisher_n_loci_RXR_half"))

## ---- Stouffer's weighted Z-method: heterodimer intervals weighted 3x half-site intervals ----
# NOTE: qnorm(1 - p) is NOT safe here -- for any p below ~1e-16, `1 - p`
# rounds to exactly 1.0 in double precision (machine epsilon), silently
# turning a real, finite Z-score into qnorm(1) = Inf. FIMO p-values from wide,
# high-information PWMs in this motif set legitimately go below that. Using
# qnorm(p, lower.tail=FALSE) computes the upper-tail quantile directly,
# without the catastrophic cancellation.
HET_WEIGHT <- 3; HALF_WEIGHT <- 1
stouffer_tbl <- interval_gene_tbl[, {
  p <- pmax(pmin(pvalue, 1 - 1e-16), 1e-300)
  z <- qnorm(p, lower.tail = FALSE)
  w <- ifelse(is_heterodimer, HET_WEIGHT, HALF_WEIGHT)
  z_combined <- sum(w * z) / sqrt(sum(w^2))
  list(stouffer_z = z_combined, stouffer_neglog10p = -pnorm(z_combined, lower.tail = FALSE, log.p = TRUE) / log(10))
}, by = gene]

## ---- DR1 hexamer evidence (script 09) ----
dr1 <- fread(file.path(out_results, "dr1_hexamer_geometry_per_gene.csv"))

## ---- baseline: this project's own pre-existing PPRE_score convention (sum/max -log10 p) ----
baseline_scores <- fread(file.path(out_results, "gene_scores_by_threshold.csv"))[
  threshold == baseline_threshold & gene_set == "nonzero" & view == "combined",
  .(gene, sum_neglog10p, max_neglog10p, n_intervals)]

## ---- assemble one wide table, all genes in the nonzero-expressed universe ----
all_genes <- data.table(gene = unique(gs$gene))
scores <- Reduce(function(a, b) merge(a, b, by = "gene", all.x = TRUE),
                  list(all_genes, baseline_scores, fisher_combined, fisher_het,
                       fisher_ppar_half, fisher_rxr_half, stouffer_tbl,
                       dr1[, .(gene, n_DR1_pairs, has_DR1)]))

fill_cols <- setdiff(names(scores), "gene")
for (col in fill_cols) scores[is.na(get(col)), (col) := 0]   # everything is now a neglog10p-style scale: no evidence = 0
scores[is.na(has_DR1), has_DR1 := FALSE]

fwrite(scores, file.path(out_results, "gene_scoring_methods.csv"))
cat(sprintf("Wrote gene_scoring_methods.csv: %d genes x %d score columns\n", nrow(scores), ncol(scores) - 1))
cat("\nColumn summary:\n"); print(summary(scores[, -"gene"]))
message("Done.")
