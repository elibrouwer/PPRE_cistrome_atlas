# Overlap deduplicated intervals (from 04) with cardiomyocyte-expressed-gene
# promoter windows (from 02), for two gene universes (nonzero-expressed,
# top20-expressed) and four motif views (heterodimer, PPAR_half, RXR_half,
# and "combined" = any family), at each threshold tier.
#
# For every gene x threshold x view: number of supporting intervals, a
# continuous PPRE score (sum/max of -log10(best_pvalue) across its
# promoter-overlapping intervals), and its rank. Rank stability vs. the
# current p<=1e-4 analysis is then reported as Spearman rho and top-100/
# top-500 overlap, per view and per gene universe, separately (not pooled)
# so a family-specific shift isn't hidden by averaging across families.

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(dplyr)
})
source("./PPRE_threshold_sensitivity/00_config.R")

interval_list <- readRDS(file.path(out_data, "deduplicated_intervals_by_threshold.rds"))

read_promoter_bed <- function(bed) {
  df <- fread(bed, header = FALSE,
              col.names = c("chrom", "start0", "end", "gene", "score", "strand"))
  GRanges(seqnames = df$chrom, ranges = IRanges(start = df$start0 + 1L, end = df$end),
          strand = df$strand, gene = df$gene)
}
gene_sets <- list(
  nonzero = read_promoter_bed(promoter_windows_nonzero_bed),
  top20   = read_promoter_bed(promoter_windows_top20_bed)
)

views <- c("heterodimer", "PPAR_half", "RXR_half", "combined")
view_filter <- function(families_present, view) {
  if (view == "combined") return(rep(TRUE, length(families_present)))
  grepl(view, families_present, fixed = TRUE)
}

promoter_overlap_rows <- list()
gene_score_rows <- list()

for (tag in names(interval_list)) {
  itv <- interval_list[[tag]]
  pt <- itv$threshold[1]

  for (view in views) {
    itv_v <- itv[view_filter(itv$families_present, view)]
    if (nrow(itv_v) == 0) next
    gr_v <- GRanges(seqnames = itv_v$chrom, ranges = IRanges(start = itv_v$start, end = itv_v$end))

    for (gs_name in names(gene_sets)) {
      gs <- gene_sets[[gs_name]]
      ov <- findOverlaps(gr_v, gs, ignore.strand = TRUE)
      if (length(ov) == 0) {
        promoter_overlap_rows[[length(promoter_overlap_rows) + 1]] <- data.table(
          threshold = pt, view = view, gene_set = gs_name,
          n_promoter_intervals = 0L, n_genes_assigned = 0L)
        next
      }
      hit_genes <- gs$gene[subjectHits(ov)]
      dt_ov <- data.table(
        interval_idx = queryHits(ov),
        gene = hit_genes,
        neglog10p = -log10(itv_v$best_pvalue[queryHits(ov)]),
        best_score = itv_v$best_score[queryHits(ov)]
      )

      promoter_overlap_rows[[length(promoter_overlap_rows) + 1]] <- data.table(
        threshold = pt, view = view, gene_set = gs_name,
        n_promoter_intervals = uniqueN(dt_ov$interval_idx),
        n_genes_assigned = uniqueN(dt_ov$gene)
      )

      gene_scores <- dt_ov[, .(
        n_intervals   = uniqueN(interval_idx),
        max_neglog10p = max(neglog10p),
        sum_neglog10p = sum(neglog10p),
        max_score     = max(best_score)
      ), by = gene]
      gene_scores[, `:=`(threshold = pt, view = view, gene_set = gs_name)]
      gene_score_rows[[length(gene_score_rows) + 1]] <- gene_scores
    }
  }
  cat("Processed", tag, "\n")
}

promoter_overlap_summary <- rbindlist(promoter_overlap_rows)
gene_scores_all <- rbindlist(gene_score_rows)

fwrite(promoter_overlap_summary, file.path(out_results, "promoter_overlap_summary_by_threshold.csv"))
fwrite(gene_scores_all, file.path(out_results, "gene_scores_by_threshold.csv"))

## ---- rank stability vs. the p<=1e-4 baseline, per view and gene universe, kept separate ----
# NOTE: loop variables are deliberately named cur_view/cur_gs, NOT view/gene_set --
# those names collide with actual columns of gene_scores_all, and inside a
# data.table `[...]` filter, get("view") resolves within data.table's own
# column-overlay evaluation frame (returning the whole column) rather than the
# calling script's scope, silently turning `view == get("view")` into an
# always-true per-row comparison. Distinct names sidestep the issue entirely.
rank_rows <- list()
for (cur_view in views) for (cur_gs in names(gene_sets)) {
  base <- gene_scores_all[threshold == baseline_threshold & view == cur_view & gene_set == cur_gs]
  if (nrow(base) == 0) next
  setorder(base, -sum_neglog10p)
  base[, rank_baseline := .I]

  for (pt in setdiff(p_thresholds, baseline_threshold)) {
    cmp <- gene_scores_all[threshold == pt & view == cur_view & gene_set == cur_gs]
    if (nrow(cmp) == 0) next
    setorder(cmp, -sum_neglog10p)
    cmp[, rank_cmp := .I]

    m <- merge(base[, .(gene, rank_baseline, sum_neglog10p_baseline = sum_neglog10p)],
               cmp[,  .(gene, rank_cmp,      sum_neglog10p_cmp      = sum_neglog10p)],
               by = "gene", all = TRUE)
    # genes absent from one side (didn't clear that threshold's promoter overlap
    # at all) are given the worst possible rank (one past the last observed gene)
    m[is.na(rank_baseline), rank_baseline := nrow(base) + 1L]
    m[is.na(rank_cmp),      rank_cmp      := nrow(cmp)  + 1L]

    rho <- suppressWarnings(cor(m$rank_baseline, m$rank_cmp, method = "spearman"))

    top_overlap <- function(n) {
      top_base <- base$gene[seq_len(min(n, nrow(base)))]
      top_cmp  <- cmp$gene[seq_len(min(n, nrow(cmp)))]
      length(intersect(top_base, top_cmp))
    }

    rank_rows[[length(rank_rows) + 1]] <- data.table(
      view = cur_view, gene_set = cur_gs,
      threshold_baseline = baseline_threshold, threshold_compare = pt,
      n_genes_baseline = nrow(base), n_genes_compare = nrow(cmp),
      spearman_rho = rho,
      top100_overlap = top_overlap(100), top100_frac = top_overlap(100) / 100,
      top500_overlap = top_overlap(500), top500_frac = top_overlap(500) / 500
    )
  }
}
rank_stability <- rbindlist(rank_rows)
fwrite(rank_stability, file.path(out_results, "rank_stability_vs_baseline.csv"))

cat("\n=== Promoter overlap summary ===\n"); print(promoter_overlap_summary)
cat("\n=== Rank stability vs. p<=1e-4 baseline ===\n"); print(rank_stability)
message("Done. Written to ", out_results)
