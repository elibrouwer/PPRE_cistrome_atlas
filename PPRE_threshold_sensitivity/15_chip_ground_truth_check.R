# The one validation that passed so far (PPAR_DR1_Q2) is itself another
# SEQUENCE-based motif prediction (TRANSFAC), not experimental ground truth --
# two motif scanners agreeing doesn't rule out both picking up the same
# sequence-composition artifact rather than real PPAR biology. This is the
# more meaningful check: does the winning method (poisson_neglog10p_heterodimer)
# actually predict real PPARA ChIP-seq occupancy (liver, from
# PPARA_chipseq_overlap/chipatlas_liver_peaks_raw.bed, already used in
# 07_chip_validation.R), any better than nothing?

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(fgsea)
})
source("./PPRE_threshold_sensitivity/00_config.R")

scores <- fread(file.path(out_results, "gene_scoring_size_normalized.csv"))

bed <- fread(promoter_windows_nonzero_bed, header = FALSE,
             col.names = c("chrom", "start0", "end", "gene", "score", "strand"))
promo_gr <- GRanges(seqnames = bed$chrom, ranges = IRanges(start = bed$start0 + 1L, end = bed$end), gene = bed$gene)

chip_bed <- file.path(repo_dir, "PPARA_chipseq_overlap", "chipatlas_liver_peaks_raw.bed")
raw <- fread(chip_bed, header = FALSE, col.names = c("chrom", "start", "end"))
std_chroms <- paste0("chr", c(1:22, "X", "Y", "M"))
raw <- raw[chrom %in% std_chroms]
chip_gr <- reduce(GRanges(raw$chrom, IRanges(raw$start + 1L, raw$end)))

ov <- findOverlaps(promo_gr, chip_gr, ignore.strand = TRUE)
chip_supported_genes <- unique(promo_gr$gene[queryHits(ov)])
cat(sprintf("%d / %d nonzero-expressed genes (%.1f%%) have a real liver PPARA ChIP peak somewhere in their promoter.\n",
            length(chip_supported_genes), length(unique(promo_gr$gene)), 100*length(chip_supported_genes)/length(unique(promo_gr$gene))))

## ---- does the winning method's ranking predict ChIP overlap, genome-wide (fgsea) ----
set.seed(1)
methods_to_check <- c("poisson_neglog10p_heterodimer", "poisson_neglog10p_combined",
                       "poisson_neglog10p_PPAR_half", "poisson_neglog10p_RXR_half")
res_list <- list()
for (m in methods_to_check) {
  stat <- setNames(scores[[m]], scores$gene)
  stat <- stat + rnorm(length(stat), 0, 1e-9 * (max(stat)-min(stat)+1e-9))
  r <- fgsea(pathways = list(chip_supported = chip_supported_genes), stats = stat,
             minSize = 5, maxSize = length(chip_supported_genes) + 100, eps = 0, scoreType = "pos")
  res_list[[m]] <- as.data.table(r)[, method := m]
}
res <- rbindlist(res_list, fill = TRUE)[, .(method, NES, pval, padj, size)]
fwrite(res, file.path(out_results, "chip_ground_truth_validation.csv"))
cat("\n=== Does each method's ranking predict REAL liver PPARA ChIP overlap? (fgsea) ===\n")
print(res)

## ---- simpler, more interpretable check: ChIP-overlap rate by score decile ----
scores[, chip_supported := gene %in% chip_supported_genes]
scores[, decile := cut(rank(-poisson_neglog10p_heterodimer), breaks = quantile(rank(-poisson_neglog10p_heterodimer), probs = 0:10/10),
                        include.lowest = TRUE, labels = FALSE)]
decile_tbl <- scores[, .(n = .N, n_chip_supported = sum(chip_supported), frac_chip_supported = mean(chip_supported)), by = decile][order(decile)]
cat("\n=== Real liver-PPARA-ChIP-overlap rate by poisson_neglog10p_heterodimer decile (1=top decile) ===\n")
print(decile_tbl)

## ---- direct comparison: is this better than the OLD threshold-based approach
## (just using p<=1e-4 heterodimer hit presence, no size correction, no Poisson)? ----
old_approach <- fread(file.path(out_results, "gene_scores_by_threshold.csv"))[
  threshold == 1e-4 & view == "heterodimer" & gene_set == "nonzero", gene]
scores[, old_approach_flag := gene %in% old_approach]
cat(sprintf("\nOld approach (binary: any heterodimer hit at p<=1e-4): %d/%d genes flagged (%.1f%%), of which %.1f%% have real ChIP support\n",
            length(old_approach), nrow(scores), 100*length(old_approach)/nrow(scores),
            100*mean(scores$chip_supported[scores$gene %in% old_approach])))
cat(sprintf("Overall ChIP-support rate (background): %.1f%%\n", 100*mean(scores$chip_supported)))
cat(sprintf("New method, top decile ChIP-support rate: %.1f%%\n", decile_tbl$frac_chip_supported[1]))
cat(sprintf("New method, bottom decile ChIP-support rate: %.1f%%\n", decile_tbl$frac_chip_supported[10]))

message("\nDone.")
