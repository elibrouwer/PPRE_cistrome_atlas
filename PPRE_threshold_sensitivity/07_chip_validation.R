# Validate FIMO threshold/family choices against real ChIP-seq peaks,
# genome-wide (not restricted to the 5 genes PPARA_chipseq_overlap covered).
#
# Ground truth: ChIP-Atlas PPARA peaks aggregated across 7 liver-context
# experiments (PPARA_chipseq_overlap/chipatlas_liver_peaks_raw.bed,
# 1,371,392 raw records -> reduced to non-redundant loci here, same
# reduce() step 06_make_liver_centrimo_fasta.R already used). This is PPARA
# antigen, liver tissue -- NOT cardiomyocyte and not PPARD/PPARG -- so per
# the project's own validation convention (Objective1_fimo_chip_validation.R),
# this is explicitly INDIRECT validation, reported as such, not proof of
# cardiac occupancy. Absence from a ChIP peak is not evidence of no binding
# (different tissue/antigen/condition).
#
# Method: same fold-enrichment-over-genome-background + one-sided binomial
# test already used in this project's Objective1_fimo_threshold_enrichment.R
# (simpler than, and deliberately not extending, PPARD_threshold_validation's
# permutation-based framework -- see README for why). Applied per threshold
# AND per motif family (heterodimer / PPAR_half / RXR_half / combined) so a
# family-specific signal isn't averaged away; RXR_half is included on the
# rationale that RXR half-sites are part of the same DR1 element PPARA
# ChIP would capture via heterodimeric co-occupancy, not because RXR itself
# was the ChIP antigen -- flagged explicitly in the output.
#
# No RXR ChIP-seq or cardiac-tissue PPAR ChIP-seq was found anywhere in this
# project during the audit (see README) -- this is the only ChIP validation
# currently possible at genome-wide scale with the data on hand.

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
})
source("./PPRE_threshold_sensitivity/00_config.R")

chip_bed <- file.path(repo_dir, "PPARA_chipseq_overlap", "chipatlas_liver_peaks_raw.bed")
stopifnot(file.exists(chip_bed))

raw <- fread(chip_bed, header = FALSE, col.names = c("chrom", "start", "end"))
std_chroms <- paste0("chr", c(1:22, "X", "Y", "M"))
raw <- raw[chrom %in% std_chroms]
chip_gr <- reduce(GRanges(raw$chrom, IRanges(raw$start + 1L, raw$end)))
cat(sprintf("ChIP-Atlas liver PPARA peaks: %s raw records -> %s non-redundant loci, %s bp total\n",
            format(nrow(raw), big.mark = ","), format(length(chip_gr), big.mark = ","),
            format(sum(width(chip_gr)), big.mark = ",")))

genome_bp <- 3.05e9   # same constant used in Objective1_fimo_threshold_enrichment.R
expected_fraction <- sum(width(chip_gr)) / genome_bp

interval_list <- readRDS(file.path(out_data, "deduplicated_intervals_by_threshold.rds"))
views <- c("heterodimer", "PPAR_half", "RXR_half", "combined")
view_filter <- function(families_present, view) {
  if (view == "combined") return(rep(TRUE, length(families_present)))
  grepl(view, families_present, fixed = TRUE)
}

rows <- list()
for (tag in names(interval_list)) {
  itv <- interval_list[[tag]]
  pt <- itv$threshold[1]
  for (view in views) {
    itv_v <- itv[view_filter(itv$families_present, view)]
    n_total <- nrow(itv_v)
    if (n_total == 0) next
    gr_v <- GRanges(seqnames = itv_v$chrom, ranges = IRanges(start = itv_v$start, end = itv_v$end))
    n_in <- length(unique(queryHits(findOverlaps(gr_v, chip_gr, ignore.strand = TRUE, minoverlap = 1))))
    obs_frac <- n_in / n_total
    fold <- obs_frac / expected_fraction
    p <- binom.test(n_in, n_total, p = expected_fraction, alternative = "greater")$p.value
    rows[[length(rows) + 1]] <- data.table(
      threshold = pt, view = view,
      n_total_intervals = n_total, n_overlapping_chip = n_in,
      observed_fraction = obs_frac, fold_enrichment = fold, binom_p = p
    )
  }
}
chip_validation <- rbindlist(rows)
fwrite(chip_validation, file.path(out_results, "chip_validation_PPARA_liver_by_threshold.csv"))
cat("\n=== Genome-wide fold enrichment vs. ChIP-Atlas liver PPARA peaks (INDIRECT validation: liver, not cardiac) ===\n")
print(chip_validation)
message("Done. Written to ", file.path(out_results, "chip_validation_PPARA_liver_by_threshold.csv"))
