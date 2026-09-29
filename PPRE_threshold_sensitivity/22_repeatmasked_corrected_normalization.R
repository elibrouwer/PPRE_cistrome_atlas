# Redo the repeat-masked-run validation (21) with a corrected Poisson
# denominator: actual non-repeat (scannable) bp per promoter window, computed
# by subtracting each window's overlap with RepeatMasker regions
# (rmsk_hg38.bed, the same annotation FIMO_results_genome_masked's masking
# was built from), instead of the nominal window bp used in 21 -- which
# implicitly assumed the whole window was scannable and wasn't, for any gene
# with repeat content in its promoter. Also computes the genome-wide rate
# using total NON-REPEAT genome bp as the denominator, for the same reason.

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(rtracklayer)
  library(readxl)
  library(ranger)
  library(pROC)
})
source("./PPRE_threshold_sensitivity/00_config.R")

## ---- 1. repeat regions, and actual scannable bp per promoter window ----
rmsk <- fread(file.path(out_data, "rmsk_hg38.bed"), header = FALSE, col.names = c("chrom","start","end"))
std_chroms <- paste0("chr", c(1:22,"X","Y","M"))
rmsk <- rmsk[chrom %in% std_chroms]
rmsk_gr <- reduce(GRanges(rmsk$chrom, IRanges(rmsk$start+1L, rmsk$end)))
cat(sprintf("RepeatMasker regions: %s non-redundant loci, %s bp total\n",
            format(length(rmsk_gr), big.mark=","), format(sum(width(rmsk_gr)), big.mark=",")))

bed <- fread(promoter_windows_nonzero_bed, header = FALSE, col.names = c("chrom","start0","end","gene","score","strand"))
promo_gr <- GRanges(bed$chrom, IRanges(bed$start0+1L, bed$end), gene = bed$gene)

# scannable bp per window = window width - overlap with repeats.
# Vectorized via coverage()+viewSums() (one RLE per chromosome, then a single
# vectorized sum per window) instead of calling intersect() once per window --
# the per-window loop this replaced was effectively O(n) GRanges calls and
# never finished in reasonable time for ~40k windows.
rmsk_gr_trim <- rmsk_gr[seqnames(rmsk_gr) %in% seqnames(promo_gr)]
seqlevels(rmsk_gr_trim) <- seqlevels(promo_gr)
cov <- coverage(rmsk_gr_trim)
common_chroms <- intersect(names(cov), seqlevels(promo_gr))
promo_gr <- promo_gr[as.character(seqnames(promo_gr)) %in% common_chroms]
bed <- bed[bed$chrom %in% common_chroms]
repeat_bp_per_window <- integer(length(promo_gr))
for (ch in common_chroms) {
  idx <- which(as.character(seqnames(promo_gr)) == ch)
  if (length(idx) == 0) next
  v <- Views(cov[[ch]], start = start(promo_gr)[idx], end = pmin(end(promo_gr)[idx], length(cov[[ch]])))
  repeat_bp_per_window[idx] <- viewSums(v)
}

window_dt <- data.table(gene = promo_gr$gene, window_bp = width(promo_gr), repeat_bp = repeat_bp_per_window)
window_dt[, scannable_bp := window_bp - repeat_bp]
gene_scannable <- window_dt[, .(promoter_bp_nominal = sum(window_bp), promoter_bp_scannable = sum(scannable_bp)), by = gene]
gene_scannable[, pct_repeat_masked := 100 * (1 - promoter_bp_scannable / promoter_bp_nominal)]

cat(sprintf("\nMedian %% of promoter bp that is repeat-masked (unscannable): %.1f%%\n",
            median(gene_scannable$pct_repeat_masked)))
cat("Genes with the most repeat-masked promoter (this is exactly where the naive-bp normalization was most wrong):\n")
print(head(gene_scannable[order(-pct_repeat_masked)], 5))

## ---- 2. genome-wide non-repeat bp (for the rate denominator) ----
genome_bp_total <- 3.05e9
genome_bp_nonrepeat <- genome_bp_total - sum(width(rmsk_gr))
cat(sprintf("\nGenome-wide: %s bp total, %s bp non-repeat (%.1f%%)\n",
            format(genome_bp_total, big.mark=","), format(genome_bp_nonrepeat, big.mark=","),
            100*genome_bp_nonrepeat/genome_bp_total))

## ---- 3. reload masked-run hit counts per gene (reuse 21's already-computed n_loci logic) ----
masked_dir <- file.path(base_dir, "FIMO_results_genome_masked")
motif_family <- fread(motif_family_csv)
hits_list <- vector("list", length(fimo_shards)); names(hits_list) <- fimo_shards
for (shard in fimo_shards) {
  f <- file.path(masked_dir, shard, "fimo.tsv")
  dt <- fread(f, sep = "\t", header = TRUE, fill = TRUE)
  dt <- dt[!is.na(start) & !startsWith(motif_id, "#")]
  hits_list[[shard]] <- dt[, .(motif_id, sequence_name = as.character(sequence_name), start, stop, pvalue = `p-value`)]
}
masked_hits <- rbindlist(hits_list); rm(hits_list); gc()
masked_hits <- merge(masked_hits, motif_family[, .(motif_id, family)], by = "motif_id", all.x = TRUE)
seqn <- ifelse(masked_hits$sequence_name %in% c("MT","M"), "chrM", paste0("chr", masked_hits$sequence_name))
keep <- seqn %in% std_chroms
masked_hits <- masked_hits[keep]; seqn <- seqn[keep]

dedup_and_count <- function(hits, tag) {
  gr <- GRanges(ifelse(hits$sequence_name %in% c("MT","M"), "chrM", paste0("chr", hits$sequence_name)),
                IRanges(hits$start, hits$stop))
  merged <- reduce(gr, ignore.strand = TRUE)
  ov <- findOverlaps(gr, merged, ignore.strand = TRUE)
  n_intervals_total <- length(merged)
  ov_promo <- findOverlaps(merged, promo_gr, ignore.strand = TRUE)
  n_loci <- data.table(gene = promo_gr$gene[subjectHits(ov_promo)])[, .N, by = gene]
  list(n_total = n_intervals_total, n_loci = n_loci)
}
res_combined <- dedup_and_count(masked_hits, "combined")
res_het <- dedup_and_count(masked_hits[family == "heterodimer"], "heterodimer")
cat(sprintf("\nMasked-run deduplicated intervals: combined=%s, heterodimer=%s\n",
            format(res_combined$n_total, big.mark=","), format(res_het$n_total, big.mark=",")))

## ---- 4. CORRECTED Poisson score: scannable bp denominator, non-repeat genome-wide rate ----
score_corrected <- function(n_loci, n_total, tag) {
  rate <- n_total / genome_bp_nonrepeat   # corrected: rate per scannable bp, not per nominal bp
  g <- merge(gene_scannable, n_loci, by = "gene", all.x = TRUE)
  g[is.na(N), N := 0L]
  g[, expected := promoter_bp_scannable * rate]   # corrected: scannable bp, not nominal bp
  g[, (paste0("masked_corrected_neglog10p_", tag)) := -ppois(N - 1L, expected, lower.tail = FALSE, log.p = TRUE) / log(10)]
  g[, c("gene", paste0("masked_corrected_neglog10p_", tag)), with = FALSE]
}
scored <- merge(score_corrected(res_combined$n_loci, res_combined$n_total, "combined"),
                 score_corrected(res_het$n_loci, res_het$n_total, "heterodimer"), by = "gene")
scored <- merge(scored, gene_scannable, by = "gene")
fwrite(scored, file.path(out_results, "repeatmasked_gene_scores_corrected.csv"))

## ---- 5. validate exactly as before (same ground truths, same CV framework) ----
raw_chip <- fread(file.path(repo_dir, "PPARA_chipseq_overlap", "chipatlas_liver_peaks_raw.bed"),
                   header = FALSE, col.names = c("chrom","start","end"))[chrom %in% std_chroms]
chip_gr <- reduce(GRanges(raw_chip$chrom, IRanges(raw_chip$start+1L, raw_chip$end)))
chip_genes <- unique(promo_gr$gene[queryHits(findOverlaps(promo_gr, chip_gr, ignore.strand = TRUE))])

hcm <- read_xlsx(file.path(base_dir, "HCM_differentially_regions.xlsx"), sheet = "Supplementary Table 2A")
gr_hcm <- GRanges(hcm$`Chromosome number`, IRanges(hcm$`Region start`, hcm$`Region end`))
pln <- read_xlsx(file.path(base_dir, "PLN_differentially_regions.xlsx"), skip = 1)
pln <- pln[!is.na(pln$`Region start†`), ]
gr_pln <- GRanges(pln$Chromosome, IRanges(pln$`Region start†`, pln$`Region end†`))
chain <- import.chain(file.path(base_dir, "hg19ToHg38.over.chain"))
cardiac_gr <- reduce(c(unlist(liftOver(gr_hcm, chain)), unlist(liftOver(gr_pln, chain))))
seqlevelsStyle(cardiac_gr) <- "UCSC"
cardiac_genes <- unique(promo_gr$gene[queryHits(findOverlaps(promo_gr, cardiac_gr, ignore.strand = TRUE))])

dat <- copy(scored)
dat[, chip_supported := as.integer(gene %in% chip_genes)]
dat[, cardiac_supported := as.integer(gene %in% cardiac_genes)]

feature_cols <- c("masked_corrected_neglog10p_combined", "masked_corrected_neglog10p_heterodimer")
set.seed(42); K <- 5
dat[, fold := sample(rep(1:K, length.out = .N))]

run_cv <- function(label_col) {
  dat[, pred_rf := NA_real_]; dat[, pred_base := NA_real_]
  for (k in 1:K) {
    train <- dat[fold != k]; test_idx <- dat[, which(fold == k)]
    rf <- ranger(x = train[, ..feature_cols], y = factor(train[[label_col]]), probability = TRUE, num.trees = 500, seed = 42)
    dat[test_idx, pred_rf := predict(rf, data = dat[test_idx, ..feature_cols])$predictions[, "1"]]
    base_model <- glm(as.formula(paste(label_col, "~ promoter_bp_scannable")), data = train, family = binomial())
    dat[test_idx, pred_base := predict(base_model, newdata = dat[test_idx], type = "response")]
  }
  list(auc_rf = as.numeric(auc(dat[[label_col]], dat$pred_rf, quiet=TRUE)),
       auc_base = as.numeric(auc(dat[[label_col]], dat$pred_base, quiet=TRUE)))
}
res_chip <- run_cv("chip_supported")
res_cardiac <- run_cv("cardiac_supported")

cat("\n=== Repeat-masked-genome run, CORRECTED (scannable-bp) normalization, CV AUC ===\n")
cat(sprintf("vs. real LIVER PPARA ChIP:      corrected score = %.3f | scannable-size-only baseline = %.3f\n", res_chip$auc_rf, res_chip$auc_base))
cat(sprintf("vs. real CARDIAC accessibility: corrected score = %.3f | scannable-size-only baseline = %.3f\n", res_cardiac$auc_rf, res_cardiac$auc_base))
cat("\n(For reference:\n")
cat("  Unmasked run:                    0.773/0.770 vs ChIP, 0.588/0.609 vs cardiac\n")
cat("  Masked run, NAIVE normalization: 0.625/0.770 vs ChIP, 0.546/0.609 vs cardiac)\n")

fwrite(dat, file.path(out_results, "repeatmasked_corrected_validated.csv"))
message("\nDone.")
