# The user's question: since no method-based scoring beats a trivial baseline
# against LIVER PPARA ChIP-seq (07, 15, 16), would a genuinely tissue-matched
# functional dataset do better? This project already has one, mostly unused:
# HCM_differentially_regions.xlsx / PLN_differentially_regions.xlsx --
# differentially-accessible chromatin regions from real human cardiac tissue
# (hypertrophic cardiomyopathy / phospholamban-mutation dilated cardiomyopathy
# patients vs. controls; hg19; used only once before, for a single motif, in
# Objective1_fimo_threshold_enrichment.R). Not PPAR-specific and not
# baseline/healthy tissue, but it IS real, human, cardiac chromatin state --
# a materially better tissue match than liver ChIP-seq.
#
# Tested the same two rigorous ways as the liver ChIP check (07/15/16):
#   (a) genome-wide fold enrichment of deduplicated intervals, by threshold
#       and motif family (same binomial method as 07)
#   (b) gene-level: does a gene's promoter overlapping a cardiac
#       differentially-accessible region correlate with its motif score, via
#       cross-validated random forest AUC vs. a promoter-size-only baseline
#       (same framework as 16, so directly comparable to the liver-ChIP AUCs
#       of 0.773 vs 0.770)

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(rtracklayer)
  library(ranger)
  library(pROC)
})
source("./PPRE_threshold_sensitivity/00_config.R")

## ---- load + liftover HCM + PLN regions (hg19 -> hg38), same as Objective1_fimo_threshold_enrichment.R ----
hcm <- read_xlsx(file.path(base_dir, "HCM_differentially_regions.xlsx"), sheet = "Supplementary Table 2A")
gr_hcm <- GRanges(seqnames = hcm$`Chromosome number`, ranges = IRanges(start = hcm$`Region start`, end = hcm$`Region end`))

pln <- read_xlsx(file.path(base_dir, "PLN_differentially_regions.xlsx"), skip = 1)
pln <- pln[!is.na(pln$`Region start†`), ]
gr_pln <- GRanges(seqnames = pln$Chromosome, ranges = IRanges(start = pln$`Region start†`, end = pln$`Region end†`))

chain <- import.chain(file.path(base_dir, "hg19ToHg38.over.chain"))
gr_hcm_38 <- unlist(liftOver(gr_hcm, chain))
gr_pln_38 <- unlist(liftOver(gr_pln, chain))
cardiac_gr <- reduce(c(gr_hcm_38, gr_pln_38))
seqlevelsStyle(cardiac_gr) <- "UCSC"
cat(sprintf("Cardiac differentially-accessible regions (HCM %d + PLN %d, hg19->hg38, merged): %d non-redundant loci, %s bp total (%.2f%% of genome)\n",
            length(gr_hcm), length(gr_pln), length(cardiac_gr), format(sum(width(cardiac_gr)), big.mark=","),
            100*sum(width(cardiac_gr))/3.05e9))

## ---- (a) genome-wide fold enrichment, by threshold and family (same method as 07) ----
interval_list <- readRDS(file.path(out_data, "deduplicated_intervals_by_threshold.rds"))
genome_bp <- 3.05e9
expected_fraction <- sum(width(cardiac_gr)) / genome_bp
views <- c("heterodimer", "PPAR_half", "RXR_half", "combined")
view_filter <- function(fp, v) if (v == "combined") rep(TRUE, length(fp)) else grepl(v, fp, fixed = TRUE)

rows <- list()
for (tag in names(interval_list)) {
  itv <- interval_list[[tag]]; pt <- itv$threshold[1]
  for (v in views) {
    itv_v <- itv[view_filter(itv$families_present, v)]
    n_total <- nrow(itv_v); if (n_total == 0) next
    gr_v <- GRanges(itv_v$chrom, IRanges(itv_v$start, itv_v$end))
    n_in <- length(unique(queryHits(findOverlaps(gr_v, cardiac_gr, ignore.strand = TRUE, minoverlap = 1))))
    obs_frac <- n_in / n_total; fold <- obs_frac / expected_fraction
    p <- binom.test(n_in, n_total, p = expected_fraction, alternative = "greater")$p.value
    rows[[length(rows)+1]] <- data.table(threshold = pt, view = v, n_total_intervals = n_total,
                                          n_overlapping_cardiac = n_in, fold_enrichment = fold, binom_p = p)
  }
}
cardiac_enrich <- rbindlist(rows)
fwrite(cardiac_enrich, file.path(out_results, "cardiac_accessibility_enrichment_by_threshold.csv"))
cat("\n=== Genome-wide fold enrichment vs. REAL cardiac differentially-accessible chromatin (HCM+PLN) ===\n")
print(cardiac_enrich)

## ---- (b) gene-level: does motif score predict cardiac-region promoter overlap? ----
bed <- fread(promoter_windows_nonzero_bed, header = FALSE, col.names = c("chrom","start0","end","gene","score","strand"))
promo_gr <- GRanges(bed$chrom, IRanges(bed$start0+1L, bed$end), gene = bed$gene)
cardiac_supported_genes <- unique(promo_gr$gene[queryHits(findOverlaps(promo_gr, cardiac_gr, ignore.strand = TRUE))])
cat(sprintf("\n%d / %d nonzero-expressed genes (%.1f%%) have a real cardiac differentially-accessible region in their promoter.\n",
            length(cardiac_supported_genes), length(unique(promo_gr$gene)), 100*length(cardiac_supported_genes)/length(unique(promo_gr$gene))))

scores <- fread(file.path(out_results, "gene_scoring_size_normalized.csv"))
dat <- copy(scores)
dat[, cardiac_supported := as.integer(gene %in% cardiac_supported_genes)]

feature_cols <- c("fisher_n_loci_combined","fisher_n_loci_heterodimer","fisher_n_loci_PPAR_half","fisher_n_loci_RXR_half",
                   "promoter_bp","n_windows",
                   "poisson_fold_combined","poisson_fold_heterodimer","poisson_fold_PPAR_half","poisson_fold_RXR_half")

set.seed(42); K <- 5
dat[, fold := sample(rep(1:K, length.out = .N))]
dat[, `:=`(pred_rf = NA_real_, pred_size_only = NA_real_)]
for (k in 1:K) {
  train <- dat[fold != k]; test_idx <- dat[, which(fold == k)]
  rf <- ranger(x = train[, ..feature_cols], y = factor(train$cardiac_supported), probability = TRUE, num.trees = 500, seed = 42)
  dat[test_idx, pred_rf := predict(rf, data = dat[test_idx, ..feature_cols])$predictions[, "1"]]
  size_model <- glm(cardiac_supported ~ promoter_bp, data = train, family = binomial())
  dat[test_idx, pred_size_only := predict(size_model, newdata = dat[test_idx], type = "response")]
}
auc_rf <- auc(dat$cardiac_supported, dat$pred_rf, quiet = TRUE)
auc_size <- auc(dat$cardiac_supported, dat$pred_size_only, quiet = TRUE)
cat(sprintf("\nRandom forest (all features), 5-fold CV AUC vs. REAL cardiac accessibility label: %.3f\n", as.numeric(auc_rf)))
cat(sprintf("Promoter-size-ONLY baseline,  5-fold CV AUC vs. REAL cardiac accessibility label: %.3f\n", as.numeric(auc_size)))
cat("(For comparison, the same test against LIVER PPARA ChIP was 0.773 vs 0.770.)\n")

## simple univariate check: heterodimer-hit presence alone (p<=1e-4), the simplest possible "logical" filter
het_genes <- fread(file.path(out_results, "gene_scores_by_threshold.csv"))[threshold==1e-4 & view=="heterodimer" & gene_set=="nonzero", gene]
dat[, has_heterodimer := gene %in% het_genes]
cat(sprintf("\nCardiac-region overlap rate, genes WITH a heterodimer motif hit (p<=1e-4): %.1f%% (n=%d)\n",
            100*mean(dat$cardiac_supported[dat$has_heterodimer]), sum(dat$has_heterodimer)))
cat(sprintf("Cardiac-region overlap rate, genes WITHOUT a heterodimer motif hit: %.1f%% (n=%d)\n",
            100*mean(dat$cardiac_supported[!dat$has_heterodimer]), sum(!dat$has_heterodimer)))
cat(sprintf("Background rate (all genes): %.1f%%\n", 100*mean(dat$cardiac_supported)))

fwrite(dat[, .(gene, pred_rf, pred_size_only, cardiac_supported, has_heterodimer)],
       file.path(out_results, "cardiac_validation_gene_scores.csv"))
message("\nDone.")
