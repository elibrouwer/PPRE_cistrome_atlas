# Does scanning the repeat-masked genome (FIMO_results_genome_masked, already
# in the project, generated to address the ALU-repeat-driven inflation this
# project's own notes flagged) produce a gene-level score that predicts real
# occupancy any better than the unmasked run? Same rigor, same two ground
# truths (liver PPARA ChIP, cardiac HCM+PLN accessibility), same
# size-normalized Poisson scoring approach that performed best among
# everything tried on the unmasked run (15/18: 0.773/0.770 vs ChIP,
# 0.588/0.609 vs cardiac).

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

masked_dir <- file.path(base_dir, "FIMO_results_genome_masked")
motif_family <- fread(motif_family_csv)

## ---- 1. load masked-run hits (p<=1e-4), same 30 motifs, same shards ----
hits_list <- vector("list", length(fimo_shards)); names(hits_list) <- fimo_shards
for (shard in fimo_shards) {
  f <- file.path(masked_dir, shard, "fimo.tsv")
  stopifnot(file.exists(f))
  dt <- fread(f, sep = "\t", header = TRUE, fill = TRUE)
  dt <- dt[!is.na(start) & !startsWith(motif_id, "#")]
  dt <- dt[, .(motif_id, sequence_name = as.character(sequence_name), start, stop, strand, score, pvalue = `p-value`)]
  hits_list[[shard]] <- dt
  cat(sprintf("%s (masked): %s rows\n", shard, format(nrow(dt), big.mark=",")))
}
masked_hits <- rbindlist(hits_list); rm(hits_list); gc()
masked_hits <- merge(masked_hits, motif_family[, .(motif_id, family)], by = "motif_id", all.x = TRUE)
cat(sprintf("\nTotal masked-run hits (p<=1e-4): %s\n", format(nrow(masked_hits), big.mark=",")))

## ---- 2. dedup into intervals (combined + heterodimer view only, to keep this fast) ----
seqn <- ifelse(masked_hits$sequence_name %in% c("MT","M"), "chrM", paste0("chr", masked_hits$sequence_name))
std_chroms <- paste0("chr", c(1:22,"X","Y","M"))
keep <- seqn %in% std_chroms
masked_hits <- masked_hits[keep]; seqn <- seqn[keep]
gr_all <- GRanges(seqn, IRanges(masked_hits$start, masked_hits$stop))

merged_combined <- reduce(gr_all, ignore.strand = TRUE)
ov_c <- findOverlaps(gr_all, merged_combined, ignore.strand = TRUE)
interval_agg_combined <- data.table(interval_id = subjectHits(ov_c), pvalue = masked_hits$pvalue[queryHits(ov_c)])[
  , .(best_pvalue = min(pvalue)), by = interval_id]
itv_combined <- data.table(chrom = as.character(seqnames(merged_combined)), start = start(merged_combined),
                            end = end(merged_combined), interval_id = seq_along(merged_combined))
itv_combined <- merge(itv_combined, interval_agg_combined, by = "interval_id")

het_hits <- masked_hits[family == "heterodimer"]
gr_het <- GRanges(ifelse(het_hits$sequence_name %in% c("MT","M"), "chrM", paste0("chr", het_hits$sequence_name)),
                   IRanges(het_hits$start, het_hits$stop))
merged_het <- reduce(gr_het, ignore.strand = TRUE)
ov_h <- findOverlaps(gr_het, merged_het, ignore.strand = TRUE)
interval_agg_het <- data.table(interval_id = subjectHits(ov_h), pvalue = het_hits$pvalue[queryHits(ov_h)])[
  , .(best_pvalue = min(pvalue)), by = interval_id]
itv_het <- data.table(chrom = as.character(seqnames(merged_het)), start = start(merged_het),
                       end = end(merged_het), interval_id = seq_along(merged_het))
itv_het <- merge(itv_het, interval_agg_het, by = "interval_id")

cat(sprintf("Masked-run deduplicated intervals: combined=%s, heterodimer=%s\n",
            format(nrow(itv_combined), big.mark=","), format(nrow(itv_het), big.mark=",")))

## ---- 3. per-gene size-normalized Poisson score (same method as 13, using the
## same isoform-unioned promoter windows for exact comparability) ----
bed <- fread(promoter_windows_nonzero_bed, header = FALSE, col.names = c("chrom","start0","end","gene","score","strand"))
gene_width <- bed[, .(promoter_bp = sum(end - start0)), by = gene]
promo_gr <- GRanges(bed$chrom, IRanges(bed$start0+1L, bed$end), gene = bed$gene)

score_gene <- function(itv_tbl, tag) {
  gr_v <- GRanges(itv_tbl$chrom, IRanges(itv_tbl$start, itv_tbl$end))
  ov <- findOverlaps(gr_v, promo_gr, ignore.strand = TRUE)
  dt_ov <- data.table(gene = promo_gr$gene[subjectHits(ov)], pvalue = itv_tbl$best_pvalue[queryHits(ov)])
  n_loci <- dt_ov[, .N, by = gene]
  genome_bp <- 3.05e9
  rate <- nrow(itv_tbl) / genome_bp
  g <- merge(gene_width, n_loci, by = "gene", all.x = TRUE)
  g[is.na(N), N := 0L]
  g[, expected := promoter_bp * rate]
  g[, (paste0("masked_poisson_neglog10p_", tag)) := -ppois(N - 1L, expected, lower.tail = FALSE, log.p = TRUE) / log(10)]
  g[, c("gene", paste0("masked_poisson_neglog10p_", tag)), with = FALSE]
}
scored <- merge(score_gene(itv_combined, "combined"), score_gene(itv_het, "heterodimer"), by = "gene")
fwrite(scored, file.path(out_results, "repeatmasked_gene_scores.csv"))

## ---- 4. validate against both ground truths, same framework as 15/18 ----
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
dat <- merge(dat, gene_width, by = "gene", all.x = TRUE)

feature_cols <- c("masked_poisson_neglog10p_combined", "masked_poisson_neglog10p_heterodimer")
set.seed(42); K <- 5
dat[, fold := sample(rep(1:K, length.out = .N))]

run_cv <- function(label_col) {
  dat[, pred_rf := NA_real_]; dat[, pred_base := NA_real_]
  for (k in 1:K) {
    train <- dat[fold != k]; test_idx <- dat[, which(fold == k)]
    rf <- ranger(x = train[, ..feature_cols], y = factor(train[[label_col]]), probability = TRUE, num.trees = 500, seed = 42)
    dat[test_idx, pred_rf := predict(rf, data = dat[test_idx, ..feature_cols])$predictions[, "1"]]
    base_model <- glm(as.formula(paste(label_col, "~ promoter_bp")), data = train, family = binomial())
    dat[test_idx, pred_base := predict(base_model, newdata = dat[test_idx], type = "response")]
  }
  list(auc_rf = as.numeric(auc(dat[[label_col]], dat$pred_rf, quiet=TRUE)),
       auc_base = as.numeric(auc(dat[[label_col]], dat$pred_base, quiet=TRUE)))
}
res_chip <- run_cv("chip_supported")
res_cardiac <- run_cv("cardiac_supported")

cat("\n=== Repeat-masked-genome run, cross-validated AUC (5-fold) ===\n")
cat(sprintf("vs. real LIVER PPARA ChIP:      masked-run score = %.3f | promoter-size-only baseline = %.3f\n", res_chip$auc_rf, res_chip$auc_base))
cat(sprintf("vs. real CARDIAC accessibility: masked-run score = %.3f | promoter-size-only baseline = %.3f\n", res_cardiac$auc_rf, res_cardiac$auc_base))
cat("(Unmasked-genome run, for reference: 0.773/0.770 vs ChIP, 0.588/0.609 vs cardiac)\n")

fwrite(dat, file.path(out_results, "repeatmasked_validated.csv"))
message("\nDone.")
