# Does a higher-order (order-1, dinucleotide) Markov background model change
# anything, relative to the order-0 background used throughout this analysis?
# Same unmasked genome, same 30 motifs, same p<=1e-4/--max-stored-scores
# 20000000 settings -- the ONLY change is --bfile genome_bg_order1.txt instead
# of genome_bg_order0.txt. Validated identically to the order-0 run (13/15/18):
# size-normalized Poisson score (nominal promoter_bp denominator -- this
# genome is unmasked, so no scannable-bp correction is needed here, unlike
# the repeat-masked comparison in 21/22), cross-validated RF AUC vs. a
# promoter-size-only baseline, against real liver ChIP and real cardiac
# accessibility.

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

order1_dir <- file.path(base_dir, "FIMO_results_bg_order1")
motif_family <- fread(motif_family_csv)
std_chroms <- paste0("chr", c(1:22,"X","Y","M"))

## ---- 1. load order-1-background hits (p<=1e-4), all 30 motifs, 8 shards ----
hits_list <- vector("list", length(fimo_shards)); names(hits_list) <- fimo_shards
for (shard in fimo_shards) {
  f <- file.path(order1_dir, shard, "fimo.tsv")
  stopifnot(file.exists(f))
  dt <- fread(f, sep = "\t", header = TRUE, fill = TRUE)
  dt <- dt[!is.na(start) & !startsWith(motif_id, "#")]
  hits_list[[shard]] <- dt[, .(motif_id, sequence_name = as.character(sequence_name), start, stop, pvalue = `p-value`)]
  cat(sprintf("%s (order-1 bg): %s rows\n", shard, format(nrow(dt), big.mark=",")))
}
hits <- rbindlist(hits_list); rm(hits_list); gc()
hits <- merge(hits, motif_family[, .(motif_id, family)], by = "motif_id", all.x = TRUE)
seqn <- ifelse(hits$sequence_name %in% c("MT","M"), "chrM", paste0("chr", hits$sequence_name))
keep <- seqn %in% std_chroms
hits <- hits[keep]
cat(sprintf("\nTotal order-1-background hits (p<=1e-4): %s\n", format(nrow(hits), big.mark=",")))

## ---- 2. dedup into intervals: combined + heterodimer views ----
dedup <- function(h) {
  gr <- GRanges(ifelse(h$sequence_name %in% c("MT","M"), "chrM", paste0("chr", h$sequence_name)),
                IRanges(h$start, h$stop))
  merged <- reduce(gr, ignore.strand = TRUE)
  ov <- findOverlaps(gr, merged, ignore.strand = TRUE)
  best_p <- data.table(interval_id = subjectHits(ov), pvalue = h$pvalue[queryHits(ov)])[, .(best_pvalue = min(pvalue)), by = interval_id]
  itv <- data.table(chrom = as.character(seqnames(merged)), start = start(merged), end = end(merged), interval_id = seq_along(merged))
  merge(itv, best_p, by = "interval_id")
}
itv_combined <- dedup(hits)
itv_het <- dedup(hits[family == "heterodimer"])
cat(sprintf("Order-1-background deduplicated intervals: combined=%s, heterodimer=%s\n",
            format(nrow(itv_combined), big.mark=","), format(nrow(itv_het), big.mark=",")))

## ---- 3. per-gene size-normalized Poisson score (nominal promoter_bp, same as 13/15/18) ----
bed <- fread(promoter_windows_nonzero_bed, header = FALSE, col.names = c("chrom","start0","end","gene","score","strand"))
gene_width <- bed[, .(promoter_bp = sum(end - start0)), by = gene]
promo_gr <- GRanges(bed$chrom, IRanges(bed$start0+1L, bed$end), gene = bed$gene)
genome_bp <- 3.05e9

score_gene <- function(itv_tbl, tag) {
  gr_v <- GRanges(itv_tbl$chrom, IRanges(itv_tbl$start, itv_tbl$end))
  ov <- findOverlaps(gr_v, promo_gr, ignore.strand = TRUE)
  n_loci <- data.table(gene = promo_gr$gene[subjectHits(ov)])[, .N, by = gene]
  rate <- nrow(itv_tbl) / genome_bp
  g <- merge(gene_width, n_loci, by = "gene", all.x = TRUE)
  g[is.na(N), N := 0L]
  g[, expected := promoter_bp * rate]
  g[, (paste0("order1_poisson_neglog10p_", tag)) := -ppois(N - 1L, expected, lower.tail = FALSE, log.p = TRUE) / log(10)]
  g[, c("gene", paste0("order1_poisson_neglog10p_", tag)), with = FALSE]
}
scored <- merge(score_gene(itv_combined, "combined"), score_gene(itv_het, "heterodimer"), by = "gene")
scored <- merge(scored, gene_width, by = "gene")
fwrite(scored, file.path(out_results, "order1_background_gene_scores.csv"))

## ---- 4. validate against both ground truths, same framework throughout ----
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

feature_cols <- c("order1_poisson_neglog10p_combined", "order1_poisson_neglog10p_heterodimer")
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

cat("\n=== Order-1-background run, cross-validated AUC (5-fold) ===\n")
cat(sprintf("vs. real LIVER PPARA ChIP:      order-1 score = %.3f | promoter-size-only baseline = %.3f\n", res_chip$auc_rf, res_chip$auc_base))
cat(sprintf("vs. real CARDIAC accessibility: order-1 score = %.3f | promoter-size-only baseline = %.3f\n", res_cardiac$auc_rf, res_cardiac$auc_base))
cat("\n(For reference:\n")
cat("  Order-0 unmasked (this analysis's baseline):  0.773/0.770 vs ChIP, 0.588/0.609 vs cardiac\n")
cat("  Repeat-masked, corrected normalization:        0.547/0.707 vs ChIP, 0.531/0.611 vs cardiac)\n")

## ---- 5. how correlated are order-0 and order-1 hits, gene-for-gene? ----
order0 <- fread(file.path(out_results, "gene_scoring_size_normalized.csv"))[, .(gene, poisson_neglog10p_heterodimer)]
cmp <- merge(dat[, .(gene, order1_poisson_neglog10p_heterodimer)], order0, by = "gene")
cat(sprintf("\nSpearman correlation, order-0 vs order-1 heterodimer score (same genes, same genome, only background changed): %.3f\n",
            cor(cmp$order1_poisson_neglog10p_heterodimer, cmp$poisson_neglog10p_heterodimer, method = "spearman")))

fwrite(dat, file.path(out_results, "order1_background_validated.csv"))
message("\nDone.")
