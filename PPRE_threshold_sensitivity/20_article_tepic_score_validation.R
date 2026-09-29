# Validate the article's own TEPIC-style gene score (19) with the same rigor
# used throughout this analysis: cross-validated random forest AUC vs. a
# trivial "how many transcripts does this gene have nearby" baseline, against
# BOTH real liver PPARA ChIP-seq (07/15/16) and real cardiac differentially-
# accessible chromatin (18) -- so it is directly comparable to every other
# method already tested.

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(ranger)
  library(pROC)
  library(rtracklayer)
})
source("./PPRE_threshold_sensitivity/00_config.R")

tepic <- fread(file.path(out_results, "article_tepic_style_gene_scores.csv"))

## ---- rebuild per-gene transcript count (the "size" baseline for this per-transcript scheme) ----
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
tx_by_gene <- transcriptsBy(txdb, by = "gene")
expr_genes <- fread(expr_nonzero_csv)[, unique(HGNC_gene_name)]
entrez_map <- AnnotationDbi::select(org.Hs.eg.db, keys = expr_genes, keytype = "SYMBOL", columns = "ENTREZID")
entrez_map <- entrez_map[!is.na(entrez_map$ENTREZID) & entrez_map$ENTREZID %in% names(tx_by_gene), ]
entrez_map <- entrez_map[!duplicated(entrez_map$SYMBOL), ]
n_tx <- lengths(tx_by_gene[entrez_map$ENTREZID])
tx_count <- data.table(gene = entrez_map$SYMBOL, n_transcripts = as.integer(n_tx))

dat <- merge(tepic, tx_count, by = "gene", all.x = TRUE)
dat[is.na(n_transcripts), n_transcripts := 1L]

## ---- real ChIP + cardiac labels (same construction as 15 / 18, using per-transcript
## TSS windows here for consistency with how the TEPIC score itself was built) ----
tx_all <- unlist(tx_by_gene[entrez_map$ENTREZID])
gene_lookup <- setNames(entrez_map$SYMBOL, entrez_map$ENTREZID)
mcols(tx_all)$gene <- gene_lookup[names(tx_all)]
tss <- ifelse(as.character(strand(tx_all)) == "-", end(tx_all), start(tx_all))
promo_tx <- GRanges(seqnames(tx_all), IRanges(tss - 3000L, tss + 1000L), gene = mcols(tx_all)$gene)

# liver PPARA ChIP
raw_chip <- fread(file.path(repo_dir, "PPARA_chipseq_overlap", "chipatlas_liver_peaks_raw.bed"),
                   header = FALSE, col.names = c("chrom","start","end"))
std_chroms <- paste0("chr", c(1:22,"X","Y","M"))
raw_chip <- raw_chip[chrom %in% std_chroms]
chip_gr <- reduce(GRanges(raw_chip$chrom, IRanges(raw_chip$start+1L, raw_chip$end)))
chip_genes <- unique(promo_tx$gene[queryHits(findOverlaps(promo_tx, chip_gr, ignore.strand = TRUE))])

# real cardiac HCM+PLN differentially-accessible regions
hcm <- readxl::read_xlsx(file.path(base_dir, "HCM_differentially_regions.xlsx"), sheet = "Supplementary Table 2A")
gr_hcm <- GRanges(hcm$`Chromosome number`, IRanges(hcm$`Region start`, hcm$`Region end`))
pln <- readxl::read_xlsx(file.path(base_dir, "PLN_differentially_regions.xlsx"), skip = 1)
pln <- pln[!is.na(pln$`Region start†`), ]
gr_pln <- GRanges(pln$Chromosome, IRanges(pln$`Region start†`, pln$`Region end†`))
chain <- import.chain(file.path(base_dir, "hg19ToHg38.over.chain"))
cardiac_gr <- reduce(c(unlist(liftOver(gr_hcm, chain)), unlist(liftOver(gr_pln, chain))))
seqlevelsStyle(cardiac_gr) <- "UCSC"
cardiac_genes <- unique(promo_tx$gene[queryHits(findOverlaps(promo_tx, cardiac_gr, ignore.strand = TRUE))])

dat[, chip_supported := as.integer(gene %in% chip_genes)]
dat[, cardiac_supported := as.integer(gene %in% cardiac_genes)]
cat(sprintf("Liver ChIP-supported: %d/%d (%.1f%%); Cardiac-accessibility-supported: %d/%d (%.1f%%)\n",
            sum(dat$chip_supported), nrow(dat), 100*mean(dat$chip_supported),
            sum(dat$cardiac_supported), nrow(dat), 100*mean(dat$cardiac_supported)))

## ---- cross-validated AUC: TEPIC score(s) vs. n_transcripts-only baseline, for both labels ----
feature_cols <- c("tepic_sum","tepic_mean","tepic_median","tepic_max","n_sites","n_transcripts")
set.seed(42); K <- 5
dat[, fold := sample(rep(1:K, length.out = .N))]

run_cv <- function(label_col) {
  dat[, pred_rf := NA_real_]; dat[, pred_base := NA_real_]
  for (k in 1:K) {
    train <- dat[fold != k]; test_idx <- dat[, which(fold == k)]
    rf <- ranger(x = train[, ..feature_cols], y = factor(train[[label_col]]), probability = TRUE, num.trees = 500, seed = 42)
    dat[test_idx, pred_rf := predict(rf, data = dat[test_idx, ..feature_cols])$predictions[, "1"]]
    base_model <- glm(as.formula(paste(label_col, "~ n_transcripts")), data = train, family = binomial())
    dat[test_idx, pred_base := predict(base_model, newdata = dat[test_idx], type = "response")]
  }
  list(auc_rf = as.numeric(auc(dat[[label_col]], dat$pred_rf, quiet = TRUE)),
       auc_base = as.numeric(auc(dat[[label_col]], dat$pred_base, quiet = TRUE)))
}

res_chip <- run_cv("chip_supported")
res_cardiac <- run_cv("cardiac_supported")

cat("\n=== Article's TEPIC-style score, cross-validated AUC (5-fold) ===\n")
cat(sprintf("vs. real LIVER PPARA ChIP:        RF (all TEPIC features) = %.3f | n_transcripts-only baseline = %.3f\n", res_chip$auc_rf, res_chip$auc_base))
cat(sprintf("vs. real CARDIAC accessibility:    RF (all TEPIC features) = %.3f | n_transcripts-only baseline = %.3f\n", res_cardiac$auc_rf, res_cardiac$auc_base))
cat("(For reference, this analysis's own size-normalized Poisson score, 15/18, scored:\n")
cat("  vs. liver ChIP:     0.773 vs 0.770 baseline\n  vs. cardiac access.: 0.588 vs 0.609 baseline)\n")

## ---- univariate: each of the 4 TEPIC summary variants alone ----
univariate <- data.table(variant = c("tepic_sum","tepic_mean","tepic_median","tepic_max"))
univariate[, auc_vs_liver_chip := sapply(variant, function(v) as.numeric(auc(dat$chip_supported, dat[[v]], quiet=TRUE)))]
univariate[, auc_vs_cardiac := sapply(variant, function(v) as.numeric(auc(dat$cardiac_supported, dat[[v]], quiet=TRUE)))]
cat("\n=== Each TEPIC summary variant alone (no CV, in-sample AUC -- upper bound) ===\n")
print(univariate)

fwrite(dat[, .(gene, tepic_sum, tepic_mean, tepic_median, tepic_max, n_sites, n_transcripts, chip_supported, cardiac_supported, pred_rf)],
       file.path(out_results, "article_tepic_score_validated.csv"))
fwrite(univariate, file.path(out_results, "article_tepic_univariate_auc.csv"))
message("\nDone.")
