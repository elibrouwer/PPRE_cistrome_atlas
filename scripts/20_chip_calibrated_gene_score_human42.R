# ChIP-calibrated gene score: the adipocyte-PPARG-trained elastic net (19) is applied to 500 bp promoter tiles (-3 kb to +1 kb of every TSS, step 250 bp, per-motif max relative FIMO score,
# 1 kb GC of the tile centre as covariate); gene score = max over the TSS tiles of the predicted binding probability (optionally x exp(-|d|/3000)). No ChIP signal enters the score at gene level.
# Evaluated like 16b/17 on the benchmark gold sets (PPARgene targets, mouse-heart PPARA-response sets; matched-background AUROC) against the hand-built near-TSS scores.
BENCH <- "C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Minor internship/R_code/benchmark"
source(file.path(BENCH, "config.R")); source(file.path(BENCH, "helpers_hits.R")); source(file.path(BENCH, "helpers_gold.R"))
suppressPackageStartupMessages({ library(GenomicRanges); library(glmnet) })
set.seed(11)
full_local <- "C:/Users/brouw/PPRE_cistrome_atlas/results/human42_pipeline_full_local"
out <- "C:/Users/brouw/PPRE_cistrome_atlas/results/human42_pipeline/chip_calibration"
fimo_root <- "C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Minor internship/PPAR_RXR_motifs_Hocomoco_Jaspar/FIMO_results_human42_unmasked"
gi <- load_genes(); tss <- readRDS(file.path(cache_dir, "tss_manuscript.rds")); K <- readRDS("C:/Users/brouw/PPRE_cistrome_atlas/results/human42_pipeline/distal_kernels/distal_kernel_scores.rds")
cal <- readRDS(file.path(full_local, "chip_calibration_features.rds")); W <- cal$W; M <- cal$M; motifs <- colnames(M)
ad <- which(W$group == "Adipocyte")
cv <- cv.glmnet(cbind(M, gc = W$gcf)[ad, ], W$y[ad], family = "binomial", alpha = 0.5, nfolds = 5)
cat("model non-zero coefficients:", sum(coef(cv, s = "lambda.1se")[-1] != 0), "\n")

# tiles
U <- unique(tss[in_universe == TRUE & ENSEMBL %in% gi$ENSEMBL, .(ENSEMBL, chr, strand, tss)]); U[, `:=`(g = match(ENSEMBL, gi$ENSEMBL), sgn = ifelse(strand == "+", 1L, -1L))]
U <- U[chr %in% paste0("chr", c(1:22, "X"))]
ds <- seq(-2750, 750, by = 250)
T <- U[rep(seq_len(.N), each = length(ds))]; T[, d := rep(ds, times = nrow(U))]; T[, pos := tss + sgn * d]; T <- T[pos > 250]; T[, id := .I]
cat("TSS:", nrow(U), " tiles:", nrow(T), "\n")
gt <- GRanges(T$chr, IRanges(T$pos - 250L, T$pos + 250L))
feat <- list()
for (f in list.files(fimo_root, "fimo\\.tsv$", recursive = TRUE, full.names = TRUE)) {
  d <- fread(f, sep = "\t", select = c("motif_id", "sequence_name", "start", "stop", "score"), colClasses = list(character = c("motif_id", "sequence_name")))
  d <- d[sequence_name %in% as.character(c(1:22, "X")) & !startsWith(motif_id, "#")]
  d[, `:=`(start = as.integer(start), stop = as.integer(stop), score = as.numeric(score))]; d[, rel := score / max(score), by = motif_id]
  ov <- findOverlaps(GRanges(paste0("chr", d$sequence_name), IRanges(d$start, d$stop)), gt)
  feat[[f]] <- data.table(id = subjectHits(ov), motif = d$motif_id[queryHits(ov)], rel = d$rel[queryHits(ov)])[, .(rel = max(rel)), by = .(id, motif)]
  cat(basename(dirname(f)), "done\n")
}
F <- rbindlist(feat)[, .(rel = max(rel)), by = .(id, motif)]
X <- matrix(0, nrow(T), length(motifs), dimnames = list(NULL, motifs)); X[cbind(F$id, match(F$motif, motifs))] <- F$rel
gcb <- fread("C:/Users/brouw/PPRE_cistrome_atlas/results/human42_pipeline/clustering/gc_bins_1kb.tsv"); gcb[, gcf := gc / pmax(length - n, 1)]
T[, bin := pos %/% 1000L]
T <- merge(T, gcb[, .(chr = chrom, bin, gcf)], by = c("chr", "bin"), all.x = TRUE, sort = FALSE)   # merge reorders rows -> re-align X via id
X <- X[T$id, , drop = FALSE]; T[is.na(gcf), gcf := median(gcb$gcf, na.rm = TRUE)]
lp <- as.numeric(predict(cv, cbind(X, gc = T$gcf), s = "lambda.1se"))
lp_nogc <- as.numeric(predict(cv, cbind(X, gc = median(W$gcf)), s = "lambda.1se"))
T[, `:=`(p = plogis(lp), p_nogc = plogis(lp_nogc), w = exp(-abs(d) / 3000))]; T[, `:=`(pw = p * w, p_nogc_w = p_nogc * w)]
nG <- nrow(gi)
tog <- function(col, fun) { r <- rep(0, nG); a <- T[, .(x = fun(get(col))), by = g]; r[a$g] <- a$x; r }
S <- list(chip_model_max = tog("p", max), chip_model_max_decay = tog("pw", max), chip_model_nogc_max = tog("p_nogc", max),
          chip_model_nogc_max_decay = tog("p_nogc_w", max), chip_model_top3_mean = { r <- rep(0, nG); a <- T[order(g, -p), .(x = mean(head(p, 3))), by = g]; r[a$g] <- a$x; r },
          prom_het = K$prom_het, prom_base = K$prom_base)
saveRDS(S, file.path(out, "chip_calibrated_gene_scores.rds"))
gold <- load_gold_all(gi); gold <- gold[c(PG_SETS, HEART_SETS)]
ev <- make_evaluator(gi, gold, ndraw = 300L)
R <- rbindlist(lapply(names(S), function(nm) ev(S[[nm]], nm)$table)); fwrite(R, file.path(out, "chip_calibrated_gene_score_auroc_by_goldset.csv"))
Wd <- dcast(R, score ~ gold_set, value.var = "auroc_matched")
Wd[, `:=`(mean_PPARgene = rowMeans(as.matrix(.SD[, PG_SETS, with = FALSE])), mean_heart = rowMeans(as.matrix(.SD[, HEART_SETS, with = FALSE])))]
Wd[, overall := (mean_PPARgene + mean_heart) / 2]; Wd <- Wd[order(-overall)]; fwrite(Wd, file.path(out, "chip_calibrated_gene_score_auroc_summary.csv"))
print(Wd[, .(score, mean_PPARgene = round(mean_PPARgene, 3), mean_heart = round(mean_heart, 3), overall = round(overall, 3))])
bc <- boot_compare(S, gold, gi, ND = 100L, B = 500L)
cmp <- function(a, b) { f <- function(sets) { d <- rowMeans(bc$ba[, a, sets] - bc$ba[, b, sets]); c(mean = mean(d), lo = unname(quantile(d, .025)), hi = unname(quantile(d, .975))) }
  data.table(a = a, b = b, PPARgene = round(f(PG_SETS)[1], 3), PPARgene_lo = round(f(PG_SETS)[2], 3), PPARgene_hi = round(f(PG_SETS)[3], 3),
             heart = round(f(HEART_SETS)[1], 3), heart_lo = round(f(HEART_SETS)[2], 3), heart_hi = round(f(HEART_SETS)[3], 3)) }
top <- setdiff(head(Wd$score, 3), c("prom_het", "prom_base"))
cb <- rbindlist(c(lapply(c("chip_model_max", "chip_model_max_decay", "chip_model_nogc_max"), function(a) cmp(a, "prom_het")), lapply(c("chip_model_max", "chip_model_max_decay"), function(a) cmp(a, "prom_base"))))
print(cb); fwrite(cb, file.path(out, "chip_calibrated_gene_score_paired_bootstrap.csv"))
