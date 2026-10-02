# Calibrating the 42 human-source motifs on PPARG ChIP-seq (ChIP-Atlas) as a test: do FIMO hits separate bound from GC-matched control windows,
# which motifs carry the signal, and does a model trained on adipocyte peaks transfer to other cell groups (incl. HUVEC, the only cardiovascular PPARG data)?
# Windows: 500 bp centred on peaks (MACS2 score >= 250, <= 20,000 peaks per group) vs 3 GC-matched control windows per peak (1 kb GC class of 2.5 %, N-free, not overlapping a peak).
# Features: per-motif maximum relative FIMO score (score / motif max) in the window, from the unmasked p <= 1e-4 hits of the 42 motifs.
# Models: covariate only (1 kb GC), all-motif max, heterodimer-motif max, elastic net on the 42 per-motif features (+ GC). Chromosome-blocked evaluation (train odd / test even chromosomes and vice versa) in adipocyte, then transfer to the other groups.
BENCH <- "C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Minor internship/R_code/benchmark"
suppressPackageStartupMessages({ library(data.table); library(GenomicRanges); library(glmnet) })
set.seed(20251001)
cache_dir <- file.path(BENCH, "cache"); STD <- paste0("chr", c(1:22, "X"))
out <- "C:/Users/brouw/PPRE_cistrome_atlas/results/human42_pipeline/chip_calibration"; dir.create(out, showWarnings = FALSE)
fimo_root <- "C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Minor internship/PPAR_RXR_motifs_Hocomoco_Jaspar/FIMO_results_human42_unmasked"
inv <- fread("C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Minor internship/PPAR_RXR_motifs_Hocomoco_Jaspar/Motifs_update_2026-09-30/motif_inventory_gates.csv", encoding = "UTF-8")
het <- inv[category == "PPAR:RXR heterodimer", motif_id]
gc <- fread("C:/Users/brouw/PPRE_cistrome_atlas/results/human42_pipeline/clustering/gc_bins_1kb.tsv"); gc[, `:=`(gcf = gc / (length - n), cls = floor(100 * gc / (length - n) / 2.5))]
gcok <- gc[n == 0 & length == 1000]

pk <- fread(file.path(cache_dir, "chipatlas_PPARG_min.tsv"), header = FALSE); setnames(pk, c("chr", "start", "end", "srx", "group", "cell", "score"))
pk <- pk[chr %in% STD & score >= 250]; pk[, group := gsub("%20", " ", group)]
cat("PPARG peaks (score>=250) per group:\n"); print(pk[, .N, by = group])
W <- list()
for (g in unique(pk$group)) {
  p <- pk[group == g]; p <- p[sample(.N, min(.N, 20000))]
  p[, mid := (start + end) %/% 2L]
  gr_p <- GRanges(pk[group == g]$chr, IRanges(pk[group == g]$start + 1L, pk[group == g]$end))
  pw <- p[, .(chr, pos = mid, y = 1L)]
  pw[, bin := pos %/% 1000L]; pw <- merge(pw, gc[, .(chr = chrom, bin, gcf, cls)], by = c("chr", "bin"), all.x = TRUE); pw <- pw[!is.na(cls)]
  cw <- rbindlist(lapply(split(pw, pw$cls), function(d) { pool <- gcok[cls == d$cls[1]]; if (!nrow(pool)) return(NULL)
    s <- pool[sample(.N, 3 * nrow(d), replace = TRUE)]; data.table(chr = s$chrom, pos = s$bin * 1000L + sample(250:749, nrow(s), TRUE), y = 0L, gcf = s$gcf, cls = s$cls) }))
  ok <- !overlapsAny(GRanges(cw$chr, IRanges(cw$pos - 250L, cw$pos + 250L)), gr_p); cw <- cw[ok]
  w <- rbind(pw[, .(chr, pos, y, gcf)], cw[, .(chr, pos, y, gcf)]); w[, group := g]; W[[g]] <- w
  cat(g, ": positives", nrow(pw), " controls", nrow(cw), "\n")
}
W <- rbindlist(W); W[, id := .I]
gw <- GRanges(W$chr, IRanges(W$pos - 250L, W$pos + 250L))

# per-motif max relative score per window
feat <- list()
for (f in list.files(fimo_root, "fimo\\.tsv$", recursive = TRUE, full.names = TRUE)) {
  d <- fread(f, sep = "\t", select = c("motif_id", "sequence_name", "start", "stop", "score"), colClasses = list(character = c("motif_id", "sequence_name")))
  d <- d[sequence_name %in% as.character(c(1:22, "X")) & !startsWith(motif_id, "#")]
  d[, `:=`(start = as.integer(start), stop = as.integer(stop), score = as.numeric(score))]
  d[, rel := score / max(score), by = motif_id]
  ov <- findOverlaps(GRanges(paste0("chr", d$sequence_name), IRanges(d$start, d$stop)), gw)
  a <- data.table(id = subjectHits(ov), motif = d$motif_id[queryHits(ov)], rel = d$rel[queryHits(ov)])[, .(rel = max(rel)), by = .(id, motif)]
  feat[[f]] <- a; cat(basename(dirname(f)), "done\n")
}
F <- rbindlist(feat); motifs <- sort(unique(inv$motif_id[inv$motif_id %in% F$motif]))
M <- matrix(0, nrow(W), length(motifs), dimnames = list(NULL, motifs)); M[cbind(F$id, match(F$motif, motifs))] <- F$rel
saveRDS(list(W = W, M = M), file.path("C:/Users/brouw/PPRE_cistrome_atlas/results/human42_pipeline_full_local", "chip_calibration_features.rds"))

auc <- function(x, y) { r <- rank(x); n1 <- sum(y == 1); n0 <- sum(y == 0); (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
hetcols <- intersect(het, motifs)
fitpred <- function(tr, te) {
  X <- cbind(M, gc = W$gcf)
  cv <- cv.glmnet(X[tr, ], W$y[tr], family = "binomial", alpha = 0.5, nfolds = 5)
  list(pred = as.numeric(predict(cv, X[te, ], s = "lambda.1se")), model = cv)
}
odd <- W$chr %in% paste0("chr", c(seq(1, 21, 2), "X")); ad <- W$group == "Adipocyte"
rows <- list()
for (dir in 1:2) {
  tr <- which(ad & (if (dir == 1) odd else !odd)); te <- which(ad & (if (dir == 1) !odd else odd))
  fp <- fitpred(tr, te)
  rows[[length(rows) + 1]] <- data.table(fold = dir, test = "Adipocyte (held-out chromosomes)", model = c("GC covariate only", "all-motif max", "heterodimer-motif max", "elastic net (42 motifs + GC)"),
    AUROC = c(auc(W$gcf[te], W$y[te]), auc(apply(M[te, ], 1, max), W$y[te]), auc(apply(M[te, hetcols], 1, max), W$y[te]), auc(fp$pred, W$y[te])))
}
full <- fitpred(which(ad), which(!ad))   # train on all adipocyte, predict everything else
for (g in setdiff(unique(W$group), "Adipocyte")) { te <- which(W$group == g); p <- full$pred[match(te, which(!ad))]
  rows[[length(rows) + 1]] <- data.table(fold = NA, test = paste0(g, " (transfer)"), model = c("GC covariate only", "all-motif max", "heterodimer-motif max", "elastic net (42 motifs + GC)"),
    AUROC = c(auc(W$gcf[te], W$y[te]), auc(apply(M[te, ], 1, max), W$y[te]), auc(apply(M[te, hetcols], 1, max), W$y[te]), auc(p, W$y[te]))) }
R <- rbindlist(rows); fwrite(R, file.path(out, "ppar_chip_calibration_auroc.csv"))
print(dcast(R[, .(AUROC = mean(AUROC)), by = .(test, model)], test ~ model, value.var = "AUROC")[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])
cf <- as.matrix(coef(full$model, s = "lambda.1se")); cf <- data.table(feature = rownames(cf), coef = cf[, 1])[feature != "(Intercept)" & coef != 0][order(-coef)]
cf <- merge(cf, inv[, .(feature = motif_id, name, category, architecture)], by = "feature", all.x = TRUE)[order(-coef)]; fwrite(cf, file.path(out, "ppar_chip_calibration_coefficients.csv"))
cat("\nnon-zero coefficients (adipocyte model):\n"); print(cf)
