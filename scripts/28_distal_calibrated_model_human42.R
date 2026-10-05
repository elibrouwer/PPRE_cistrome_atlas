# Idea 3: distal gene scores from the ChIP-calibrated model. The elastic net trained on adipocyte PPARG peaks (19/20) is applied genome-wide to 500 bp bins
# (linear predictor = intercept + sum_m coef_m x best relative score of motif m in the bin + coef_GC x GC of the 1 kb window); p = plogis(.).
# Gene score = max over bins within a window of p x distance kernel. Variants: promoter only (-3k/+1k, tau 3 kb), +-100 kb with tau 10 / 20 / 50 kb, nearest-TSS exclusive (tau 20 kb), without the GC term, and rank averages.
# Evaluated on the gold registry (27a), in particular the distal-peak sets (family D) that a promoter score cannot see.
BENCH <- "C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Minor internship/R_code/benchmark"
source(file.path(BENCH, "config.R")); source(file.path(BENCH, "helpers_hits.R")); source(file.path(BENCH, "helpers_gold.R"))
suppressPackageStartupMessages(library(glmnet))
set.seed(28)
full_local <- "C:/Users/brouw/PPRE_cistrome_atlas/results/human42_pipeline_full_local"; out <- "C:/Users/brouw/PPRE_cistrome_atlas/results/human42_pipeline/score_exploration"
fimo_root <- "C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Minor internship/PPAR_RXR_motifs_Hocomoco_Jaspar/FIMO_results_human42_unmasked"
cal <- readRDS(file.path(full_local, "chip_calibration_features.rds")); W <- cal$W; M <- cal$M; ad <- which(W$group == "Adipocyte")
cv <- cv.glmnet(cbind(M, gc = W$gcf)[ad, ], W$y[ad], family = "binomial", alpha = 0.5, nfolds = 5)
cf <- as.matrix(coef(cv, s = "lambda.1se")); b0 <- cf[1, 1]; cgc <- cf["gc", 1]; cm <- cf[setdiff(rownames(cf), c("(Intercept)", "gc")), 1]; cm <- cm[cm != 0]
cat("non-zero motif coefficients:", length(cm), " GC coef:", round(cgc, 3), " intercept:", round(b0, 3), "\n")
acc <- list()
for (f in list.files(fimo_root, "fimo\\.tsv$", recursive = TRUE, full.names = TRUE)) {
  d <- fread(f, sep = "\t", select = c("motif_id", "sequence_name", "start", "stop", "score"), colClasses = list(character = c("motif_id", "sequence_name")))
  d <- d[sequence_name %in% as.character(c(1:22, "X")) & !startsWith(motif_id, "#")]; d[, `:=`(start = as.integer(start), stop = as.integer(stop), score = as.numeric(score))]
  d[, rel := score / max(score), by = motif_id]; d <- d[motif_id %in% names(cm)]; if (!nrow(d)) next
  d[, `:=`(chr = paste0("chr", sequence_name), bin = ((start + stop) %/% 2L) %/% 500L)]
  a <- d[, .(rel = max(rel)), by = .(chr, bin, motif_id)]; a[, contrib := cm[motif_id] * rel]; acc[[f]] <- a[, .(contrib = sum(contrib)), by = .(chr, bin)]; cat(basename(dirname(f)), "done\n") }
Bn <- rbindlist(acc)[, .(contrib = sum(contrib)), by = .(chr, bin)]; Bn[, pos := bin * 500L + 250L]
gcb <- fread("C:/Users/brouw/PPRE_cistrome_atlas/results/human42_pipeline/clustering/gc_bins_1kb.tsv"); gcb[, gcf := gc / pmax(length - n, 1)]
Bn[, bin1k := pos %/% 1000L]; Bn <- merge(Bn, gcb[, .(chr = chrom, bin1k = bin, gcf)], by = c("chr", "bin1k"), all.x = TRUE); Bn[is.na(gcf), gcf := median(gcb$gcf, na.rm = TRUE)]
Bn[, `:=`(p = plogis(b0 + contrib + cgc * gcf), p_nogc = plogis(b0 + contrib + cgc * median(W$gcf)))]
cat("bins with a model-relevant hit:", nrow(Bn), "  p quantiles:", round(quantile(Bn$p, c(.5, .9, .99, .999)), 3), "\n")
saveRDS(Bn[, .(chr, pos, p, p_nogc)], file.path(full_local, "calibrated_bins500.rds"))

gi <- load_genes(); tss <- readRDS(file.path(cache_dir, "tss_manuscript.rds")); nG <- nrow(gi)
U <- unique(tss[in_universe == TRUE & ENSEMBL %in% gi$ENSEMBL, .(ENSEMBL, chr, strand, tss)]); U[, `:=`(g = match(ENSEMBL, gi$ENSEMBL), lo = tss - 100000L, hi = tss + 100000L, sgn = ifelse(strand == "+", 1L, -1L))]
J <- Bn[U, on = .(chr, pos >= lo, pos <= hi), nomatch = 0L, allow.cartesian = TRUE, .(g = i.g, d = (x.pos - i.tss) * i.sgn, p, p_nogc)]
cat("joined rows:", nrow(J), "\n")
gagg <- function(D, val, fun = "max") { D[, vv := val]; a <- if (fun == "max") D[, .(x = max(vv)), by = g] else D[, .(x = sum(vv)), by = g]; D[, vv := NULL]; r <- rep(0, nG); r[a$g] <- a$x; r }   # val is stored as a column so it is grouped with the rows (a global vector inside j would be recycled over the whole table)
sc <- function(D, col, k) gagg(D, D[[col]] * k, "max")
P <- J[d >= -3000 & d <= 1000]; kP <- exp(-abs(P$d) / 3000)
S <- list(cal_prom = sc(P, "p", kP), cal_prom_nogc = sc(P, "p_nogc", kP))
for (tau in c(10000, 20000, 50000)) { k <- exp(-abs(J$d) / tau); S[[sprintf("cal_wide_tau%dk", tau / 1000)]] <- sc(J, "p", k); S[[sprintf("cal_wide_nogc_tau%dk", tau / 1000)]] <- sc(J, "p_nogc", k) }
top3 <- function(D, k) { D[, v := p * k]; r <- rep(0, nG); a <- D[order(g, -v), .(x = mean(head(v, 3))), by = g]; r[a$g] <- a$x; D[, v := NULL]; r }
S$cal_wide_top3_tau20k <- top3(J, exp(-abs(J$d) / 20000))
Tall <- unique(tss[, .(ENSEMBL, chr, strand, tss)]); Tall[, `:=`(g = match(ENSEMBL, gi$ENSEMBL), sgn = ifelse(strand == "+", 1L, -1L))]; setorder(Tall, chr, tss)
nearest <- function(chr_, pos_) { T <- Tall[chr == chr_]; i <- findInterval(pos_, T$tss); i0 <- pmax(i, 1L); i1 <- pmin(i + 1L, nrow(T)); pick <- ifelse(abs(pos_ - T$tss[i0]) <= abs(T$tss[i1] - pos_), i0, i1); list(g = T$g[pick], d = (pos_ - T$tss[pick]) * T$sgn[pick]) }
E <- rbindlist(lapply(unique(Bn$chr), function(ch) { b <- Bn[chr == ch]; n <- nearest(ch, b$pos); data.table(g = n$g, d = n$d, p = b$p) })); E <- E[!is.na(g) & abs(d) <= 100000]
S$cal_excl_tau20k <- sc(E, "p", exp(-abs(E$d) / 20000))
S$cal_rankavg_prom_wide20 <- (frank(S$cal_prom) + frank(S$cal_wide_tau20k)) / 2
saveRDS(S, file.path(out, "distal_calibrated_scores.rds"))
reg <- readRDS(file.path(out, "gold_registry.rds")); reg <- reg[lengths(reg) >= 15]
Sub <- readRDS("C:/Users/brouw/PPRE_cistrome_atlas/results/human42_pipeline/distal_kernels/motif_subset_scores.rds"); S$ref_PPAR14_prom <- Sub$PPAR_containing_14; S$ref_het3_prom <- Sub$heterodimer_3; S$ref_all42_prom <- Sub$all_42
ev <- NULL; for (s in 1:40) { set.seed(1500 + s); ev <- tryCatch(make_evaluator(gi, reg, ndraw = 200L), error = function(e) NULL); if (!is.null(ev)) break }; stopifnot(!is.null(ev))
R <- rbindlist(lapply(names(S), function(nm) ev(S[[nm]], nm)$table)); R[, family := substr(gold_set, 1, 1)]; fwrite(R, file.path(out, "distal_calibrated_auroc_by_goldset.csv"))
fam <- R[, .(mean_AUROC = mean(auroc_matched)), by = .(score, family)]; Wf <- dcast(fam, score ~ family, value.var = "mean_AUROC"); Wf[, overall := rowMeans(.SD), .SDcols = c("A", "B", "C", "D", "E")]
setnames(Wf, c("A", "B", "C", "D", "E"), c("A literature", "B PPARgene+heart", "C promoter-peak", "D distal-peak", "E PPARD induced"))
print(Wf[order(-overall), lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)], nrows = 40); fwrite(Wf, file.path(out, "distal_calibrated_family_means.csv"))
