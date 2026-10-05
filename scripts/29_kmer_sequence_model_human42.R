# Idea 6 (exploration): a sequence model on top of the motif hits. Canonical 6-mer composition of the 500 bp window (2,080 features, elastic net) trained on PPARG ChIP-seq peaks vs GC-matched controls
# (gkm-SVM-like; no deep-learning library is installed, a CNN would need PyTorch). Questions:
#  (1) peak level: does a 6-mer model beat the 42 motif features (and does motif + 6-mer beat both) on held-out chromosomes (adipocyte) and in transfer (other PPARG cell groups, PPARD myofibroblast)?
#  (2) gene level: applied to -3 kb/+1 kb promoter tiles, does the 6-mer gene score beat the motif-based gene scores on the gold registry (27a)?
BENCH <- "C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Minor internship/R_code/benchmark"
source(file.path(BENCH, "config.R")); source(file.path(BENCH, "helpers_hits.R")); source(file.path(BENCH, "helpers_gold.R"))
suppressPackageStartupMessages({ library(glmnet); library(Matrix); library(Biostrings); library(BSgenome.Hsapiens.UCSC.hg38); library(GenomicRanges); library(rtracklayer) })
set.seed(29)
full_local <- "C:/Users/brouw/PPRE_cistrome_atlas/results/human42_pipeline_full_local"; out <- "C:/Users/brouw/PPRE_cistrome_atlas/results/human42_pipeline/score_exploration"
K <- 6L
all6 <- mkAllStrings(c("A", "C", "G", "T"), K); rc6 <- as.character(reverseComplement(DNAStringSet(all6))); canon <- pmin(all6, rc6); ucan <- sort(unique(canon)); cmap <- match(canon, ucan)
fold <- sparseMatrix(i = seq_along(all6), j = cmap, x = 1, dims = c(length(all6), length(ucan)))
kfeat <- function(gr) {                        # canonical 6-mer frequencies of the 500 bp windows
  seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38, gr); f <- oligonucleotideFrequency(seqs, K, as.prob = FALSE); n <- rowSums(f)
  X <- as(Matrix(f, sparse = TRUE), "CsparseMatrix") %*% fold; X <- Diagonal(x = 1 / pmax(n, 1)) %*% X; colnames(X) <- ucan; X * 500 }   # counts per 500 bp
win <- function(chr, pos) GRanges(chr, IRanges(pmax(pos - 250L, 1L), pos + 249L))
auc <- function(x, y) { r <- rank(x); n1 <- sum(y == 1); n0 <- sum(y == 0); (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0) }

# ---- (1) peak level ----
cal <- readRDS(file.path(full_local, "chip_calibration_features.rds")); W <- cal$W; M <- cal$M
ok <- W$chr %in% seqnames(BSgenome.Hsapiens.UCSC.hg38)
Xk <- kfeat(win(W$chr, W$pos)); cat("6-mer matrix:", dim(Xk), "\n")
Xm <- cbind(M, gc = W$gcf)
fit <- function(X, y, tr) { cv <- cv.glmnet(X[tr, ], y[tr], family = "binomial", alpha = 0.5, nfolds = 3); cv }
odd <- W$chr %in% paste0("chr", c(seq(1, 21, 2), "X")); ad <- W$group == "Adipocyte"
mods <- list(`motif (42 + GC)` = Xm, `6-mer` = Xk, `motif + 6-mer` = cbind(Xm, Xk))
rows <- list()
for (dir in 1:2) { tr <- which(ad & (if (dir == 1) odd else !odd)); te <- which(ad & (if (dir == 1) !odd else odd))
  for (nm in names(mods)) { cvm <- fit(mods[[nm]], W$y, tr); rows[[length(rows) + 1]] <- data.table(test = "Adipocyte (held-out chromosomes)", model = nm, fold = dir, AUROC = auc(as.numeric(predict(cvm, mods[[nm]][te, ], s = "lambda.1se")), W$y[te])) }
  cat("fold", dir, "done\n") }
full <- lapply(mods, function(X) fit(X, W$y, which(ad)))
for (g in setdiff(unique(W$group), "Adipocyte")) { te <- which(W$group == g)
  for (nm in names(mods)) rows[[length(rows) + 1]] <- data.table(test = paste0(g, " (transfer)"), model = nm, fold = NA, AUROC = auc(as.numeric(predict(full[[nm]], mods[[nm]][te, ], s = "lambda.1se")), W$y[te])) }
# PPARD myofibroblast peaks (hg19 -> hg38), GC-matched controls
pk <- fread(file.path(BENCH, "external_data/PPARD_myofibroblast/PPAR_DMSO_peaks.tsv"), check.names = FALSE)
g19 <- GRanges(pk$chr, IRanges(pk$start + 1L, pk$end)); lf <- liftOver(g19, import.chain(chain_file)); gr <- unlist(lf[lengths(lf) == 1]); gr <- gr[seqnames(gr) %in% paste0("chr", c(1:22, "X"))]
gc <- fread("C:/Users/brouw/PPRE_cistrome_atlas/results/human42_pipeline/clustering/gc_bins_1kb.tsv"); gc[, `:=`(gcf = gc / (length - n), cls = floor(100 * gc / (length - n) / 2.5))]; gcok <- gc[n == 0 & length == 1000]
mid <- (start(gr) + end(gr)) %/% 2L; pw <- data.table(chr = as.character(seqnames(gr)), pos = mid, y = 1L); pw[, bin := pos %/% 1000L]; pw <- merge(pw, gc[, .(chr = chrom, bin, gcf, cls)], by = c("chr", "bin"), all.x = TRUE); pw <- pw[!is.na(cls)]
cw <- rbindlist(lapply(split(pw, pw$cls), function(d) { pool <- gcok[cls == d$cls[1]]; if (!nrow(pool)) return(NULL); s <- pool[sample(.N, 3 * nrow(d), replace = TRUE)]; data.table(chr = s$chrom, pos = s$bin * 1000L + sample(250:749, nrow(s), TRUE), y = 0L) }))
cw <- cw[!overlapsAny(GRanges(cw$chr, IRanges(cw$pos - 250L, cw$pos + 250L)), gr)]; Wd <- rbind(pw[, .(chr, pos, y)], cw)
Xd <- kfeat(win(Wd$chr, Wd$pos)); rows[[length(rows) + 1]] <- data.table(test = "PPARD myofibroblast (transfer)", model = "6-mer", fold = NA, AUROC = auc(as.numeric(predict(full[["6-mer"]], Xd, s = "lambda.1se")), Wd$y))
R <- rbindlist(rows); fwrite(R, file.path(out, "kmer_model_peak_level_auroc.csv"))
print(dcast(R[, .(AUROC = mean(AUROC)), by = .(test, model)], test ~ model, value.var = "AUROC")[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])
cf <- as.matrix(coef(full[["6-mer"]], s = "lambda.1se")); cf <- data.table(kmer = rownames(cf), coef = cf[, 1])[kmer != "(Intercept)" & coef != 0][order(-coef)]
cat("non-zero 6-mers:", nrow(cf), "\ntop positive:\n"); print(head(cf, 15)); cat("top negative:\n"); print(tail(cf, 8)); fwrite(cf, file.path(out, "kmer_model_coefficients.csv"))

# ---- (2) gene level: promoter tiles ----
gi <- load_genes(); tss <- readRDS(file.path(cache_dir, "tss_manuscript.rds"))
U <- unique(tss[in_universe == TRUE & ENSEMBL %in% gi$ENSEMBL, .(ENSEMBL, chr, strand, tss)]); U <- U[chr %in% paste0("chr", c(1:22, "X"))]; U[, `:=`(g = match(ENSEMBL, gi$ENSEMBL), sgn = ifelse(strand == "+", 1L, -1L))]
ds <- seq(-2750, 750, by = 250); T <- U[rep(seq_len(.N), each = length(ds))]; T[, d := rep(ds, times = nrow(U))]; T[, pos := tss + sgn * d]; T <- T[pos > 500]
cat("tiles:", nrow(T), "\n"); T[, lp := NA_real_]
chunks <- split(seq_len(nrow(T)), ceiling(seq_len(nrow(T)) / 100000)); k6 <- full[["6-mer"]]
for (ci in seq_along(chunks)) { i <- chunks[[ci]]; Xt <- kfeat(win(T$chr[i], T$pos[i])); set(T, i, "lp", as.numeric(predict(k6, Xt, s = "lambda.1se", type = "link"))); cat("chunk", ci, "of", length(chunks), "\n") }
T[, `:=`(p = plogis(lp), w = exp(-abs(d) / 3000))]; nG <- nrow(gi)
tog <- function(col, fun) { r <- rep(0, nG); a <- T[, .(x = fun(get(col))), by = g]; r[a$g] <- a$x; r }
T[, pw := p * w]
S <- list(kmer_max = tog("p", max), kmer_max_decay = tog("pw", max), kmer_top3 = { r <- rep(0, nG); a <- T[order(g, -p), .(x = mean(head(p, 3))), by = g]; r[a$g] <- a$x; r })
saveRDS(S, file.path(out, "kmer_gene_scores.rds"))
Sub <- readRDS("C:/Users/brouw/PPRE_cistrome_atlas/results/human42_pipeline/distal_kernels/motif_subset_scores.rds"); C <- readRDS("C:/Users/brouw/PPRE_cistrome_atlas/results/human42_pipeline/chip_calibration/chip_calibrated_gene_scores.rds")
S$ref_PPAR14 <- Sub$PPAR_containing_14; S$ref_het3 <- Sub$heterodimer_3; S$ref_chip_top3 <- C$chip_model_top3_mean; S$ref_all42 <- Sub$all_42
reg <- readRDS(file.path(out, "gold_registry.rds")); reg <- reg[lengths(reg) >= 15]
ev <- NULL; for (s in 1:40) { set.seed(1700 + s); ev <- tryCatch(make_evaluator(gi, reg, ndraw = 200L), error = function(e) NULL); if (!is.null(ev)) break }; stopifnot(!is.null(ev))
Rg <- rbindlist(lapply(names(S), function(nm) ev(S[[nm]], nm)$table)); Rg[, family := substr(gold_set, 1, 1)]; fwrite(Rg, file.path(out, "kmer_gene_level_auroc_by_goldset.csv"))
fam <- Rg[, .(mean_AUROC = mean(auroc_matched)), by = .(score, family)]; Wf <- dcast(fam, score ~ family, value.var = "mean_AUROC"); Wf[, overall := rowMeans(.SD), .SDcols = c("A", "B", "C", "D", "E")]
setnames(Wf, c("A", "B", "C", "D", "E"), c("A literature", "B PPARgene+heart", "C promoter-peak", "D distal-peak", "E PPARD induced"))
print(Wf[order(-overall), lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)]); fwrite(Wf, file.path(out, "kmer_gene_level_family_means.csv"))
