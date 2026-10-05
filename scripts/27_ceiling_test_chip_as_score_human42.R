# Ceiling test: how well does REAL ChIP-seq binding predict the gold sets? Gene score = strongest ChIP peak (MACS score) in the promoter window (-3 kb/+1 kb) or within +-10 kb of the TSS (0 if none).
# If real binding only reaches AUROC about 0.6 on PPARgene targets and mouse-heart PPARA-response genes, then the motif scores (about 0.55-0.6) are near what these gold standards allow.
# ChIP sets: ChIP-Atlas PPARA / PPARG / RXRA by cell group (score >= 250) and PPARD in human myofibroblasts (E-MTAB-371, lifted hg19 -> hg38).
BENCH <- "C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Minor internship/R_code/benchmark"
source(file.path(BENCH, "config.R")); source(file.path(BENCH, "helpers_hits.R")); source(file.path(BENCH, "helpers_gold.R"))
suppressPackageStartupMessages({ library(GenomicRanges); library(rtracklayer) })
out <- "C:/Users/brouw/PPRE_cistrome_atlas/results/human42_pipeline/ceiling_test"; dir.create(out, showWarnings = FALSE)
gi <- load_genes(); tss <- readRDS(file.path(cache_dir, "tss_manuscript.rds")); u <- tss[in_universe == TRUE]
mk <- function(up, dn) { w <- GRanges(u$chr, IRanges(ifelse(u$strand == "+", u$tss - up, u$tss - dn), ifelse(u$strand == "+", u$tss + dn, u$tss + up))); w <- trim(suppressWarnings(w)); w$ENSEMBL <- u$ENSEMBL; w }
Wn <- list(prom = mk(3000, 1000), `10kb` = mk(10000, 10000))
rd <- function(f) { d <- fread(file.path(cache_dir, f), header = FALSE); setnames(d, c("chr", "start", "end", "srx", "group", "cell", "score")); d[chr %in% STD_CHR & score >= 250][, group := gsub("%20", " ", group)] }
sets <- list()
for (tf in c("PPARA", "PPARG", "RXRA")) { d <- rd(paste0("chipatlas_", tf, "_min.tsv")); for (g in unique(d$group)) { x <- d[group == g]; if (nrow(x) >= 1000) sets[[paste(tf, g)]] <- GRanges(x$chr, IRanges(x$start + 1L, x$end), score = x$score) } }
pk <- fread(file.path(BENCH, "external_data/PPARD_myofibroblast/PPAR_DMSO_peaks.tsv"), check.names = FALSE)
g19 <- GRanges(pk$chr, IRanges(pk$start + 1L, pk$end), score = pk$`-log(p-value)`); lf <- liftOver(g19, import.chain(chain_file)); g38 <- unlist(lf[lengths(lf) == 1]); sets[["PPARD myofibroblast"]] <- g38[seqnames(g38) %in% STD_CHR]
chipscore <- function(gr, w) { ov <- findOverlaps(w, gr); a <- data.table(g = match(w$ENSEMBL[queryHits(ov)], gi$ENSEMBL), s = mcols(gr)$score[subjectHits(ov)])[!is.na(g), .(x = max(s)), by = g]; r <- rep(0, nrow(gi)); r[a$g] <- a$x; r }
S <- list(); for (n in names(sets)) for (wn in names(Wn)) S[[paste(n, wn, sep = " | ")]] <- chipscore(sets[[n]], Wn[[wn]])
Sub <- readRDS("C:/Users/brouw/PPRE_cistrome_atlas/results/human42_pipeline/distal_kernels/motif_subset_scores.rds")
S$`MOTIF: PPAR-containing 14` <- Sub$PPAR_containing_14; S$`MOTIF: heterodimer 3` <- Sub$heterodimer_3; S$`MOTIF: all 42` <- Sub$all_42
gold <- load_gold_all(gi); gold <- gold[c(PG_SETS, HEART_SETS)]
ev <- NULL; for (s in 1:30) { set.seed(1100 + s); ev <- tryCatch(make_evaluator(gi, gold, ndraw = 300L), error = function(e) NULL); if (!is.null(ev)) break }
R <- rbindlist(lapply(names(S), function(nm) ev(S[[nm]], nm)$table)); fwrite(R, file.path(out, "ceiling_test_auroc_by_goldset.csv"))
W <- dcast(R, score ~ gold_set, value.var = "auroc_matched"); W[, `:=`(mean_PPARgene = rowMeans(as.matrix(.SD[, PG_SETS, with = FALSE])), mean_heart = rowMeans(as.matrix(.SD[, HEART_SETS, with = FALSE])))]; W[, overall := (mean_PPARgene + mean_heart) / 2]
W <- W[order(-overall)]; fwrite(W, file.path(out, "ceiling_test_auroc_summary.csv"))
print(W[, .(score, PPARgene = round(mean_PPARgene, 3), heart = round(mean_heart, 3), overall = round(overall, 3))], nrows = 60)
