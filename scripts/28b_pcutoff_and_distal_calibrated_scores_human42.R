# Two ideas to raise the AUROC of the promoter score (42 motifs):
#  (1) stricter hit filter: only hits with p <= 1e-4 / 3e-5 / 1e-5 / 1e-6 (proxy for a shuffled-motif FDR cut), for all 42 / PPAR-containing 14 / heterodimer 3 motifs, -3 kb/+1 kb window, exp(-|d|/3000), max of the relative FIMO score;
#  (2) distal ChIP-calibrated model: adipocyte-PPARG elastic net applied to 500 bp bins genome-wide; bin score = predicted binding probability above the baseline probability (excess, >= 0);
#      gene score = max or sum over bins within +-100 kb of any TSS, weighted exp(-|d|/tau), tau = 3, 10, 20, 50 kb; plus a promoter-only version (-3 kb/+1 kb, tau 3 kb).
# Evaluated on five gold panels (gold_panels_human42.R): PPARgene, mouse-heart, PPARG promoter-peak genes, PPARD promoter-peak genes, literature PPRE genes (matched-background AUROC).
BENCH <- "C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Minor internship/R_code/benchmark"
source(file.path(BENCH, "config.R")); source(file.path(BENCH, "helpers_hits.R")); source(file.path(BENCH, "helpers_gold.R"))
suppressPackageStartupMessages({ library(GenomicRanges); library(rtracklayer) })
source("C:/Users/brouw/PPRE_cistrome_atlas/scripts/gold_panels_human42.R")
out <- "C:/Users/brouw/PPRE_cistrome_atlas/results/human42_pipeline/score_improvement"; dir.create(out, showWarnings = FALSE)
full_local <- "C:/Users/brouw/PPRE_cistrome_atlas/results/human42_pipeline_full_local"
gi <- load_genes(); tss <- readRDS(file.path(cache_dir, "tss_manuscript.rds"))
B <- readRDS(file.path(full_local, "fimo_human42_bins500_pcuts_lp.rds")); B[, bpos := bin * 500L + 250L]; B[, pos := bpos]; mc <- readRDS(file.path(full_local, "adipocyte_model_coefs.rds"))
p0 <- plogis(mc$b0 + mc$bgc * mc$gc_med); B[, ex := pmax(plogis(mc$b0 + mc$bgc * mc$gc_med + lp) - p0, 0)]
cat("bins:", nrow(B), " bins with excess > 0:", sum(B$ex > 0), " max excess:", round(max(B$ex), 3), "\n")
U <- unique(tss[in_universe == TRUE & ENSEMBL %in% gi$ENSEMBL, .(ENSEMBL, chr, strand, tss)]); U[, `:=`(g = match(ENSEMBL, gi$ENSEMBL), sgn = ifelse(strand == "+", 1L, -1L))]
nG <- nrow(gi); S <- list()
vec <- function(a) { r <- rep(0, nG); r[a$g] <- a$x; r }
# (1) promoter window, all p cut-offs and subsets
U3 <- copy(U)[, `:=`(lo = tss - 3000L, hi = tss + 3000L)]
JP <- B[U3, on = .(chr, pos >= lo, pos <= hi), nomatch = 0L, allow.cartesian = TRUE]; JP[, d := (bpos - tss) * sgn]; JP <- JP[d >= -3000 & d <= 1000]; JP[, k := exp(-abs(d) / 3000)]
for (col in grep("^(all|pp14|het3)_", names(B), value = TRUE)) S[[col]] <- vec(JP[, .(x = max(get(col) * k)), by = g])
S$cal_prom_max <- vec(JP[, .(x = max(ex * k)), by = g]); S$cal_prom_sum <- vec(JP[, .(x = sum(ex * k)), by = g]); rm(JP)
# (2) distal calibrated
Bx <- B[ex > 0, .(chr, bpos, ex)]; U1 <- copy(U)[, `:=`(lo = tss - 100000L, hi = tss + 100000L)]
JD <- Bx[U1, on = .(chr, bpos >= lo, bpos <= hi), nomatch = 0L, allow.cartesian = TRUE, .(g = i.g, d = (x.bpos - i.tss) * i.sgn, ex)]
cat("distal rows:", nrow(JD), "\n")
for (tau in c(3000, 10000, 20000, 50000)) { k <- exp(-abs(JD$d) / tau); tg <- paste0("cal_distal_tau", tau / 1000, "k")
  S[[paste0(tg, "_max")]] <- vec(JD[, .(x = max(ex * k)), by = g]); S[[paste0(tg, "_sum")]] <- vec(JD[, .(x = sum(ex * k)), by = g]) }
rm(JD, B)
saveRDS(S, file.path(out, "pcutoff_and_distal_scores.rds"))
panels <- build_gold_panels(gi, tss)
cat("panel sizes:\n"); print(lapply(panels, function(p) lengths(p)))
R <- eval_panels(S, panels, gi); fwrite(R, file.path(out, "pcutoff_and_distal_auroc_by_goldset.csv"))
PM <- panel_means(R); fwrite(PM, file.path(out, "pcutoff_and_distal_auroc_panel_means.csv"))
print(PM[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)], nrows = 60)

# paired bootstrap of the overall (mean over the five panels) difference for selected comparisons
boot_pair <- function(a, b) {
  per <- lapply(names(panels), function(pn) { bc <- NULL; for (s in 1:30) { set.seed(3000 + s); bc <- tryCatch(boot_compare(S[c(a, b)], panels[[pn]], gi, ND = 100L, B = 300L), error = function(e) NULL); if (!is.null(bc)) break }
    rowMeans(bc$ba[, a, , drop = FALSE] - bc$ba[, b, , drop = FALSE]) })
  d <- Reduce(`+`, per) / length(per); data.table(a = a, b = b, mean_diff = mean(d), lo = unname(quantile(d, .025)), hi = unname(quantile(d, .975)), per_panel = paste(round(sapply(per, mean), 3), collapse = " / ")) }
best_cut <- function(sub) { m <- PM[grepl(paste0("^", sub, "_"), score) & !grepl("distal|cal_", score)]; m$score[which.max(m$mean_all_panels)] }
bp <- c(best_cut("pp14"), best_cut("het3"), PM[grepl("^cal_distal", score)]$score[1])
cmp <- rbindlist(list(boot_pair(best_cut("pp14"), "pp14_1e-4"), boot_pair(best_cut("het3"), "het3_1e-4"), boot_pair(PM[grepl("^cal_distal", score)]$score[1], "pp14_1e-4"), boot_pair("cal_prom_max", "pp14_1e-4")))
cat("\nbootstrap (panels: PPARgene / heart / PPARG promoter / PPARD promoter / literature):\n"); print(cmp[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)]); fwrite(cmp, file.path(out, "pcutoff_and_distal_paired_bootstrap.csv"))
