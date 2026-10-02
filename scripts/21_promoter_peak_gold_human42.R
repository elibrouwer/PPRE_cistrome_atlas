# Test the promoter-based gene scores against a gold standard that matches what they measure: genes with a ChIP-seq peak close to the TSS (promoter window -3 kb/+1 kb; also +-1 kb).
# Peak sets: ChIP-Atlas PPARG, PPARA and RXRA (MACS2 score >= 250), per cell group. Positives = genes with >= 1 peak in the window (<= 800 sampled), negatives = promoter-length/GC matched genes (5 per positive).
# Scores: near-TSS FIMO scores (all-motif, heterodimer-only), the adipocyte-PPARG-calibrated model, plus baselines (promoter GC, cardiomyocyte expression).
# NOTE: the calibrated model was trained on adipocyte PPARG peaks, so the PPARG adipocyte rows are in-sample for that score; all other groups are transfer tests.
BENCH <- "C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Minor internship/R_code/benchmark"
source(file.path(BENCH, "config.R")); source(file.path(BENCH, "helpers_hits.R")); source(file.path(BENCH, "helpers_gold.R"))
suppressPackageStartupMessages(library(GenomicRanges))
set.seed(3)
out <- "C:/Users/brouw/PPRE_cistrome_atlas/results/human42_pipeline/chip_calibration"
gi <- load_genes(); tss <- readRDS(file.path(cache_dir, "tss_manuscript.rds"))
K <- readRDS("C:/Users/brouw/PPRE_cistrome_atlas/results/human42_pipeline/distal_kernels/distal_kernel_scores.rds")
C <- readRDS(file.path(out, "chip_calibrated_gene_scores.rds"))
J <- readRDS(file.path(res_dir, "J_candidate_gene_table.rds"))[match(gi$ENSEMBL, ENSEMBL)]
ex <- ifelse(is.na(J$log10_expression), log10(1e-4), J$log10_expression)
S <- list(prom_base = K$prom_base, prom_het = K$prom_het, chip_model_max = C$chip_model_max, chip_model_top3_mean = C$chip_model_top3_mean, `promoter GC` = gi$gc, `cardiomyocyte expression` = ex)
rd <- function(f) { d <- fread(file.path(cache_dir, f), header = FALSE); setnames(d, c("chr", "start", "end", "srx", "group", "cell", "score")); d[chr %in% STD_CHR & score >= 250][, group := gsub("%20", " ", group)] }
P <- list(PPARG = rd("chipatlas_PPARG_min.tsv"), PPARA = rd("chipatlas_PPARA_min.tsv"), RXRA = rd("chipatlas_RXRA_min.tsv"))
u <- tss[in_universe == TRUE]
mkwin <- function(up, dn) { w <- GRanges(u$chr, IRanges(ifelse(u$strand == "+", u$tss - up, u$tss - dn), ifelse(u$strand == "+", u$tss + dn, u$tss + up))); w <- trim(suppressWarnings(w)); w$ENSEMBL <- u$ENSEMBL; w }
W <- list(`promoter -3k/+1k` = mkwin(3000, 1000), `+-1kb` = mkwin(1000, 1000))
lab <- list()
for (tf in names(P)) for (g in unique(P[[tf]]$group)) { d <- P[[tf]][group == g]; if (nrow(d) < 300) next
  gr <- GRanges(d$chr, IRanges(d$start + 1L, d$end))
  for (wn in names(W)) { h <- unique(W[[wn]]$ENSEMBL[overlapsAny(W[[wn]], gr)]); if (length(h) >= 30) lab[[paste(tf, g, wn, sep = " | ")]] <- h } }
comp <- data.table(label = names(lab), n_genes = lengths(lab)); print(comp); fwrite(comp, file.path(out, "promoter_peak_gold_labels.csv"))
lab <- lapply(lab, function(g) if (length(g) > 800) sample(g, 800) else g)
ev <- make_evaluator(gi, lab, ndraw = 300L)
R <- rbindlist(lapply(names(S), function(nm) ev(S[[nm]], nm)$table)); fwrite(R, file.path(out, "promoter_peak_gold_auroc.csv"))
Wd <- dcast(R, gold_set + n_gold ~ score, value.var = "auroc_matched")
cols <- c("prom_base", "prom_het", "chip_model_max", "chip_model_top3_mean", "promoter GC", "cardiomyocyte expression")
print(Wd[, c("gold_set", "n_gold", cols), with = FALSE][, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])
