# Ceiling test (as 27) with the literature-curated PPRE gene set (literature_PPRE_human_prioritised.csv) as the gold standard.
# Gene scores: real ChIP-seq binding (strongest peak in the promoter window or +-10 kb) for ChIP-Atlas PPARA/PPARG/RXRA groups and PPARD myofibroblasts, versus the motif scores (heterodimer 3, PPAR-containing 14, all 42, ChIP-calibrated model).
# Gold sets: tier 1 (TSS-referenced, in -3 kb/+1 kb window), tier 1+2 (same), all tiers 1-7 (any position; a gene counts if it has a literature PPRE anywhere), and all tiers excluding rows flagged weak/negative.
BENCH <- "C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Minor internship/R_code/benchmark"
source(file.path(BENCH, "config.R")); source(file.path(BENCH, "helpers_hits.R")); source(file.path(BENCH, "helpers_gold.R"))
suppressPackageStartupMessages({ library(GenomicRanges); library(rtracklayer) })
source("C:/Users/brouw/PPRE_cistrome_atlas/scripts/gold_panels_human42.R")
out <- "C:/Users/brouw/PPRE_cistrome_atlas/results/human42_pipeline/ceiling_test"; dir.create(out, showWarnings = FALSE)
gi <- load_genes(); tss <- readRDS(file.path(cache_dir, "tss_manuscript.rds")); u <- tss[in_universe == TRUE]
mk <- function(up, dn) { w <- GRanges(u$chr, IRanges(ifelse(u$strand == "+", u$tss - up, u$tss - dn), ifelse(u$strand == "+", u$tss + dn, u$tss + up))); w <- trim(suppressWarnings(w)); w$ENSEMBL <- u$ENSEMBL; w }
Wn <- list(prom = mk(3000, 1000), `10kb` = mk(10000, 10000))
rd <- function(f) { d <- fread(file.path(cache_dir, f), header = FALSE); setnames(d, c("chr", "start", "end", "srx", "group", "cell", "score")); d[chr %in% STD_CHR & score >= 250][, group := gsub("%20", " ", group)] }
sets <- list()
for (tf in c("PPARA", "PPARG", "RXRA")) { d <- rd(paste0("chipatlas_", tf, "_min.tsv")); for (g in unique(d$group)) { x <- d[group == g]; if (nrow(x) >= 1000) sets[[paste(tf, g)]] <- GRanges(x$chr, IRanges(x$start + 1L, x$end), score = x$score) } }
pk <- fread(file.path(BENCH, "external_data/PPARD_myofibroblast/PPAR_DMSO_peaks.tsv"), check.names = FALSE)
lf <- liftOver(GRanges(pk$chr, IRanges(pk$start + 1L, pk$end), score = pk$`-log(p-value)`), import.chain(chain_file)); g38 <- unlist(lf[lengths(lf) == 1]); sets[["PPARD myofibroblast"]] <- g38[seqnames(g38) %in% STD_CHR]
chipscore <- function(gr, w) { ov <- findOverlaps(w, gr); a <- data.table(g = match(w$ENSEMBL[queryHits(ov)], gi$ENSEMBL), s = mcols(gr)$score[subjectHits(ov)])[!is.na(g), .(x = max(s)), by = g]; r <- rep(0, nrow(gi)); r[a$g] <- a$x; r }
S <- list(); for (n in names(sets)) for (wn in names(Wn)) S[[paste("ChIP:", n, "|", wn)]] <- chipscore(sets[[n]], Wn[[wn]])
Sub <- readRDS("C:/Users/brouw/PPRE_cistrome_atlas/results/human42_pipeline/distal_kernels/motif_subset_scores.rds"); C <- readRDS("C:/Users/brouw/PPRE_cistrome_atlas/results/human42_pipeline/chip_calibration/chip_calibrated_gene_scores.rds")
S$`MOTIF: heterodimer 3` <- Sub$heterodimer_3; S$`MOTIF: PPAR-containing 14` <- Sub$PPAR_containing_14; S$`MOTIF: all 42` <- Sub$all_42; S$`MOTIF: ChIP-calibrated model (top-3 tiles)` <- C$chip_model_top3_mean
# literature gold sets
lit <- fread("C:/Users/brouw/Downloads/literature_PPRE_human_prioritised.csv")[is_human == "yes"]
cand <- function(s) unique(toupper(trimws(unlist(strsplit(gsub("[()/,]", " ", s), "\\s+")))))
lit[, ENSEMBL := vapply(Gene_symbol, function(s) { m <- gi[toupper(gi$SYMBOL) %in% cand(s)]; if (nrow(m)) m$ENSEMBL[1] else NA_character_ }, character(1))]
cat("literature rows:", nrow(lit), " mapped:", sum(!is.na(lit$ENSEMBL)), "\n"); print(lit[, .(rows = .N, mapped_genes = uniqueN(ENSEMBL[!is.na(ENSEMBL)])), by = tier][order(tier)])
inwin <- lit[window_final == "yes" & ref_class == "TSS" & !is.na(ENSEMBL)]
G <- list(`tier 1 (in window)` = unique(inwin[tier == 1, ENSEMBL]), `tier 1+2 (in window)` = unique(inwin[tier %in% 1:2, ENSEMBL]),
          `all tiers 1-7, any position` = unique(lit[!is.na(ENSEMBL), ENSEMBL]), `all tiers, not flagged weak/negative` = unique(lit[!is.na(ENSEMBL) & flag_negative_or_weak != "yes", ENSEMBL]))
print(sapply(G, length))
ev <- robust_evaluator(gi, G, 1500)
R <- rbindlist(lapply(names(S), function(nm) ev(S[[nm]], nm)$table)); fwrite(R, file.path(out, "ceiling_test_literature_geneset_auroc.csv"))
W <- dcast(R, score ~ gold_set, value.var = "auroc_matched"); W[, mean_over_sets := rowMeans(as.matrix(.SD)), .SDcols = names(G)]; W <- W[order(-mean_over_sets)]
print(W[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)], nrows = 60)
fwrite(W, file.path(out, "ceiling_test_literature_geneset_auroc_summary.csv"))
