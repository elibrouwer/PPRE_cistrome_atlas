# Helper (sourced by 28b and later scripts): builds the gold-standard panels used throughout the 42-motif score validation.
#   PPARgene       4 PPARgene sets (human-derived, mouse-only, PPARA, PPARG)            (benchmark load_gold_all)
#   heart          4 mouse-heart PPARA-response sets (GEO GSE30495, GSE33101, GSE30553)  (benchmark)
#   PPARG_promoter genes with a promoter (-3 kb/+1 kb) PPARG ChIP-Atlas peak, per cell group (score >= 250, >= 300 peaks, >= 30 genes; <= 800 sampled)
#   PPARD_promoter genes with a promoter PPARD peak, human myofibroblasts (E-MTAB-371, hg19 -> hg38)
#   literature     literature-curated human PPREs, tier 1+2, TSS-referenced, within the promoter window (literature_PPRE_human_prioritised.csv)
# Requires config.R / helpers_hits.R / helpers_gold.R of ./benchmark to be sourced and gi, tss to exist.
build_gold_panels <- function(gi, tss, seed = 1) {
  set.seed(seed)
  g <- load_gold_all(gi)
  panels <- list(PPARgene = g[PG_SETS], heart = g[HEART_SETS])
  u <- tss[in_universe == TRUE]
  w <- GRanges(u$chr, IRanges(ifelse(u$strand == "+", u$tss - 3000, u$tss - 1000), ifelse(u$strand == "+", u$tss + 1000, u$tss + 3000))); w <- trim(suppressWarnings(w)); w$ENSEMBL <- u$ENSEMBL
  d <- fread(file.path(cache_dir, "chipatlas_PPARG_min.tsv"), header = FALSE); setnames(d, c("chr", "start", "end", "srx", "group", "cell", "score"))
  d <- d[chr %in% STD_CHR & score >= 250][, group := gsub("%20", " ", group)]
  pg <- list(); for (gr_ in unique(d$group)) { x <- d[group == gr_]; if (nrow(x) < 300) next
    h <- unique(w$ENSEMBL[overlapsAny(w, GRanges(x$chr, IRanges(x$start + 1L, x$end)))]); h <- intersect(h, gi$ENSEMBL); if (length(h) >= 30) pg[[paste("PPARG", gr_)]] <- if (length(h) > 800) sample(h, 800) else h }
  panels$PPARG_promoter <- pg
  pk <- fread(file.path(BENCH, "external_data/PPARD_myofibroblast/PPAR_DMSO_peaks.tsv"), check.names = FALSE)
  lf <- liftOver(GRanges(pk$chr, IRanges(pk$start + 1L, pk$end)), import.chain(chain_file)); g38 <- unlist(lf[lengths(lf) == 1]); g38 <- g38[seqnames(g38) %in% STD_CHR]
  h <- intersect(unique(w$ENSEMBL[overlapsAny(w, g38)]), gi$ENSEMBL); panels$PPARD_promoter <- list(`PPARD myofibroblast` = if (length(h) > 800) sample(h, 800) else h)
  lit <- fread("C:/Users/brouw/Downloads/literature_PPRE_human_prioritised.csv")[tier %in% 1:2 & is_human == "yes" & window_final == "yes" & ref_class == "TSS"]
  cand <- function(s) unique(toupper(trimws(unlist(strsplit(gsub("[()/,]", " ", s), "\\s+")))))
  e <- unique(unlist(lapply(lit$Gene_symbol, function(s) { m <- gi[toupper(gi$SYMBOL) %in% cand(s)]; if (nrow(m)) m$ENSEMBL[1] else NULL })))
  panels$literature <- list(`literature tier 1+2` = e)
  panels
}
# retry-with-new-seed wrapper around the benchmark evaluator (the matcher can run out of candidate genes for large gold sets)
robust_evaluator <- function(gi, G, seed0 = 100) { for (s in 1:30) { set.seed(seed0 + s); ev <- tryCatch(make_evaluator(gi, G, ndraw = 300L), error = function(e) NULL); if (!is.null(ev)) return(ev) }; stop("evaluator failed") }
eval_panels <- function(S, panels, gi) {
  rbindlist(lapply(names(panels), function(pn) { ev <- robust_evaluator(gi, panels[[pn]]); R <- rbindlist(lapply(names(S), function(nm) ev(S[[nm]], nm)$table)); R[, panel := pn]; R }))
}
panel_means <- function(R) { m <- R[, .(auroc = mean(auroc_matched)), by = .(panel, score)]; w <- dcast(m, score ~ panel, value.var = "auroc"); cols <- setdiff(names(w), "score"); w[, mean_all_panels := rowMeans(as.matrix(.SD)), .SDcols = cols]; w[order(-mean_all_panels)] }
