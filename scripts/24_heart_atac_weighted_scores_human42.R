# Does restricting / weighting motif hits by human heart open chromatin help? ENCODE ENCSR117PYB (ATAC-seq, adult left-ventricle tissue, GRCh38), IDR thresholded peaks ENCFF174HPZ (26,918 peaks, 30 Mb).
# Scores (500 bp bins of the 42-motif FIMO hits; het = the 3 PPAR:RXR heterodimer motifs, all = all 42; value = relative FIMO score):
#   prom_het / prom_base            reference: -3 kb/+1 kb, exp(-|d|/3000), max
#   atac_prom_het / atac_prom_all   same, only bins inside heart ATAC peaks
#   atac_wide_het / atac_wide_all   +-100 kb, only ATAC bins, exp(-|d|/20 kb), max
#   rankavg_het                     mean rank of prom_het and atac_wide_het
# Evaluation as in 21/22: PPARgene + mouse-heart gold sets and genes with a promoter ChIP peak (PPARG/PPARA/RXRA, ChIP-Atlas), matched-background AUROC.
BENCH <- "C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Minor internship/R_code/benchmark"
source(file.path(BENCH, "config.R")); source(file.path(BENCH, "helpers_hits.R")); source(file.path(BENCH, "helpers_gold.R"))
suppressPackageStartupMessages(library(GenomicRanges))
set.seed(9)
out <- "C:/Users/brouw/PPRE_cistrome_atlas/results/human42_pipeline/distal_kernels"
gi <- load_genes(); tss <- readRDS(file.path(cache_dir, "tss_manuscript.rds"))
B <- readRDS("C:/Users/brouw/PPRE_cistrome_atlas/results/human42_pipeline_full_local/fimo_human42_bins500.rds"); B[, pos := bin * 500L + 250L]
atac <- fread(file.path(BENCH, "external_data/heart_chromatin/ENCFF174HPZ.bed.gz"), select = 1:3, col.names = c("chr", "start", "end"))
ga <- reduce(GRanges(atac$chr, IRanges(atac$start + 1L, atac$end)))
B[, atac := countOverlaps(GRanges(chr, IRanges(pos, pos)), ga) > 0]
cat("ATAC peaks (merged):", length(ga), " bins in ATAC:", sum(B$atac), "of", nrow(B), " (", round(100 * mean(B$atac), 2), "%)\n")
U <- unique(tss[in_universe == TRUE & ENSEMBL %in% gi$ENSEMBL, .(ENSEMBL, chr, strand, tss)]); U[, `:=`(g = match(ENSEMBL, gi$ENSEMBL), lo = tss - 100000L, hi = tss + 100000L, sgn = ifelse(strand == "+", 1L, -1L))]
J <- B[U, on = .(chr, pos >= lo, pos <= hi), nomatch = 0L, allow.cartesian = TRUE, .(g = i.g, d = (x.pos - i.tss) * i.sgn, max_rel, het_rel, atac)]
nG <- nrow(gi); sc <- function(D, col, k) { r <- rep(0, nG); a <- D[, .(x = max(get(col) * k)), by = g]; r[a$g] <- a$x; r }
P <- J[d >= -3000 & d <= 1000]; kP <- exp(-abs(P$d) / 3000); PA <- P[atac == TRUE]; kPA <- exp(-abs(PA$d) / 3000)
JA <- J[atac == TRUE]; kW <- exp(-abs(JA$d) / 20000)
S <- list(prom_base = sc(P, "max_rel", kP), prom_het = sc(P, "het_rel", kP), atac_prom_all = sc(PA, "max_rel", kPA), atac_prom_het = sc(PA, "het_rel", kPA),
          atac_wide_all = sc(JA, "max_rel", kW), atac_wide_het = sc(JA, "het_rel", kW))
S$rankavg_het <- (frank(S$prom_het) + frank(S$atac_wide_het)) / 2
saveRDS(S, file.path(out, "heart_atac_weighted_scores.rds"))

gold <- load_gold_all(gi); gold <- gold[c(PG_SETS, HEART_SETS)]
rd <- function(f) { d <- fread(file.path(cache_dir, f), header = FALSE); setnames(d, c("chr", "start", "end", "srx", "group", "cell", "score")); d[chr %in% STD_CHR & score >= 250][, group := gsub("%20", " ", group)] }
Pk <- list(PPARG = rd("chipatlas_PPARG_min.tsv"), PPARA = rd("chipatlas_PPARA_min.tsv"), RXRA = rd("chipatlas_RXRA_min.tsv"))
u <- tss[in_universe == TRUE]; w <- GRanges(u$chr, IRanges(ifelse(u$strand == "+", u$tss - 3000, u$tss - 1000), ifelse(u$strand == "+", u$tss + 1000, u$tss + 3000))); w <- trim(suppressWarnings(w)); w$ENSEMBL <- u$ENSEMBL
pk <- list()
for (tf in names(Pk)) for (g in unique(Pk[[tf]]$group)) { d <- Pk[[tf]][group == g]; if (nrow(d) < 300) next
  h <- unique(w$ENSEMBL[overlapsAny(w, GRanges(d$chr, IRanges(d$start + 1L, d$end)))]); if (length(h) >= 30) pk[[paste(tf, g)]] <- if (length(h) > 800) sample(h, 800) else h }
runset <- function(G, tag) { ev <- NULL; for (s in 1:20) { set.seed(100 + s); ev <- tryCatch(make_evaluator(gi, G, ndraw = 300L), error = function(e) NULL); if (!is.null(ev)) break }; stopifnot(!is.null(ev)); R <- rbindlist(lapply(names(S), function(nm) ev(S[[nm]], nm)$table)); R[, panel := tag]; R }
R <- rbind(runset(gold, "PPARgene + heart"), runset(pk, "promoter-peak genes")); fwrite(R, file.path(out, "heart_atac_weighted_auroc_by_goldset.csv"))
W <- dcast(R, panel + gold_set + n_gold ~ score, value.var = "auroc_matched"); cols <- names(S)
sm <- W[, lapply(.SD, mean), by = panel, .SDcols = cols]
pg <- W[panel == "promoter-peak genes" & grepl("^PPARG", gold_set), lapply(.SD, mean), .SDcols = cols][, panel := "PPARG promoter-peak (5 groups)"]
pa <- W[panel == "promoter-peak genes" & !grepl("^PPARG", gold_set), lapply(.SD, mean), .SDcols = cols][, panel := "PPARA+RXRA promoter-peak"]
sm <- rbind(sm, pg, pa, fill = TRUE); fwrite(sm, file.path(out, "heart_atac_weighted_auroc_panel_means.csv"))
print(sm[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])
bc <- NULL; for (s in 1:20) { set.seed(200 + s); bc <- tryCatch(boot_compare(S, gold, gi, ND = 100L, B = 400L), error = function(e) NULL); if (!is.null(bc)) break }
d <- function(a, b) { x <- rowMeans(bc$ba[, a, , drop = FALSE] - bc$ba[, b, , drop = FALSE]); data.table(panel = "PPARgene + heart", a = a, b = b, diff = mean(x), lo = unname(quantile(x, .025)), hi = unname(quantile(x, .975))) }
bp <- NULL; for (s in 1:20) { set.seed(300 + s); bp <- tryCatch(boot_compare(S, pk[grepl("^PPARG", names(pk))], gi, ND = 100L, B = 400L), error = function(e) NULL); if (!is.null(bp)) break }
d2 <- function(a, b) { x <- rowMeans(bp$ba[, a, , drop = FALSE] - bp$ba[, b, , drop = FALSE]); data.table(panel = "PPARG promoter-peak", a = a, b = b, diff = mean(x), lo = unname(quantile(x, .025)), hi = unname(quantile(x, .975))) }
cb <- rbind(d("atac_prom_het", "prom_het"), d("atac_wide_het", "prom_het"), d("rankavg_het", "prom_het"), d("atac_prom_all", "prom_base"), d("atac_wide_all", "prom_base"),
            d2("atac_prom_het", "prom_het"), d2("atac_wide_het", "prom_het"), d2("rankavg_het", "prom_het"), d2("atac_prom_all", "prom_base"), d2("atac_wide_all", "prom_base"))
print(cb[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)]); fwrite(cb, file.path(out, "heart_atac_weighted_paired_bootstrap.csv"))
