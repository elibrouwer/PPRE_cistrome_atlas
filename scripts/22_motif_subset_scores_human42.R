# Which motif subset should the promoter score use? Defined by biology before testing:
#   all     all 42 motifs
#   pp      PPAR-containing motifs (11 PPAR half-site/single + 3 PPAR:RXR heterodimer = 14)
#   het     the 3 PPAR:RXR heterodimer (DR1) motifs
#   ppar    the 11 PPAR-only motifs
#   rxr     the 28 RXR-only motifs (RXR partners with many receptors, so these are not PPAR-specific)
# Score (same as prom_base/prom_het in 16b): per gene, max over hits within -3 kb/+1 kb of a TSS of (FIMO score / best score of its motif) x exp(-|d|/3000), 500 bp bins.
# Evaluated on (a) PPARgene + mouse-heart gold sets, (b) genes with a promoter PPARG/PPARA/RXRA ChIP peak (ChIP-Atlas, per cell group), matched-background AUROC.
BENCH <- "C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Minor internship/R_code/benchmark"
source(file.path(BENCH, "config.R")); source(file.path(BENCH, "helpers_hits.R")); source(file.path(BENCH, "helpers_gold.R"))
suppressPackageStartupMessages(library(GenomicRanges))
set.seed(5)
out <- "C:/Users/brouw/PPRE_cistrome_atlas/results/human42_pipeline/distal_kernels"
fimo_root <- "C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Minor internship/PPAR_RXR_motifs_Hocomoco_Jaspar/FIMO_results_human42_unmasked"
inv <- fread("C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Minor internship/PPAR_RXR_motifs_Hocomoco_Jaspar/Motifs_update_2026-09-30/motif_inventory_gates.csv", encoding = "UTF-8")
cat_of <- setNames(inv$category, inv$motif_id)
gi <- load_genes(); tss <- readRDS(file.path(cache_dir, "tss_manuscript.rds"))
U <- unique(tss[in_universe == TRUE & ENSEMBL %in% gi$ENSEMBL, .(ENSEMBL, chr, strand, tss)]); U[, `:=`(g = match(ENSEMBL, gi$ENSEMBL), sgn = ifelse(strand == "+", 1L, -1L), lo = tss - 3000L, hi = tss + 3000L)]
acc <- list()
for (f in list.files(fimo_root, "fimo\\.tsv$", recursive = TRUE, full.names = TRUE)) {
  d <- fread(f, sep = "\t", select = c("motif_id", "sequence_name", "start", "stop", "score"), colClasses = list(character = c("motif_id", "sequence_name")))
  d <- d[sequence_name %in% as.character(c(1:22, "X")) & !startsWith(motif_id, "#")]
  d[, `:=`(start = as.integer(start), stop = as.integer(stop), score = as.numeric(score))]; d[, rel := score / max(score), by = motif_id]
  cg <- cat_of[d$motif_id]
  d[, `:=`(pos = (start + stop) %/% 2L, r_all = rel, r_pp = fifelse(cg %in% c("PPAR", "PPAR:RXR heterodimer"), rel, 0), r_het = fifelse(cg == "PPAR:RXR heterodimer", rel, 0),
           r_ppar = fifelse(cg == "PPAR", rel, 0), r_rxr = fifelse(cg == "RXR", rel, 0), chr = paste0("chr", sequence_name))]
  d[, bin := pos %/% 500L]
  acc[[f]] <- d[, .(r_all = max(r_all), r_pp = max(r_pp), r_het = max(r_het), r_ppar = max(r_ppar), r_rxr = max(r_rxr)), by = .(chr, bin)]; cat(basename(dirname(f)), "done\n")
}
B <- rbindlist(acc)[, lapply(.SD, max), by = .(chr, bin), .SDcols = c("r_all", "r_pp", "r_het", "r_ppar", "r_rxr")]; B[, pos := bin * 500L + 250L]
J <- B[U, on = .(chr, pos >= lo, pos <= hi), nomatch = 0L, allow.cartesian = TRUE, .(g = i.g, d = (x.pos - i.tss) * i.sgn, r_all, r_pp, r_het, r_ppar, r_rxr)]
J <- J[d >= -3000 & d <= 1000]; k <- exp(-abs(J$d) / 3000); nG <- nrow(gi)
gagg <- function(D, val, fun = "max") { D[, vv := val]; a <- if (fun == "max") D[, .(x = max(vv)), by = g] else D[, .(x = sum(vv)), by = g]; D[, vv := NULL]; r <- rep(0, nG); r[a$g] <- a$x; r }   # val is stored as a column so it is grouped with the rows (a global vector inside j would be recycled over the whole table)
sc <- function(col) gagg(J, J[[col]] * k, "max")
S <- list(all_42 = sc("r_all"), PPAR_containing_14 = sc("r_pp"), heterodimer_3 = sc("r_het"), PPAR_only_11 = sc("r_ppar"), RXR_only_28 = sc("r_rxr"))
saveRDS(S, file.path(out, "motif_subset_scores.rds"))

# gold standards
gold <- load_gold_all(gi); gold <- gold[c(PG_SETS, HEART_SETS)]
rd <- function(f) { d <- fread(file.path(cache_dir, f), header = FALSE); setnames(d, c("chr", "start", "end", "srx", "group", "cell", "score")); d[chr %in% STD_CHR & score >= 250][, group := gsub("%20", " ", group)] }
P <- list(PPARG = rd("chipatlas_PPARG_min.tsv"), PPARA = rd("chipatlas_PPARA_min.tsv"), RXRA = rd("chipatlas_RXRA_min.tsv"))
u <- tss[in_universe == TRUE]; w <- GRanges(u$chr, IRanges(ifelse(u$strand == "+", u$tss - 3000, u$tss - 1000), ifelse(u$strand == "+", u$tss + 1000, u$tss + 3000))); w <- trim(suppressWarnings(w)); w$ENSEMBL <- u$ENSEMBL
pk <- list()
for (tf in names(P)) for (g in unique(P[[tf]]$group)) { d <- P[[tf]][group == g]; if (nrow(d) < 300) next
  h <- unique(w$ENSEMBL[overlapsAny(w, GRanges(d$chr, IRanges(d$start + 1L, d$end)))]); if (length(h) >= 30) pk[[paste(tf, g)]] <- if (length(h) > 800) sample(h, 800) else h }
runset <- function(G, tag) {
  ev <- make_evaluator(gi, G, ndraw = 300L); R <- rbindlist(lapply(names(S), function(nm) ev(S[[nm]], nm)$table)); R[, panel := tag]; R }
R1 <- runset(gold, "PPARgene + heart"); R2 <- runset(pk, "promoter-peak genes")
R <- rbind(R1, R2); fwrite(R, file.path(out, "motif_subset_auroc_by_goldset.csv"))
W <- dcast(R, panel + gold_set + n_gold ~ score, value.var = "auroc_matched")
cols <- names(S); W <- W[, c("panel", "gold_set", "n_gold", cols), with = FALSE]
print(W[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)], nrows = 60)
sm <- W[, lapply(.SD, mean), by = panel, .SDcols = cols]; print(sm[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])
fwrite(W, file.path(out, "motif_subset_auroc_wide.csv")); fwrite(sm, file.path(out, "motif_subset_auroc_panel_means.csv"))
# paired bootstrap: PPAR-containing (14) vs all (42) and vs heterodimer (3); PPARG promoter-peak sets and PPARgene+heart sets
bcmp <- function(Gset, tag) { bc <- boot_compare(S, Gset, gi, ND = 100L, B = 400L); sets <- names(Gset)
  f <- function(a, b) { d <- rowMeans(bc$ba[, a, sets, drop = FALSE] - bc$ba[, b, sets, drop = FALSE]); data.table(panel = tag, a = a, b = b, diff = mean(d), lo = unname(quantile(d, .025)), hi = unname(quantile(d, .975))) }
  rbind(f("PPAR_containing_14", "all_42"), f("heterodimer_3", "all_42"), f("PPAR_containing_14", "heterodimer_3"), f("PPAR_only_11", "all_42")) }
cb <- rbind(bcmp(gold, "PPARgene + heart"), bcmp(pk[grepl("^PPARG", names(pk))], "PPARG promoter-peak genes"))
print(cb[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)]); fwrite(cb, file.path(out, "motif_subset_paired_bootstrap.csv"))
