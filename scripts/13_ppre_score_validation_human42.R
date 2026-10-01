# Robustness and validation of the PPRE score for the 42 human-source motifs.
# Score variants: (a) original per-hit score * decay, aggregated over raw hits (sum / max / median);
#                 (b) locus-merged: hits of one gene within the same 100 bp TSS-relative window (any motif, any strand) are collapsed to one locus (best hit), then max / sum over loci;
#                 (c) normalised: -log10(FIMO p) instead of the raw FIMO score, which is not comparable across motifs.
# Decay constants 1000 / 2000 / 3000 / 5000 bp and no decay.
# Validation: AUC of each variant for separating genes of independent, motif-free GO gene sets (fatty acid beta-oxidation, fatty acid catabolism,
# lipid catabolism) from the other expressed genes, with a bootstrap CI and a label-permutation null.
suppressMessages({library(data.table); library(ggplot2); library(org.Hs.eg.db); library(AnnotationDbi)})
root <- "C:/Users/brouw/PPRE_cistrome_atlas/results"
ck <- file.path(root, "human42_pipeline_full_local/ChIPseeker_output")
out <- file.path(root, "human42_pipeline/ppre_score_validation"); dir.create(out, showWarnings = FALSE)
universe <- unique(fread(file.path(root, "human42_pipeline/fimo_analysis/Bedfiles_output/PPRE_gtex_ranked_scores_abs_3k_test.bed"))$ensembl)
cat("expressed genes with promoter hit:", length(universe), "\n")

hits <- rbindlist(lapply(list.files(ck, "\\.txt$", full.names = TRUE), function(f)
  fread(f, sep = " ", select = c("seqnames", "start", "end", "score", "pvalue", "motif_id", "distanceToTSS", "ENSEMBL", "SYMBOL"))))
hits <- hits[ENSEMBL %in% universe]
hits[, `:=`(nlp = -log10(pmax(pvalue, 1e-300)), absd = abs(distanceToTSS))]
cat("promoter hits in expressed genes:", nrow(hits), "\n")

# locus definition: fixed 100 bp windows relative to the TSS within a gene (chaining overlapping hits collapses the whole promoter into one locus,
# because 42 motifs at p <= 1e-4 tile it continuously)
hits[, locus := floor(distanceToTSS / 100)]
cat("loci:", uniqueN(hits[, .(ENSEMBL, locus)]), "of", nrow(hits), "hits\n")

genes <- data.table(ENSEMBL = universe)
variants <- list()
for (dc in c(1000, 2000, 3000, 5000, Inf)) {
  w <- if (is.infinite(dc)) 1 else exp(-hits$absd / dc)
  hits[, `:=`(s_raw = score * w, s_nlp = nlp * w)]
  lab <- if (is.infinite(dc)) "nodecay" else paste0("d", dc)
  for (sc in c("raw", "nlp")) {
    col <- paste0("s_", sc)
    r <- hits[, .(sum = sum(get(col)), max = max(get(col)), median = median(get(col))), by = ENSEMBL]
    l <- hits[, .(v = max(get(col))), by = .(ENSEMBL, locus)][, .(locus_max = max(v), locus_sum = sum(v)), by = ENSEMBL]
    m <- merge(r, l, by = "ENSEMBL")
    for (k in c("sum", "max", "median", "locus_max", "locus_sum"))
      variants[[paste(sc, lab, k, sep = "|")]] <- m[[k]][match(universe, m$ENSEMBL)]
  }
}
V <- as.data.table(variants); V[, ENSEMBL := universe]
fwrite(V, file.path(out, "ppre_score_variants_per_gene.csv"))

# independent gene sets
go2g <- function(id) {
  e <- unique(unlist(as.list(org.Hs.egGO2ALLEGS[id]))); ens <- unique(unlist(as.list(org.Hs.egENSEMBL[intersect(e, mappedkeys(org.Hs.egENSEMBL))]))); intersect(ens, universe)
}
sets <- list("fatty acid beta-oxidation (GO:0006635)" = go2g("GO:0006635"), "fatty acid catabolic process (GO:0009062)" = go2g("GO:0009062"),
             "lipid catabolic process (GO:0016042)" = go2g("GO:0016042"))
print(sapply(sets, length))
auc <- function(x, y) { r <- rank(x, na.last = "keep"); ok <- !is.na(x); r <- rank(x[ok]); yy <- y[ok]; n1 <- sum(yy); n0 <- sum(!yy); (sum(r[yy]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
set.seed(1)
res <- rbindlist(lapply(names(sets), function(sn) {
  y <- universe %in% sets[[sn]]
  rbindlist(lapply(names(variants), function(vn) {
    x <- variants[[vn]]; a <- auc(x, y)
    bs <- replicate(200, { i <- sample(seq_along(y), replace = TRUE); auc(x[i], y[i]) })
    nul <- replicate(200, auc(x, sample(y)))
    data.table(set = sn, n_genes = sum(y), variant = vn, AUC = a, lo = quantile(bs, .025, na.rm = TRUE), hi = quantile(bs, .975, na.rm = TRUE),
               null_sd = sd(nul), z = (a - 0.5) / sd(nul))
  }))
}))
res[, c("score", "decay", "agg") := tstrsplit(variant, "|", fixed = TRUE)]
fwrite(res, file.path(out, "ppre_score_validation_auc.csv"))

# rank stability across decay constants (reference = original definition: raw, d3000, median) and top-5000 overlap
ref <- variants[["raw|d3000|median"]]
stab <- rbindlist(lapply(names(variants), function(vn) data.table(variant = vn,
  spearman_vs_original = cor(ref, variants[[vn]], method = "spearman", use = "complete.obs"),
  top5000_overlap = length(intersect(universe[order(-ref)][1:5000], universe[order(-variants[[vn]])][1:5000])) / 5000)))
fwrite(stab, file.path(out, "ppre_score_rank_stability.csv"))
cat("\nAUC, fatty acid beta-oxidation (best 12):\n"); print(res[set == names(sets)[1]][order(-AUC)][1:12, .(variant, AUC, lo, hi, z)])
cat("\noriginal (raw|d3000|median):\n"); print(res[variant == "raw|d3000|median", .(set, AUC, lo, hi, z)])

p <- ggplot(res[decay %in% c("d1000", "d2000", "d3000", "d5000", "nodecay")], aes(factor(decay, c("d1000", "d2000", "d3000", "d5000", "nodecay")), AUC, colour = paste(score, agg))) +
  geom_hline(yintercept = 0.5, linetype = 2) + geom_point(position = position_dodge(.6)) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, position = position_dodge(.6)) +
  facet_wrap(~set, ncol = 1) + labs(x = "decay constant", y = "AUC vs other expressed genes", colour = "FIMO value / aggregation") +
  theme_bw(base_size = 9) + theme(panel.grid.minor = element_blank())
ggsave(file.path(out, "ppre_score_validation_auc.png"), p, width = 8, height = 9, dpi = 160)
cat("done\n")
