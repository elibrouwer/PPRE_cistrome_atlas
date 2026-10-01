# Alternative gene-level scores built from the 42-motif FIMO promoter hits, compared by AUC against independent GO gene sets.
# Candidate scores (all with an exp(-|d|/3000) distance weight w unless stated):
#   base_raw_max        original max of FIMO score * w
#   rel_max             FIMO score / motif max observed score  (relative score, comparable across motifs)
#   pct_max             within-motif percentile of the score
#   z_max               within-motif z-score
#   aff_sum             TRAP-like affinity sum: sum 2^(score - motif max) * w
#   cl_sum / cl_noisyor best hit per r>=0.98 motif cluster (16 clusters), then summed / noisy-OR (removes the redundancy of 42 correlated motifs)
#   cl_strong_n         number of motif clusters with a hit in the top 10% of its motif within 1 kb of the TSS
#   het_max             rel score of the PPAR:RXR heterodimer motifs only
#   ppar_rxr_gm         geometric mean of the best PPAR-family and best RXR-family hit (both partners needed for the heterodimer)
#   top3_loci           mean of the three best 100 bp loci
#   gauss_max / flat1k_max   Gaussian (sigma 1500) and flat (+-1 kb) distance kernels
#   gcadj_rel_max       rel_max ranked within promoter-GC deciles
#   resid_rel_max       rel_max residual after regressing on promoter GC and number of hits
#   ens_rank            mean rank of rel_max, cl_sum and ppar_rxr_gm
# Evaluation: crude AUC, GC x length stratified AUC, and an honest "best-of" estimate (choose the score on a random half of the genes, evaluate on the other half).
suppressMessages({library(data.table); library(org.Hs.eg.db); library(AnnotationDbi); library(ggplot2)})
root <- "C:/Users/brouw/PPRE_cistrome_atlas/results"
ck <- file.path(root, "human42_pipeline_full_local/ChIPseeker_output")
outd <- file.path(root, "human42_pipeline/ppre_score_validation")
inv <- fread("C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Minor internship/PPAR_RXR_motifs_Hocomoco_Jaspar/Motifs_update_2026-09-30/motif_inventory_gates.csv", encoding = "UTF-8")
inv <- inv[, .(motif_id, category, cluster = cluster_r098)]
cov <- fread(file.path(outd, "gene_covariates_gc_length_expression.csv"))
cov[, loglen := log10(gene_len)]
hits <- rbindlist(lapply(list.files(ck, "\\.txt$", full.names = TRUE), function(f)
  fread(f, sep = " ", select = c("seqnames", "start", "end", "score", "pvalue", "motif_id", "distanceToTSS", "ENSEMBL"))))
hits <- hits[ENSEMBL %in% cov$ENSEMBL]
hits <- merge(hits, inv, by = "motif_id")
cat("hits:", nrow(hits), " motifs:", uniqueN(hits$motif_id), " clusters:", uniqueN(hits$cluster), "\n")
hits[, absd := abs(distanceToTSS)]
hits[, `:=`(mmax = max(score), mmu = mean(score), msd = sd(score)), by = motif_id]
hits[, `:=`(rel = score / mmax, z = (score - mmu) / msd, pct = frank(score) / .N, aff = 2^(score - mmax)), by = motif_id]
hits[, w := exp(-absd / 3000)]
genes <- cov$ENSEMBL
S <- list()
put <- function(nm, dt, col) S[[nm]] <<- dt[[col]][match(genes, dt$ENSEMBL)]
agg <- function(expr) hits[, .(v = eval(expr)), by = ENSEMBL]

put("base_raw_max", agg(quote(max(score * w))), "v")
put("rel_max", agg(quote(max(rel * w))), "v")
put("pct_max", agg(quote(max(pct * w))), "v")
put("z_max", agg(quote(max(z * w))), "v")
put("aff_sum", agg(quote(sum(aff * w))), "v")
cm <- hits[, .(r = max(rel * w), p = max(pct * w)), by = .(ENSEMBL, cluster)]
put("cl_sum", cm[, .(v = sum(r)), by = ENSEMBL], "v")
put("cl_noisyor", cm[, .(v = 1 - prod(1 - p)), by = ENSEMBL], "v")
put("cl_strong_n", hits[absd <= 1000 & pct >= 0.9, .(v = uniqueN(cluster)), by = ENSEMBL], "v"); S$cl_strong_n[is.na(S$cl_strong_n)] <- 0
put("het_max", hits[category == "PPAR:RXR heterodimer", .(v = max(rel * w)), by = ENSEMBL], "v"); S$het_max[is.na(S$het_max)] <- 0
pp <- hits[category %in% c("PPAR", "PPAR:RXR heterodimer"), .(p = max(rel * w)), by = ENSEMBL]
rx <- hits[category %in% c("RXR", "PPAR:RXR heterodimer"), .(r = max(rel * w)), by = ENSEMBL]
pr <- merge(pp, rx, all = TRUE, by = "ENSEMBL"); pr[is.na(p), p := 0]; pr[is.na(r), r := 0]; pr[, v := sqrt(p * r)]
put("ppar_rxr_gm", pr, "v"); S$ppar_rxr_gm[is.na(S$ppar_rxr_gm)] <- 0
hits[, locus := floor(distanceToTSS / 100)]
lo <- hits[, .(v = max(rel * w)), by = .(ENSEMBL, locus)][order(ENSEMBL, -v)][, .(v = mean(head(v, 3))), by = ENSEMBL]
put("top3_loci", lo, "v")
put("gauss_max", agg(quote(max(rel * exp(-(absd / 1500)^2)))), "v")
put("flat1k_max", agg(quote(max(rel * (absd <= 1000)))), "v")
nh <- hits[, .N, by = ENSEMBL]; nhg <- nh$N[match(genes, nh$ENSEMBL)]
dd <- data.table(rel = S$rel_max, gc = cov$promoter_gc, lnh = log(nhg))
S$resid_rel_max <- residuals(lm(rel ~ gc + I(gc^2) + lnh, data = dd))
gcdec <- cut(cov$promoter_gc, quantile(cov$promoter_gc, 0:10 / 10), include.lowest = TRUE, labels = FALSE)
S$gcadj_rel_max <- ave(S$rel_max, gcdec, FUN = function(x) frank(x) / length(x))
S$ens_rank <- (frank(S$rel_max) + frank(S$cl_sum) + frank(S$ppar_rxr_gm)) / 3
S <- S[c("base_raw_max", "rel_max", "pct_max", "z_max", "aff_sum", "cl_sum", "cl_noisyor", "cl_strong_n", "het_max", "ppar_rxr_gm", "top3_loci", "gauss_max", "flat1k_max", "gcadj_rel_max", "resid_rel_max", "ens_rank")]
fwrite(cbind(ENSEMBL = genes, as.data.table(S)), file.path(outd, "ppre_alternative_scores_per_gene.csv"))

go2g <- function(id) { e <- unique(unlist(as.list(org.Hs.egGO2ALLEGS[id]))); ens <- unique(unlist(as.list(org.Hs.egENSEMBL[intersect(e, mappedkeys(org.Hs.egENSEMBL))]))); intersect(ens, genes) }
sets <- list("FA beta-oxidation" = go2g("GO:0006635"), "FA catabolic" = go2g("GO:0009062"), "lipid catabolic" = go2g("GO:0016042"),
             "TCA cycle" = go2g("GO:0006099"), "OXPHOS" = go2g("GO:0006119"))
print(sapply(sets, length))
Y <- lapply(sets, function(s) genes %in% s)
strat <- paste(cut(cov$promoter_gc, quantile(cov$promoter_gc, 0:5 / 5), include.lowest = TRUE, labels = FALSE), cut(cov$loglen, quantile(cov$loglen, 0:3 / 3), include.lowest = TRUE, labels = FALSE))
auc <- function(x, y) { r <- rank(x); n1 <- sum(y); n0 <- sum(!y); (sum(r[y]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
sauc <- function(x, y, s) { t <- data.table(x, y, s)[, .(a = if (sum(y) > 0 && sum(!y) > 0) auc(x, y) else NA_real_, w = sum(y) * sum(!y)), by = s][!is.na(a)]; sum(t$a * t$w) / sum(t$w) }
res <- rbindlist(lapply(names(S), function(sn) rbindlist(lapply(names(Y), function(gs)
  data.table(score = sn, set = gs, AUC = auc(S[[sn]], Y[[gs]]), AUC_strat = sauc(S[[sn]], Y[[gs]], strat))))))
fwrite(res, file.path(outd, "ppre_alternative_scores_auc.csv"))
w <- dcast(res, score ~ set, value.var = "AUC"); w[, mean_AUC := rowMeans(.SD), .SDcols = names(sets)]
ws <- res[, .(mean_AUC_strat = mean(AUC_strat)), by = score]; w <- merge(w, ws, by = "score")[order(-mean_AUC)]
print(w[, lapply(.SD, function(v) if (is.numeric(v)) round(v, 3) else v)])

# honest best-of: choose on half A (by mean AUC over sets), evaluate on half B
set.seed(7); n <- length(genes); pick <- character(); evalB <- numeric(); baseB <- numeric()
for (i in 1:100) {
  a <- sample(n, n / 2); b <- setdiff(seq_len(n), a)
  mA <- sapply(S, function(x) mean(sapply(Y, function(y) auc(x[a], y[a]))))
  best <- names(which.max(mA)); pick <- c(pick, best)
  evalB <- c(evalB, mean(sapply(Y, function(y) auc(S[[best]][b], y[b]))))
  baseB <- c(baseB, mean(sapply(Y, function(y) auc(S$base_raw_max[b], y[b]))))
}
cat("\nbest-of selected on half A, evaluated on half B (100 splits): mean AUC", round(mean(evalB), 4), " baseline", round(mean(baseB), 4),
    " diff", round(mean(evalB - baseB), 4), " P(diff>0)", mean(evalB > baseB), "\n"); print(sort(table(pick), decreasing = TRUE))
# bootstrap for the top 3 vs baseline (mean AUC over sets)
top <- w$score[1:3]; cat("\npaired bootstrap of mean AUC difference vs baseline (200 resamples):\n")
for (sn in setdiff(top, "base_raw_max")) { d <- replicate(200, { i <- sample(n, replace = TRUE); mean(sapply(Y, function(y) auc(S[[sn]][i], y[i]))) - mean(sapply(Y, function(y) auc(S$base_raw_max[i], y[i]))) })
  cat(sn, ": diff", round(mean(d), 4), " 95% CI", round(quantile(d, c(.025, .975)), 4), "\n") }
p <- ggplot(res, aes(reorder(score, AUC, mean), AUC, colour = set)) + geom_hline(yintercept = 0.5, linetype = 2) + geom_point() + coord_flip() +
  labs(x = NULL, y = "AUC vs other expressed genes") + theme_bw(base_size = 9)
ggsave(file.path(outd, "ppre_alternative_scores_auc.png"), p, width = 7, height = 6, dpi = 160)
