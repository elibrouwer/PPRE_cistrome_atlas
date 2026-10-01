# Does the PPRE score still separate fatty-acid-oxidation genes after controlling for promoter GC, gene length and expression?
# Covariates: promoter GC (mean of the 1 kb windows covering TSS-3000..TSS+1000, from gc_bins_1kb.tsv, N-aware), log10 gene length, GTEx expression.
# (1) logistic regression: gene-set membership ~ scaled score + covariates; (2) stratified AUC (strata = GC quintile x length tertile) with a within-stratum permutation null.
suppressMessages({library(data.table); library(org.Hs.eg.db); library(AnnotationDbi)})
root <- "C:/Users/brouw/PPRE_cistrome_atlas/results"
V <- fread(file.path(root, "human42_pipeline/ppre_score_validation/ppre_score_variants_per_gene.csv"))
bed <- fread(file.path(root, "human42_pipeline/fimo_analysis/Bedfiles_output/PPRE_gtex_ranked_scores_abs_3k_test.bed"))
cov <- unique(bed[, .(ENSEMBL = ensembl, expr = GTEx_expression)], by = "ENSEMBL")
# TSS and gene length from the per-motif ChIPseeker files (one row per gene)
ck <- file.path(root, "human42_pipeline_full_local/ChIPseeker_output")
gi <- unique(rbindlist(lapply(list.files(ck, "\\.txt$", full.names = TRUE)[1:6], function(f)
  fread(f, sep = " ", select = c("geneChr", "geneStart", "geneEnd", "geneLength", "geneStrand", "ENSEMBL")))), by = "ENSEMBL")
gi <- gi[ENSEMBL %in% V$ENSEMBL]
# make sure every gene has a TSS: read more files for the missing ones
miss <- setdiff(V$ENSEMBL, gi$ENSEMBL)
for (f in list.files(ck, "\\.txt$", full.names = TRUE)[-(1:6)]) {
  if (!length(miss)) break
  x <- fread(f, sep = " ", select = c("geneChr", "geneStart", "geneEnd", "geneLength", "geneStrand", "ENSEMBL"))[ENSEMBL %in% miss]
  gi <- rbind(gi, unique(x, by = "ENSEMBL")); miss <- setdiff(miss, gi$ENSEMBL)
}
cat("genes without TSS info:", length(miss), "\n")
gi[, tss := ifelse(geneStrand == 1, geneStart, geneEnd)]
gi[, chrom := ifelse(grepl("^chr", geneChr), geneChr, paste0("chr", geneChr))]
gc <- fread(file.path(root, "human42_pipeline/clustering/gc_bins_1kb.tsv"))
gcl <- split(gc, gc$chrom)
pgc <- function(ch, tss, strand) {
  a <- max(0, if (strand == 1) tss - 3000 else tss - 1000); b <- if (strand == 1) tss + 1000 else tss + 3000
  g <- gcl[[ch]]; if (is.null(g)) return(NA_real_)
  w <- g[bin >= a %/% 1000 & bin <= b %/% 1000]; if (!nrow(w)) return(NA_real_)
  sum(w$gc) / sum(w$length - w$n)
}
gi[, promoter_gc := mapply(pgc, chrom, tss, geneStrand)]
D <- merge(merge(V, gi[, .(ENSEMBL, promoter_gc, gene_len = geneLength)], by = "ENSEMBL"), cov, by = "ENSEMBL")
D <- D[!is.na(promoter_gc)]
D[, `:=`(loglen = log10(gene_len), lexpr = log1p(expr))]
cat("genes analysed:", nrow(D), "\n")
cat("correlation of max score with covariates:\n"); print(round(cor(D[, .(s = `raw|d3000|max`, promoter_gc, loglen, lexpr)], method = "spearman")[1, -1], 3))

go2g <- function(id) { e <- unique(unlist(as.list(org.Hs.egGO2ALLEGS[id]))); ens <- unique(unlist(as.list(org.Hs.egENSEMBL[intersect(e, mappedkeys(org.Hs.egENSEMBL))]))); intersect(ens, D$ENSEMBL) }
sets <- list("fatty acid beta-oxidation" = go2g("GO:0006635"), "fatty acid catabolic process" = go2g("GO:0009062"), "lipid catabolic process" = go2g("GO:0016042"))
D[, strat := paste(cut(promoter_gc, quantile(promoter_gc, 0:5 / 5), include.lowest = TRUE, labels = FALSE), cut(loglen, quantile(loglen, 0:3 / 3), include.lowest = TRUE, labels = FALSE))]
auc <- function(x, y) { r <- rank(x); n1 <- sum(y); n0 <- sum(!y); (sum(r[y]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
strat_auc <- function(x, y, s) { tab <- data.table(x, y, s)[, .(a = if (sum(y) > 0 && sum(!y) > 0) auc(x, y) else NA_real_, w = sum(y) * sum(!y)), by = s][!is.na(a)]; sum(tab$a * tab$w) / sum(tab$w) }
set.seed(1)
vars <- c("raw|d3000|max", "raw|d3000|sum", "raw|d3000|median", "nlp|d3000|max")
res <- rbindlist(lapply(names(sets), function(sn) {
  y <- D$ENSEMBL %in% sets[[sn]]
  rbindlist(lapply(vars, function(vn) {
    x <- D[[vn]]
    fit <- glm(y ~ scale(x) + scale(promoter_gc) + scale(loglen) + scale(lexpr), data = cbind(D, y = y, x = x), family = binomial)
    co <- summary(fit)$coefficients["scale(x)", ]
    sa <- strat_auc(x, y, D$strat)
    nul <- replicate(300, { yp <- y; for (s in unique(D$strat)) { i <- which(D$strat == s); yp[i] <- sample(y[i]) }; strat_auc(x, yp, D$strat) })
    data.table(set = sn, n_genes = sum(y), variant = vn, crude_AUC = auc(x, y), stratified_AUC = sa, strat_z = (sa - mean(nul)) / sd(nul), strat_p_perm = (sum(nul >= sa) + 1) / 301,
               glm_OR_per_SD = exp(co[1]), glm_p = co[4])
  }))
}))
print(res[, lapply(.SD, function(v) if (is.numeric(v)) round(v, 4) else v)])
out <- file.path(root, "human42_pipeline/ppre_score_validation"); fwrite(res, file.path(out, "ppre_score_confounder_control.csv"))
fwrite(D[, .(ENSEMBL, promoter_gc, gene_len, expr)], file.path(out, "gene_covariates_gc_length_expression.csv"))
