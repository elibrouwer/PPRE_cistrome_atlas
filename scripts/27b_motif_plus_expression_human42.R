# Idea 2: combine promoter motif evidence with cardiomyocyte expression (GTEx snRNA-seq log10 expression) and ask whether the combination beats expression alone (incremental value) and motif alone.
# Combination = weighted rank average, w = weight of the motif score (0.25 / 0.5 / 0.75). Evaluated on the gold registry (27a) with the matched-background AUROC (promoter length / GC) and paired bootstraps per gold family.
BENCH <- "C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Minor internship/R_code/benchmark"
source(file.path(BENCH, "config.R")); source(file.path(BENCH, "helpers_hits.R")); source(file.path(BENCH, "helpers_gold.R"))
out <- "C:/Users/brouw/PPRE_cistrome_atlas/results/human42_pipeline/score_exploration"
gi <- load_genes(); reg <- readRDS(file.path(out, "gold_registry.rds")); reg <- reg[lengths(reg) >= 15]
Sub <- readRDS("C:/Users/brouw/PPRE_cistrome_atlas/results/human42_pipeline/distal_kernels/motif_subset_scores.rds"); C <- readRDS("C:/Users/brouw/PPRE_cistrome_atlas/results/human42_pipeline/chip_calibration/chip_calibrated_gene_scores.rds")
J <- readRDS(file.path(res_dir, "J_candidate_gene_table.rds"))[match(gi$ENSEMBL, ENSEMBL)]
expr <- ifelse(is.na(J$log10_expression), log10(1e-4), J$log10_expression)
rk <- function(x) frank(x) / length(x)
S <- list(PPAR14 = Sub$PPAR_containing_14, het3 = Sub$heterodimer_3, all42 = Sub$all_42, chip_top3 = C$chip_model_top3_mean, expression = expr)
for (m in c("PPAR14", "het3", "chip_top3")) for (w in c(0.25, 0.5, 0.75)) S[[sprintf("%s+expr (w=%.2f)", m, w)]] <- w * rk(S[[m]]) + (1 - w) * rk(expr)
ev <- NULL; for (s in 1:40) { set.seed(1100 + s); ev <- tryCatch(make_evaluator(gi, reg, ndraw = 200L), error = function(e) NULL); if (!is.null(ev)) break }; stopifnot(!is.null(ev))
R <- rbindlist(lapply(names(S), function(nm) ev(S[[nm]], nm)$table)); R[, family := substr(gold_set, 1, 1)]; fwrite(R, file.path(out, "motif_plus_expression_auroc_by_goldset.csv"))
fam <- R[, .(mean_AUROC = mean(auroc_matched)), by = .(score, family)]; Wf <- dcast(fam, score ~ family, value.var = "mean_AUROC"); Wf[, overall := rowMeans(.SD), .SDcols = c("A", "B", "C", "D", "E")]
setnames(Wf, c("A", "B", "C", "D", "E"), c("A literature", "B PPARgene+heart", "C promoter-peak", "D distal-peak", "E PPARD induced"))
print(Wf[order(-overall), lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)]); fwrite(Wf, file.path(out, "motif_plus_expression_family_means.csv"))
# paired bootstrap: combination vs expression alone, and vs motif alone (PPAR14), per family
fams <- split(names(reg), substr(names(reg), 1, 1)); rowsb <- list()
for (f in names(fams)) { G <- reg[fams[[f]]]; G <- G[lengths(G) >= 15]; if (!length(G)) next; bc <- NULL
  for (s in 1:40) { set.seed(1300 + s); bc <- tryCatch(boot_compare(S[c("PPAR14", "expression", "PPAR14+expr (w=0.50)", "het3", "het3+expr (w=0.50)", "chip_top3", "chip_top3+expr (w=0.50)")], G, gi, ND = 60L, B = 300L), error = function(e) NULL); if (!is.null(bc)) break }
  if (is.null(bc)) next
  cmp <- function(a, b) { x <- rowMeans(bc$ba[, a, , drop = FALSE] - bc$ba[, b, , drop = FALSE]); data.table(family = f, a = a, b = b, diff = mean(x), lo = unname(quantile(x, .025)), hi = unname(quantile(x, .975))) }
  rowsb[[f]] <- rbind(cmp("PPAR14+expr (w=0.50)", "expression"), cmp("PPAR14+expr (w=0.50)", "PPAR14"), cmp("het3+expr (w=0.50)", "expression"), cmp("chip_top3+expr (w=0.50)", "expression"), cmp("PPAR14", "expression")) }
cb <- rbindlist(rowsb); print(cb[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)], nrows = 40); fwrite(cb, file.path(out, "motif_plus_expression_paired_bootstrap.csv"))
