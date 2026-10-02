# Are plain FIMO hit counts per gene (fig7-fig12 definitions, 42 motifs) as good as FIMO-score-based gene scores?
# Same gold sets and matched-background AUROC evaluator as 16b (./benchmark). Restricted to the 15,703 cardiomyocyte-expressed genes of the hit-count tables.
BENCH <- "C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Minor internship/R_code/benchmark"
source(file.path(BENCH, "config.R")); source(file.path(BENCH, "helpers_hits.R")); source(file.path(BENCH, "helpers_gold.R"))
hp <- "C:/Users/brouw/PPRE_cistrome_atlas/results/human42_pipeline/hits_per_gene"
out <- "C:/Users/brouw/PPRE_cistrome_atlas/results/human42_pipeline/distal_kernels"
gi0 <- load_genes()
tabs <- list(
  count_raw          = c("hits_per_gene_combined_baseline.csv", "n_hits"),
  count_corrected    = c("hits_per_gene_size_corrected.csv", "fold_corrected"),
  count_fixed_window = c("hits_per_gene_fixed_window.csv", "n_hits"),
  rep_count_raw      = c("hits_per_gene_repeat_filtered.csv", "n_hits"),
  rep_count_corrected = c("hits_per_gene_size_corrected_repeat_filtered.csv", "poisson_fold"),
  rep_count_fixed    = c("hits_per_gene_fixed_window_repeat_filtered.csv", "n_hits"))
D <- lapply(tabs, function(a) { d <- fread(file.path(hp, a[1])); setNames(d[[a[2]]], d$gene) })
genes <- Reduce(intersect, c(lapply(D, names), list(gi0$SYMBOL)))
gi <- gi0[match(genes, gi0$SYMBOL)]; gi <- gi[!duplicated(gi$SYMBOL)]
cat("genes evaluated:", nrow(gi), "\n")
S <- lapply(D, function(v) as.numeric(v[gi$SYMBOL]))
K <- readRDS(file.path(out, "distal_kernel_scores.rds")); idx <- match(gi$ENSEMBL, gi0$ENSEMBL)
S$score_prom_base <- K$prom_base[idx]; S$score_prom_het <- K$prom_het[idx]
gold <- load_gold_all(gi); gold <- gold[c(PG_SETS, HEART_SETS)]; gold <- lapply(gold, function(g) intersect(g, gi$ENSEMBL))
cat("gold sizes:", paste(names(gold), lengths(gold), collapse = "; "), "\n")
ev <- make_evaluator(gi, gold, ndraw = 300L)
R <- rbindlist(lapply(names(S), function(nm) ev(S[[nm]], nm)$table)); fwrite(R, file.path(out, "hit_count_vs_score_auroc_by_goldset.csv"))
W <- dcast(R, score ~ gold_set, value.var = "auroc_matched")
W[, `:=`(mean_PPARgene = rowMeans(as.matrix(.SD[, PG_SETS, with = FALSE])), mean_heart = rowMeans(as.matrix(.SD[, HEART_SETS, with = FALSE])))]
W[, overall := (mean_PPARgene + mean_heart) / 2]; W <- W[order(-overall)]; fwrite(W, file.path(out, "hit_count_vs_score_auroc_summary.csv"))
print(W[, .(score, mean_PPARgene = round(mean_PPARgene, 3), mean_heart = round(mean_heart, 3), overall = round(overall, 3))])
bc <- boot_compare(S[c("count_corrected", "count_raw", "score_prom_base", "score_prom_het", "rep_count_corrected")], gold, gi, ND = 100L, B = 500L)
cmp <- function(a, b) { f <- function(sets) { d <- rowMeans(bc$ba[, a, sets] - bc$ba[, b, sets]); c(mean = mean(d), lo = unname(quantile(d, .025)), hi = unname(quantile(d, .975))) }
  cbind(a = a, b = b, PPARgene = t(round(f(PG_SETS), 3)), heart = t(round(f(HEART_SETS), 3))) }
print(rbind(cmp("count_corrected", "score_prom_base"), cmp("count_corrected", "score_prom_het"), cmp("count_raw", "score_prom_base"), cmp("rep_count_corrected", "count_corrected")))
