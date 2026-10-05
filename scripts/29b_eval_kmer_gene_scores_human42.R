# Gene-level evaluation of the 6-mer sequence model scores saved by 29 (the job was interrupted after the tile predictions were written).
BENCH <- "C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Minor internship/R_code/benchmark"
source(file.path(BENCH, "config.R")); source(file.path(BENCH, "helpers_hits.R")); source(file.path(BENCH, "helpers_gold.R"))
out <- "C:/Users/brouw/PPRE_cistrome_atlas/results/human42_pipeline/score_exploration"
gi <- load_genes(); S <- readRDS(file.path(out, "kmer_gene_scores.rds"))
Sub <- readRDS("C:/Users/brouw/PPRE_cistrome_atlas/results/human42_pipeline/distal_kernels/motif_subset_scores.rds"); C <- readRDS("C:/Users/brouw/PPRE_cistrome_atlas/results/human42_pipeline/chip_calibration/chip_calibrated_gene_scores.rds")
J <- readRDS(file.path(res_dir, "J_candidate_gene_table.rds"))[match(gi$ENSEMBL, ENSEMBL)]
S$ref_PPAR14 <- Sub$PPAR_containing_14; S$ref_het3 <- Sub$heterodimer_3; S$ref_chip_top3 <- C$chip_model_top3_mean; S$ref_all42 <- Sub$all_42; S$expression <- ifelse(is.na(J$log10_expression), log10(1e-4), J$log10_expression)
rk <- function(x) frank(x) / length(x); S$kmer_top3_plus_expr <- 0.5 * rk(S$kmer_top3) + 0.5 * rk(S$expression)
reg <- readRDS(file.path(out, "gold_registry.rds")); reg <- reg[lengths(reg) >= 15]
ev <- NULL; for (s in 1:40) { set.seed(1700 + s); ev <- tryCatch(make_evaluator(gi, reg, ndraw = 200L), error = function(e) NULL); if (!is.null(ev)) break }; stopifnot(!is.null(ev))
Rg <- rbindlist(lapply(names(S), function(nm) ev(S[[nm]], nm)$table)); Rg[, family := substr(gold_set, 1, 1)]; fwrite(Rg, file.path(out, "kmer_gene_level_auroc_by_goldset.csv"))
fam <- Rg[, .(mean_AUROC = mean(auroc_matched)), by = .(score, family)]; Wf <- dcast(fam, score ~ family, value.var = "mean_AUROC"); Wf[, overall := rowMeans(.SD), .SDcols = c("A", "B", "C", "D", "E")]
setnames(Wf, c("A", "B", "C", "D", "E"), c("A literature", "B PPARgene+heart", "C promoter-peak", "D distal-peak", "E PPARD induced"))
print(Wf[order(-overall), lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)]); fwrite(Wf, file.path(out, "kmer_gene_level_family_means.csv"))
print(Rg[family == "C", .(gold_set, score, auroc = round(auroc_matched, 3))][score %in% c("kmer_top3", "ref_PPAR14", "expression")][order(gold_set, score)], nrows = 50)
