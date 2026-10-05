# Literature-curated human PPREs (literature_PPRE_human_prioritised.csv; tier 1 = strongest evidence, tier 2 = weaker/abstract-only) as gold standard for the 42-motif score.
# (a) gene level: genes carrying a tier-1 (and tier 1+2) PPRE; promoter scores vs promoter-length/GC matched genes (matched-background AUROC, as 16b-22).
# (b) site level: motif evidence (best relative FIMO score in the 500 bp bin(s) around the validated PPRE, TSS-relative, all TSS of the gene) vs the same relative coordinates placed at random genes.
# Rows used: ref_class == TSS and numbering within -3 kb/+1 kb (window_final == yes). Gene symbols are matched to the universe incl. the alias in parentheses.
BENCH <- "C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Minor internship/R_code/benchmark"
source(file.path(BENCH, "config.R")); source(file.path(BENCH, "helpers_hits.R")); source(file.path(BENCH, "helpers_gold.R"))
set.seed(25)
out <- "C:/Users/brouw/PPRE_cistrome_atlas/results/human42_pipeline/literature_ppre"; dir.create(out, showWarnings = FALSE)
lit <- fread("C:/Users/brouw/Downloads/literature_PPRE_human_prioritised.csv")
lit <- lit[tier %in% 1:2 & is_human == "yes"]
gi <- load_genes(); tss <- readRDS(file.path(cache_dir, "tss_manuscript.rds"))
sym_candidates <- function(s) unique(toupper(trimws(unlist(strsplit(gsub("[()/,]", " ", s), "\\s+")))))
lit[, row := .I]
map <- rbindlist(lapply(seq_len(nrow(lit)), function(i) { cand <- sym_candidates(lit$Gene_symbol[i]); m <- gi[toupper(gi$SYMBOL) %in% cand]; if (!nrow(m)) return(NULL); data.table(row = i, ENSEMBL = m$ENSEMBL[1]) }))
lit <- merge(lit, map, by = "row", all.x = TRUE)
cat("rows tier1+2:", nrow(lit), " mapped to universe:", sum(!is.na(lit$ENSEMBL)), "\n"); print(unique(lit[is.na(ENSEMBL), .(tier, Gene_symbol)]))
lit[, rel := ifelse(is.na(end_num), start_num, (start_num + end_num) / 2)]
use <- lit[!is.na(ENSEMBL) & window_final == "yes" & ref_class == "TSS" & !is.na(rel) & rel >= -3000 & rel <= 1000]
cat("sites used (in window, TSS-referenced):", nrow(use), " genes:", uniqueN(use$ENSEMBL), " tier1 genes:", uniqueN(use[tier == 1, ENSEMBL]), " tier1+2 genes:", uniqueN(use$ENSEMBL), "\n")
gold <- list(`tier 1` = unique(use[tier == 1, ENSEMBL]), `tier 1+2` = unique(use$ENSEMBL))
gold$`tier 1 + 2, PPARG only` <- unique(use[toupper(Subtype) == "PPARG", ENSEMBL]); gold$`tier 1 + 2, PPARA/PPARD only` <- unique(use[toupper(Subtype) %in% c("PPARA", "PPARD"), ENSEMBL])
print(sapply(gold, length))

# ---- (a) gene-level ----
K <- readRDS("C:/Users/brouw/PPRE_cistrome_atlas/results/human42_pipeline/distal_kernels/distal_kernel_scores.rds")
Sub <- readRDS("C:/Users/brouw/PPRE_cistrome_atlas/results/human42_pipeline/distal_kernels/motif_subset_scores.rds")
C <- readRDS("C:/Users/brouw/PPRE_cistrome_atlas/results/human42_pipeline/chip_calibration/chip_calibrated_gene_scores.rds")
J <- readRDS(file.path(res_dir, "J_candidate_gene_table.rds"))[match(gi$ENSEMBL, ENSEMBL)]
S <- list(all_42 = Sub$all_42, PPAR_containing_14 = Sub$PPAR_containing_14, heterodimer_3 = Sub$heterodimer_3, RXR_only_28 = Sub$RXR_only_28, chip_model_top3_mean = C$chip_model_top3_mean,
          `promoter GC` = gi$gc, `cardiomyocyte expression` = ifelse(is.na(J$log10_expression), log10(1e-4), J$log10_expression))
ev <- NULL; for (s in 1:30) { set.seed(500 + s); ev <- tryCatch(make_evaluator(gi, gold, ndraw = 300L), error = function(e) NULL); if (!is.null(ev)) break }
R <- rbindlist(lapply(names(S), function(nm) ev(S[[nm]], nm)$table)); fwrite(R, file.path(out, "literature_ppre_gene_level_auroc.csv"))
print(dcast(R, gold_set + n_gold ~ score, value.var = "auroc_matched")[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])
bc <- NULL; for (s in 1:30) { set.seed(700 + s); bc <- tryCatch(boot_compare(S[c("all_42", "PPAR_containing_14", "heterodimer_3")], gold[1:2], gi, ND = 100L, B = 500L), error = function(e) NULL); if (!is.null(bc)) break }
if (!is.null(bc)) { f <- function(a, b, st) { x <- bc$ba[, a, st] - bc$ba[, b, st]; data.table(gold = st, a = a, b = b, diff = mean(x), lo = unname(quantile(x, .025)), hi = unname(quantile(x, .975))) }
  cb <- rbindlist(lapply(names(gold)[1:2], function(st) rbind(f("PPAR_containing_14", "all_42", st), f("heterodimer_3", "all_42", st), f("heterodimer_3", "PPAR_containing_14", st))))
  print(cb[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)]); fwrite(cb, file.path(out, "literature_ppre_gene_level_paired_bootstrap.csv")) }

# ---- (b) site level ----
B <- readRDS("C:/Users/brouw/PPRE_cistrome_atlas/results/human42_pipeline_full_local/fimo_human42_bins500.rds"); setkey(B, chr, bin)
U <- unique(tss[in_universe == TRUE & ENSEMBL %in% gi$ENSEMBL, .(ENSEMBL, chr, strand, tss)]); U[, sgn := ifelse(strand == "+", 1L, -1L)]; setkey(U, ENSEMBL)
stat_at <- function(ens, rel) {   # max over the gene's TSSs of the best bin value around the TSS-relative position
  u <- U[.(ens), nomatch = 0L]; if (!nrow(u)) return(c(NA_real_, NA_real_))
  pos <- u$tss + u$sgn * as.integer(round(rel)); bins <- c(pos %/% 500L - 0L, pos %/% 500L, pos %/% 500L)   # bin containing the position
  q <- rbindlist(lapply(c(-1L, 0L, 1L), function(o) data.table(chr = u$chr, bin = pos %/% 500L + o)))   # +-1 bin = about +-500 bp around the site
  v <- B[q, nomatch = 0L]; if (!nrow(v)) return(c(0, 0)); c(max(v$max_rel), max(v$het_rel)) }
use[, ui := .I]; obs <- t(mapply(stat_at, use$ENSEMBL, use$rel)); use[, `:=`(obs_all = obs[, 1], obs_het = obs[, 2])]
univ <- unique(U$ENSEMBL); NNULL <- 300
nul <- lapply(seq_len(nrow(use)), function(i) { g <- sample(setdiff(univ, use$ENSEMBL[i]), NNULL); t(vapply(g, function(e) stat_at(e, use$rel[i]), numeric(2))) })
use[, `:=`(pct_all = vapply(seq_len(.N), function(i) mean(nul[[i]][, 1] < obs_all[i], na.rm = TRUE) + 0.5 * mean(nul[[i]][, 1] == obs_all[i], na.rm = TRUE), numeric(1)),
           pct_het = vapply(seq_len(.N), function(i) mean(nul[[i]][, 2] < obs_het[i], na.rm = TRUE) + 0.5 * mean(nul[[i]][, 2] == obs_het[i], na.rm = TRUE), numeric(1)))]
sm <- function(d, lab) data.table(set = lab, n_sites = nrow(d), mean_pct_all = mean(d$pct_all), top25_all = mean(d$pct_all >= .75), mean_pct_het = mean(d$pct_het), top25_het = mean(d$pct_het >= .75),
  median_obs_het = median(d$obs_het), median_null_het = median(unlist(lapply(nul[d$ui], function(m) m[, 2]))))
res <- rbind(sm(use[tier == 1], "tier 1"), sm(use, "tier 1+2"), sm(use[toupper(Subtype) == "PPARG"], "tier 1+2 PPARG"), sm(use[toupper(Subtype) %in% c("PPARA", "PPARD")], "tier 1+2 PPARA/PPARD"))
print(res[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)]); fwrite(res, file.path(out, "literature_ppre_site_level_summary.csv"))
fwrite(use[, .(tier, Gene_symbol, Subtype, PMID, rel, obs_all, obs_het, pct_all, pct_het, evidence_score)], file.path(out, "literature_ppre_site_level_per_site.csv"))
# null expectation of the summary: mean percentile should be 0.5 and top-25% fraction 0.25 if there is no signal
