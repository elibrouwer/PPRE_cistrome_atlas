# Test each gene-scoring method (from 11) against independent, community-
# curated PPAR gene sets (from 10) using fgsea (Korotkevich et al., the
# standard GSEA implementation in R/Bioconductor -- the same enrichment
# framework this project already used elsewhere, e.g. Objective1_results/
# GSEA_nonzero_PPRE_score_*_GSEA_results.csv). Only methods that clear an
# explicit, pre-stated bar are then combined into a single consensus ranking
# via RobustRankAggreg (Kolde, Laur, Adler & Vilo 2012, Bioinformatics 28(4)
# -- rank-aggregation method purpose-built for combining several independent
# gene rankings, ~1000+ citations), rather than combining all candidate
# methods uncritically. The consensus list is then validated the same way, to
# check it actually outperforms (or at least matches) the best individual
# method -- the empirical "did aggregating help" check.

suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
  library(RobustRankAggreg)
  library(ggplot2)
})
source("./PPRE_threshold_sensitivity/00_config.R")

scores <- fread(file.path(out_results, "gene_scoring_methods.csv"))
gene_sets <- readRDS(file.path(out_data, "reference_gene_sets.rds"))

project_genes <- gene_sets[["PROJECT_CONFIRMED_PPARA_FAO_GENES"]]
gsea_sets <- gene_sets[setdiff(names(gene_sets), "PROJECT_CONFIRMED_PPARA_FAO_GENES")]

methods <- c(
  "sum_neglog10p", "max_neglog10p",
  "fisher_neglog10p_combined", "fisher_neglog10p_heterodimer",
  "fisher_neglog10p_PPAR_half", "fisher_neglog10p_RXR_half",
  "stouffer_neglog10p", "n_DR1_pairs"
)

set.seed(1)
make_stat <- function(col) {
  x <- scores[[col]]
  names(x) <- scores$gene
  # break exact ties (common for n_DR1_pairs, an integer 0-16) with negligible
  # jitter so fgsea's ranking is well-defined without altering real ordering
  x + rnorm(length(x), 0, 1e-9 * (max(x) - min(x) + 1e-9))
}

## ---- fgsea for every (method, external gene set) pair ----
gsea_rows <- list()
for (m in methods) {
  stat <- make_stat(m)
  res <- fgsea(pathways = gsea_sets, stats = stat, minSize = 5, maxSize = 500, eps = 0)
  res <- as.data.table(res)[, .(pathway, NES, pval, padj, size)]
  res[, method := m]
  gsea_rows[[m]] <- res
}
gsea_tbl <- rbindlist(gsea_rows)
setcolorder(gsea_tbl, c("method", "pathway", "NES", "pval", "padj", "size"))
fwrite(gsea_tbl, file.path(out_results, "gsea_validation_by_method.csv"))
cat("=== fgsea NES / padj, every method x every external gene set ===\n")
print(gsea_tbl[order(pathway, padj)], nrow = Inf)

## ---- project-confirmed genes (n=5, too small for fgsea's permutation null):
## report percentile rank directly instead ----
pct_rank <- function(col) {
  r <- rank(-scores[[col]], ties.method = "average")
  100 * (1 - (r - 1) / (nrow(scores) - 1))   # 100 = top of the list
}
confirmed_tbl <- data.table(gene = project_genes)
for (m in methods) confirmed_tbl[[m]] <- pct_rank(m)[match(project_genes, scores$gene)]
fwrite(confirmed_tbl, file.path(out_results, "project_confirmed_genes_percentile_by_method.csv"))
cat("\n=== Percentile rank (100=top) of the 5 project-confirmed PPARA FAO genes, per method ===\n")
print(confirmed_tbl)

## ---- select methods that pass an explicit bar: significant (padj<0.05) enrichment
## for at least one of the 3 substantial community gene sets (excluding the small
## SANDERSON set, which is under-powered for a strict FDR cutoff at n=15) ----
bar_sets <- c("PPAR_DR1_Q2", "KEGG_PPAR_SIGNALING_PATHWAY", "REACTOME_REGULATION_OF_LIPID_METABOLISM_BY_PPARALPHA")
passing <- gsea_tbl[pathway %in% bar_sets & padj < 0.05, unique(method)]
failing <- setdiff(methods, passing)
cat(sprintf("\nMethods passing the bar (padj<0.05 for >=1 of %s): %s\n",
            paste(bar_sets, collapse = ", "), paste(passing, collapse = ", ")))
cat(sprintf("Methods NOT passing (excluded from aggregation): %s\n", paste(failing, collapse = ", ")))
stopifnot(length(passing) >= 2)  # aggregation needs at least 2 independent rankings

## ---- RobustRankAggreg: combine the passing methods' full gene rankings into
## one consensus list ----
glist <- lapply(passing, function(m) scores$gene[order(-scores[[m]])])
rra <- aggregateRanks(glist = glist, N = nrow(scores))
setDT(rra)
setnames(rra, c("Name", "Score"), c("gene", "RRA_score"))
rra <- rra[order(RRA_score)]
rra[, consensus_rank := .I]
fwrite(rra, file.path(out_results, "consensus_enriched_gene_list_RRA.csv"))

## ---- validate the consensus list the same way ----
consensus_stat <- setNames(-rra$RRA_score + rnorm(nrow(rra), 0, 1e-12), rra$gene)  # RRA_score: lower=better -> flip sign
consensus_res <- fgsea(pathways = gsea_sets, stats = consensus_stat, minSize = 5, maxSize = 500, eps = 0)
consensus_res <- as.data.table(consensus_res)[, .(pathway, NES, pval, padj, size)]
consensus_res[, method := "RRA_consensus"]
fwrite(consensus_res, file.path(out_results, "gsea_validation_consensus.csv"))

cat("\n=== Consensus (RRA) list validation ===\n")
print(consensus_res)

## ---- comparison plot: NES per method (incl. consensus) per gene set, bar_sets only ----
plot_tbl <- rbindlist(list(gsea_tbl[pathway %in% bar_sets], consensus_res[pathway %in% bar_sets]), fill = TRUE)
plot_tbl[, significant := padj < 0.05]
p <- ggplot(plot_tbl, aes(x = method, y = NES, fill = significant)) +
  geom_col() +
  geom_hline(yintercept = 0, color = "grey40") +
  facet_wrap(~pathway, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c(`TRUE` = "#2C7BB6", `FALSE` = "grey70")) +
  coord_flip() +
  labs(x = NULL, y = "fgsea normalized enrichment score (NES)",
       title = "Which gene-scoring method best recovers community-curated PPAR gene sets?",
       subtitle = "blue = significant at padj<0.05; RRA_consensus = RobustRankAggreg combination of the passing methods") +
  theme_minimal(base_size = 10)
ggsave(file.path(out_figures, "fig5_method_validation_NES.png"), p, width = 9, height = 9, dpi = 150, bg = "white")

cat("\n=== Top 30 of the final consensus enriched gene list ===\n")
print(head(rra, 30))
message("\nDone. Written: gsea_validation_by_method.csv, project_confirmed_genes_percentile_by_method.csv, consensus_enriched_gene_list_RRA.csv, gsea_validation_consensus.csv, fig5_method_validation_NES.png")
