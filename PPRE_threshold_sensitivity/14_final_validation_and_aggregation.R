# Re-run the method validation + RobustRankAggreg consensus (from 12), now
# using the size-normalized Poisson enrichment scores (13) instead of the
# raw/Fisher/Stouffer scores from 11, which turned out to be dominated by
# promoter-window-size variance (spearman rho=0.862 with promoter_bp; see 13).
# Both the naive (11) and corrected (13) results are validated side by side
# against the same community gene sets so the improvement is quantified, not
# just asserted.

suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
  library(RobustRankAggreg)
  library(ggplot2)
})
source("./PPRE_threshold_sensitivity/00_config.R")

scores_naive <- fread(file.path(out_results, "gene_scoring_methods.csv"))
scores_norm  <- fread(file.path(out_results, "gene_scoring_size_normalized.csv"))
dr1          <- fread(file.path(out_results, "dr1_hexamer_geometry_per_gene.csv"))
gene_sets    <- readRDS(file.path(out_data, "reference_gene_sets.rds"))
project_genes <- gene_sets[["PROJECT_CONFIRMED_PPARA_FAO_GENES"]]
gsea_sets    <- gene_sets[setdiff(names(gene_sets), "PROJECT_CONFIRMED_PPARA_FAO_GENES")]

scores <- Reduce(function(a,b) merge(a,b,by="gene",all.x=TRUE),
                  list(scores_naive[, .(gene, sum_neglog10p_NAIVE = sum_neglog10p)],
                       scores_norm,
                       dr1[, .(gene, n_DR1_pairs)]))
scores[, DR1_pairs_per_kb := n_DR1_pairs / (promoter_bp / 1000)]

methods <- c(
  "sum_neglog10p_NAIVE",                          # uncorrected baseline, for comparison only
  "poisson_neglog10p_combined", "poisson_neglog10p_heterodimer",
  "poisson_neglog10p_PPAR_half", "poisson_neglog10p_RXR_half",
  "DR1_pairs_per_kb"
)

set.seed(1)
make_stat <- function(col) {
  x <- scores[[col]]; names(x) <- scores$gene
  x + rnorm(length(x), 0, 1e-9 * (max(x, na.rm=TRUE) - min(x, na.rm=TRUE) + 1e-9))
}

gsea_rows <- list()
for (m in methods) {
  stat <- make_stat(m)
  res <- as.data.table(fgsea(pathways = gsea_sets, stats = stat, minSize = 5, maxSize = 500,
                              eps = 0, scoreType = "pos"))[, .(pathway, NES, pval, padj, size)]
  res[, method := m]
  gsea_rows[[m]] <- res
}
gsea_tbl <- rbindlist(gsea_rows)
setcolorder(gsea_tbl, c("method","pathway","NES","pval","padj","size"))
fwrite(gsea_tbl, file.path(out_results, "gsea_validation_size_normalized.csv"))
cat("=== fgsea NES/padj: naive (uncorrected) vs. size-normalized Poisson methods ===\n")
print(gsea_tbl[order(pathway, padj)], nrow = Inf)

## ---- project-confirmed genes: percentile rank, naive vs. corrected ----
pct_rank <- function(col) { r <- rank(-scores[[col]], ties.method="average", na.last="keep")
  100 * (1 - (r-1)/(nrow(scores)-1)) }
confirmed_tbl <- data.table(gene = project_genes)
for (m in methods) confirmed_tbl[[m]] <- pct_rank(m)[match(project_genes, scores$gene)]
fwrite(confirmed_tbl, file.path(out_results, "project_confirmed_genes_percentile_FINAL.csv"))
cat("\n=== Percentile rank (100=top) of project-confirmed genes, naive vs. corrected ===\n")
print(confirmed_tbl)

## ---- select methods passing the bar (padj<0.05 for >=1 of the 3 substantial sets),
## EXCLUDING the naive baseline by design (kept only for the comparison table above,
## not eligible for the consensus list since it's the flawed method) ----
bar_sets <- c("PPAR_DR1_Q2", "KEGG_PPAR_SIGNALING_PATHWAY", "REACTOME_REGULATION_OF_LIPID_METABOLISM_BY_PPARALPHA")
eligible <- setdiff(methods, "sum_neglog10p_NAIVE")
passing <- gsea_tbl[method %in% eligible & pathway %in% bar_sets & padj < 0.05, unique(method)]
cat(sprintf("\nCorrected methods passing the bar: %s\n", paste(passing, collapse=", ")))
cat(sprintf("Corrected methods NOT passing: %s\n", paste(setdiff(eligible, passing), collapse=", ")))
stopifnot(length(passing) >= 2)

glist <- lapply(passing, function(m) scores$gene[order(-scores[[m]])])
rra <- setDT(aggregateRanks(glist = glist, N = nrow(scores)))
setnames(rra, c("Name","Score"), c("gene","RRA_score"))
rra <- rra[order(RRA_score)][, consensus_rank := .I]
fwrite(rra, file.path(out_results, "consensus_enriched_gene_list_FINAL.csv"))

consensus_stat <- setNames(-rra$RRA_score + rnorm(nrow(rra), 0, 1e-12), rra$gene)
consensus_res <- as.data.table(fgsea(pathways = gsea_sets, stats = consensus_stat, minSize = 5,
                                      maxSize = 500, eps = 0, scoreType = "pos"))[, .(pathway, NES, pval, padj, size)]
consensus_res[, method := "RRA_consensus_FINAL"]
fwrite(consensus_res, file.path(out_results, "gsea_validation_consensus_FINAL.csv"))
cat("\n=== FINAL consensus (size-corrected RRA) list validation ===\n"); print(consensus_res)

## ---- where do the project-confirmed genes land in the final consensus? ----
rra[, is_confirmed := gene %in% project_genes]
cat("\n=== Project-confirmed genes' rank in the final consensus list (of 15,703) ===\n")
print(rra[gene %in% project_genes])

## ---- comparison plot ----
plot_tbl <- rbindlist(list(gsea_tbl[pathway %in% bar_sets], consensus_res[pathway %in% bar_sets]), fill = TRUE)
plot_tbl[, significant := padj < 0.05]
plot_tbl[, method := factor(method, levels = c("sum_neglog10p_NAIVE", eligible, "RRA_consensus_FINAL"))]
p <- ggplot(plot_tbl, aes(x = method, y = NES, fill = significant)) +
  geom_col() + geom_hline(yintercept = 0, color = "grey40") +
  facet_wrap(~pathway, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c(`TRUE`="#2C7BB6", `FALSE`="grey70")) +
  coord_flip() +
  labs(x = NULL, y = "fgsea NES", title = "Naive (size-biased) vs. size-normalized scoring methods",
       subtitle = "blue = padj<0.05; NAIVE = uncorrected baseline kept for comparison only, excluded from the consensus") +
  theme_minimal(base_size = 10)
ggsave(file.path(out_figures, "fig6_final_method_validation_NES.png"), p, width = 9, height = 9, dpi = 150, bg = "white")

cat("\n=== Top 30 of the FINAL, size-corrected consensus enriched gene list ===\n")
print(head(rra, 30))
message("\nDone.")
