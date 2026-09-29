# Fair test of the nearest-neighbor "comparables" method from 16: the
# centroid was built FROM the 5 known genes and then those same 5 genes were
# checked against it -- with n=5 defining a 13-dimensional centroid, ranking
# those 5 genes highly is close to guaranteed by construction (circularity),
# not evidence the method generalizes to an unseen gene. Leave-one-out fixes
# this: for each of the 5 genes, build the centroid from the OTHER 4 only,
# then check where the held-out gene lands.

suppressPackageStartupMessages(library(data.table))
source("./PPRE_threshold_sensitivity/00_config.R")

dat <- fread(file.path(out_results, "ml_and_similarity_gene_scores.csv"))
scores <- fread(file.path(out_results, "gene_scoring_size_normalized.csv"))
dr1    <- fread(file.path(out_results, "dr1_hexamer_geometry_per_gene.csv"))[, .(gene, n_DR1_pairs)]
expr   <- fread(expr_nonzero_csv)[, .(gene = HGNC_gene_name, expression_value)]
expr   <- unique(expr, by = "gene")

full <- Reduce(function(a,b) merge(a,b,by="gene",all.x=TRUE), list(scores, dr1, expr))
full[is.na(n_DR1_pairs), n_DR1_pairs := 0]
full[, DR1_pairs_per_kb := n_DR1_pairs / (promoter_bp/1000)]
full[is.na(expression_value), expression_value := 0]

feature_cols <- c("fisher_n_loci_combined","fisher_n_loci_heterodimer","fisher_n_loci_PPAR_half","fisher_n_loci_RXR_half",
                   "promoter_bp","n_windows",
                   "poisson_fold_combined","poisson_fold_heterodimer","poisson_fold_PPAR_half","poisson_fold_RXR_half",
                   "n_DR1_pairs","DR1_pairs_per_kb","expression_value")

z_mat <- scale(as.matrix(full[, ..feature_cols]))
rownames(z_mat) <- full$gene
project_genes <- c("CPT1A","CPT1B","HADHA","HADHB","ACADVL")

pct_rank <- function(x) { r <- rank(-x, ties.method="average"); 100*(1-(r-1)/(length(x)-1)) }

loo_results <- data.table(gene = project_genes, held_out_percentile = NA_real_, all5_percentile = NA_real_)

# for comparison: the (circular) all-5 percentile, recomputed here identically to 16
centroid_all5 <- colMeans(z_mat[rownames(z_mat) %in% project_genes, , drop = FALSE])
dist_all5 <- sqrt(rowSums(sweep(z_mat, 2, centroid_all5, "-")^2))
sim_all5 <- -dist_all5[match(full$gene, rownames(z_mat))]
loo_results$all5_percentile <- pct_rank(sim_all5)[match(project_genes, full$gene)]

for (g in project_genes) {
  others <- setdiff(project_genes, g)
  centroid_loo <- colMeans(z_mat[rownames(z_mat) %in% others, , drop = FALSE])
  dist_loo <- sqrt(rowSums(sweep(z_mat, 2, centroid_loo, "-")^2))
  sim_loo <- -dist_loo[match(full$gene, rownames(z_mat))]
  pct <- pct_rank(sim_loo)
  loo_results[gene == g, held_out_percentile := pct[match(g, full$gene)]]
}

cat("=== Leave-one-out fairness check: similarity-to-known-positives method ===\n")
cat("'all5_percentile' = circular (built from all 5, including itself) -- optimistic.\n")
cat("'held_out_percentile' = fair (centroid built from the OTHER 4 only) -- the real test.\n\n")
print(loo_results)

cat(sprintf("\nMean percentile, circular (all-5) version: %.1f\n", mean(loo_results$all5_percentile)))
cat(sprintf("Mean percentile, fair leave-one-out version: %.1f\n", mean(loo_results$held_out_percentile)))

fwrite(loo_results, file.path(out_results, "similarity_method_leave_one_out_check.csv"))
message("\nDone.")
