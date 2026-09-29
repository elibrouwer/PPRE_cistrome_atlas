# Gene-level contribution of half-site motifs vs. full heterodimer motifs.
#
# Interval-level heterodimer/half-site/multi-family counts are already in
# deduplicated_interval_summary_by_threshold.csv (from 04). This script adds
# the more biologically relevant gene-level question: of the genes assigned a
# candidate PPRE in their promoter (the "combined" view from 05), how many
# have that support from a full heterodimer motif at all, vs. only from a
# PPAR or RXR half-site motif with no heterodimer hit anywhere in the same
# promoter? Heterodimer motifs are expected to be much rarer hits than
# half-site motifs (a DR1 heterodimer match requires both half-sites AND the
# correct 1bp spacer simultaneously), so this quantifies how much of the
# candidate gene set would be lost if half-site-only evidence were excluded.

suppressPackageStartupMessages(library(data.table))
source("./PPRE_threshold_sensitivity/00_config.R")

gene_scores <- fread(file.path(out_results, "gene_scores_by_threshold.csv"))

rows <- list()
for (pt in p_thresholds) for (gs_name in unique(gene_scores$gene_set)) {
  het   <- gene_scores[threshold == pt & gene_set == gs_name & view == "heterodimer", gene]
  ppar  <- gene_scores[threshold == pt & gene_set == gs_name & view == "PPAR_half",   gene]
  rxr   <- gene_scores[threshold == pt & gene_set == gs_name & view == "RXR_half",    gene]
  comb  <- gene_scores[threshold == pt & gene_set == gs_name & view == "combined",    gene]
  half_only <- setdiff(comb, het)   # in the combined gene set but not via any heterodimer hit

  rows[[length(rows) + 1]] <- data.table(
    threshold = pt, gene_set = gs_name,
    n_genes_combined            = length(comb),
    n_genes_with_heterodimer    = length(het),
    n_genes_half_site_only      = length(half_only),
    frac_half_site_only         = if (length(comb) > 0) length(half_only) / length(comb) else NA_real_,
    n_genes_PPAR_half_any       = length(ppar),
    n_genes_RXR_half_any        = length(rxr),
    n_genes_PPAR_and_RXR_half   = length(intersect(ppar, rxr))
  )
}
family_contribution <- rbindlist(rows)
fwrite(family_contribution, file.path(out_results, "family_contribution_by_threshold.csv"))
cat("=== Gene-level contribution: heterodimer vs. half-site-only evidence ===\n")
print(family_contribution)
message("Done. Written to ", file.path(out_results, "family_contribution_by_threshold.csv"))
