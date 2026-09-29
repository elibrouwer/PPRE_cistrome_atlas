# Assemble the final "results by cutoff and motif family" table, and the
# summary plots: hit counts, promoter overlap, gene-ranking stability, and
# the ChIP fold-enrichment tradeoff across thresholds.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})
source("./PPRE_threshold_sensitivity/00_config.R")

raw_family   <- fread(file.path(out_results, "raw_hits_per_family_by_threshold.csv"))
interval_sum <- fread(file.path(out_results, "deduplicated_interval_summary_by_threshold.csv"))
promo_sum    <- fread(file.path(out_results, "promoter_overlap_summary_by_threshold.csv"))
fam_contrib  <- fread(file.path(out_results, "family_contribution_by_threshold.csv"))
chip_val     <- fread(file.path(out_results, "chip_validation_PPARA_liver_by_threshold.csv"))
rank_stab    <- fread(file.path(out_results, "rank_stability_vs_baseline.csv"))

thr_label <- function(x) sprintf("p<=%.0e", x)

## ---- master results table: cutoff x family ----
# chip_val's n_total_intervals is already the per-(threshold, family) deduplicated
# interval count (same family/view filter as everywhere else in the pipeline), so
# it is reused here rather than merging the genome-wide-only interval_sum table
# (that one has no per-family breakdown -- kept as its own separate summary CSV).
master <- merge(raw_family, chip_val, by.x = c("threshold", "family"), by.y = c("threshold", "view"), all = TRUE)
setnames(master, c("n_raw_hits", "n_total_intervals"), c("n_raw_motif_calls", "n_deduplicated_intervals"))
setnames(master, "family", "motif_family")
fwrite(master, file.path(out_results, "MASTER_results_by_cutoff_and_family.csv"))

## ---- Plot 1: raw hits and deduplicated intervals per family, by threshold ----
p1_data <- melt(raw_family[, .(threshold, family, n_raw_hits)],
                 id.vars = c("threshold", "family"), value.name = "n")
p1_data[, threshold_lab := factor(thr_label(threshold), levels = thr_label(sort(p_thresholds, decreasing = TRUE)))]
p1 <- ggplot(p1_data, aes(x = threshold_lab, y = n, fill = family)) +
  geom_col(position = "dodge") +
  scale_y_log10(labels = comma) +
  labs(x = "FIMO p-value threshold", y = "raw motif calls (log scale)",
       title = "Raw FIMO hits per motif family, by threshold",
       subtitle = "JASPAR/HOCOMOCO/CIS-BP overlapping calls NOT yet deduplicated") +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")
ggsave(file.path(out_figures, "fig1_raw_hits_per_family_by_threshold.png"), p1, width = 8, height = 5.5, dpi = 150, bg = "white")

## ---- Plot 2: promoter-overlapping intervals and genes assigned, combined view ----
p2_data <- promo_sum[view == "combined"]
p2_data[, threshold_lab := factor(thr_label(threshold), levels = thr_label(sort(p_thresholds, decreasing = TRUE)))]
p2 <- ggplot(p2_data, aes(x = threshold_lab, y = n_genes_assigned, fill = gene_set)) +
  geom_col(position = "dodge") +
  geom_text(aes(label = comma(n_genes_assigned)), position = position_dodge(width = 0.9), vjust = -0.4, size = 3) +
  labs(x = "FIMO p-value threshold", y = "cardiomyocyte-expressed genes with a promoter-overlapping candidate PPRE",
       title = "Genes assigned a candidate PPRE, by threshold",
       subtitle = "any motif family (combined view); nonzero = any detectable expression, top20 = top 20% by mean expression") +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")
ggsave(file.path(out_figures, "fig2_genes_assigned_by_threshold.png"), p2, width = 8, height = 5.5, dpi = 150, bg = "white")

## ---- Plot 3: rank stability vs. p<=1e-4 baseline ----
p3_data <- rank_stab[gene_set == "nonzero"]
p3_data[, threshold_lab := factor(thr_label(threshold_compare), levels = thr_label(sort(setdiff(p_thresholds, baseline_threshold), decreasing = TRUE)))]
p3 <- ggplot(p3_data, aes(x = threshold_lab, y = spearman_rho, fill = view)) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  labs(x = sprintf("compared to p<=%.0e baseline", baseline_threshold), y = "Spearman rank correlation",
       title = "Gene-ranking stability vs. the current p<=1e-4 analysis",
       subtitle = "nonzero-expressed gene universe; 1.0 = identical ranking") +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")
ggsave(file.path(out_figures, "fig3_rank_stability_vs_baseline.png"), p3, width = 8, height = 5.5, dpi = 150, bg = "white")

## ---- Plot 4: ChIP fold enrichment vs. threshold (indirect, liver PPARA) ----
chip_val[, threshold_lab := factor(thr_label(threshold), levels = thr_label(sort(p_thresholds, decreasing = TRUE)))]
p4 <- ggplot(chip_val, aes(x = threshold_lab, y = fold_enrichment, fill = view)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  geom_text(aes(label = comma(n_total_intervals)), position = position_dodge(width = 0.8), vjust = -0.4, size = 2.6) +
  labs(x = "FIMO p-value threshold", y = "fold enrichment in ChIP-Atlas liver PPARA peaks\n(vs. genome-wide chance)",
       title = "Genome-wide fold enrichment vs. real PPARA ChIP-seq, by threshold and family",
       subtitle = "INDIRECT validation: liver tissue, PPARA antigen only -- not cardiomyocyte, not PPARD/PPARG/RXR ChIP.\nlabels = deduplicated intervals kept at that threshold/family") +
  theme_minimal(base_size = 10) + theme(legend.position = "bottom")
ggsave(file.path(out_figures, "fig4_chip_fold_enrichment_by_threshold.png"), p4, width = 9, height = 6, dpi = 150, bg = "white")

cat("Wrote MASTER_results_by_cutoff_and_family.csv and 4 figures to\n  ", out_results, "\n  ", out_figures, "\n")
