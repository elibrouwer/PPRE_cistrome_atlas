# Histogram of FIMO motif hits per gene (all 30 motifs, combined view,
# p<=1e-4 baseline threshold, cardiomyocyte-expressed gene universe).
#
# x-axis = number of promoter-overlapping motif-hit intervals for a gene
# (n_intervals from 05_promoter_overlap_gene_ranking.R), y-axis = number of
# genes with that many hits.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})
source("./PPRE_threshold_sensitivity/00_config.R")

gene_scores <- fread(file.path(out_results, "gene_scores_by_threshold.csv"))
d <- gene_scores[threshold == baseline_threshold & view == "combined" & gene_set == "nonzero"]

cat(sprintf("Genes: %s | hits per gene: min=%d, median=%.0f, mean=%.1f, max=%d\n",
            format(nrow(d), big.mark = ","), min(d$n_intervals), median(d$n_intervals),
            mean(d$n_intervals), max(d$n_intervals)))

# Project-confirmed PPARA fatty-acid-oxidation marker genes (same set used in
# 10_reference_gene_sets.R / 16_ml_and_similarity_scoring.R's project_genes)
marker_genes <- c("CPT1A", "CPT1B", "HADHA", "HADHB", "ACADVL")
markers <- d[gene %in% marker_genes, .(gene, n_intervals)]
missing_markers <- setdiff(marker_genes, markers$gene)
if (length(missing_markers) > 0)
  warning("Marker gene(s) not found in this threshold/view/gene_set subset: ", paste(missing_markers, collapse = ", "))
setorder(markers, n_intervals)

binwidth <- 5
max_count <- max(hist(d$n_intervals, breaks = seq(0, max(d$n_intervals) + binwidth, by = binwidth), plot = FALSE)$counts)
# stagger label heights top-to-bottom in x order so nearby markers don't overlap
markers[, y_label := max_count * seq(1.22, 0.9, length.out = .N)[frank(n_intervals, ties.method = "first")]]

p <- ggplot(d, aes(x = n_intervals)) +
  geom_histogram(binwidth = binwidth, fill = "#3B7FB6", color = "white", linewidth = 0.2) +
  geom_vline(data = markers, aes(xintercept = n_intervals), color = "#D6604D", linetype = "dashed", linewidth = 0.5) +
  geom_point(data = markers, aes(x = n_intervals, y = y_label), color = "#B2182B", size = 1.2) +
  geom_text(data = markers, aes(x = n_intervals, y = y_label, label = sprintf("%s (%d)", gene, n_intervals)),
            color = "#B2182B", size = 3, hjust = -0.08, vjust = 0.35) +
  scale_x_continuous(breaks = pretty_breaks(10)) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.22))) +
  coord_cartesian(clip = "off") +
  labs(x = "FIMO motif hits per gene (promoter-overlapping intervals, 30 motifs combined)",
       y = "number of genes",
       title = "Distribution of FIMO motif hits per gene") +
  theme_minimal(base_size = 11)

ggsave(file.path(out_figures, "fig7_hits_per_gene_histogram.png"), p, width = 8, height = 5.5, dpi = 150, bg = "white")
fwrite(d[, .(gene, n_hits = n_intervals)][order(-n_hits)],
       file.path(out_results, "hits_per_gene_combined_baseline.csv"))
message("Wrote fig7_hits_per_gene_histogram.png and hits_per_gene_combined_baseline.csv")
