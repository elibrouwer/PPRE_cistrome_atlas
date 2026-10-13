# Remaining manuscript figures for the 42-motif analysis: score distribution and known genes (Fig. 8 replacement), site count vs promoter length (Fig. 3A replacement),
# and PNG versions (300 dpi) of the rerun heatmaps, Venn, H3K27ac overlap and ORA plots.
suppressPackageStartupMessages({library(ggplot2); library(ggrepel); library(pdftools)})
base <- "C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Minor internship"
res <- file.path(base, "PPAR_RXR_motifs_Hocomoco_Jaspar", "Human42_analysis", "Objective1_results"); out <- file.path(base, "R_code", "Claude outputs", "PPRE_boxplots_png")
s <- read.csv(file.path(base, "R_code", "results", "PPRE_score_v4", "PPRE_gene_scores_sitecount.csv"), stringsAsFactors = FALSE)
known <- c("CPT1A", "CPT1B", "HADHA", "HADHB", "ACADVL"); recep <- c("PPARA", "PPARD", "PPARG", "RXRA", "RXRB", "RXRG")
s$grp <- ifelse(s$SYMBOL %in% known, "Known PPRE-responsive", ifelse(s$SYMBOL %in% recep, "PPAR / RXR gene", "other"))
cols <- c("Known PPRE-responsive" = "#D96C6C", "PPAR / RXR gene" = "#4F86D9")
th <- theme_classic() + theme(legend.title = element_blank())
sv <- function(p, n, w, h) { ggsave(file.path(out, paste0(n, ".png")), p, width = w, height = h, dpi = 300); ggsave(file.path(res, "Statistics_plots", "single_plots", paste0(n, ".pdf")), p, width = w, height = h) }
# Fig. 8A: distribution of the number of PPRE sites per promoter with the reference genes
ref <- s[s$grp != "other", ]
p8a <- ggplot(s, aes(n_sites)) + geom_histogram(binwidth = 2, fill = "grey80", colour = "white") + geom_vline(data = ref, aes(xintercept = n_sites, colour = grp), linewidth = 0.5) +
  geom_text_repel(data = ref, aes(x = n_sites, y = Inf, label = SYMBOL, colour = grp), angle = 90, vjust = 1, size = 2.6, segment.size = 0.2, direction = "x", show.legend = FALSE, nudge_y = 0) +
  scale_colour_manual(values = cols) + coord_cartesian(xlim = c(0, 80)) + labs(x = "PPRE sites per promoter (n)", y = "Genes (n)") + th + theme(legend.position = "bottom")
sv(p8a, "fig8A_sites_per_promoter_distribution", 6, 3.6)
# Fig. 8B: percentile (within promoter length x GC group) of the reference genes
r2 <- ref[order(ref$PPRE_score_percentile), ]; r2$SYMBOL <- factor(r2$SYMBOL, levels = r2$SYMBOL)
p8b <- ggplot(r2, aes(PPRE_score_percentile, SYMBOL, colour = grp)) + geom_vline(xintercept = 50, linetype = "dashed", colour = "grey50") + geom_point(size = 3) +
  scale_colour_manual(values = cols) + scale_x_continuous(limits = c(0, 100)) + labs(x = "PPRE score (percentile within promoter length x GC group)", y = NULL) + th + theme(legend.position = "bottom")
sv(p8b, "fig8B_known_genes_percentile", 5, 3.6)
# Fig. 3A: sites vs promoter length (the size dependence the percentile removes)
p3a <- ggplot(s, aes(promoter_length / 1000, n_sites + 1)) + geom_point(colour = "grey70", size = 0.4, alpha = 0.4) + geom_smooth(method = "gam", colour = "black", se = FALSE, linewidth = 0.6) +
  geom_point(data = ref, aes(colour = grp), size = 2) + geom_text_repel(data = ref, aes(label = SYMBOL, colour = grp), size = 2.6, show.legend = FALSE, max.overlaps = 30) +
  scale_colour_manual(values = cols) + scale_x_log10() + scale_y_log10() + labs(x = "Promoter length (kb, log scale)", y = "PPRE sites per promoter + 1 (log scale)") + th + theme(legend.position = "bottom")
sv(p3a, "fig3A_sites_vs_promoter_length", 6, 4)
cat(sprintf("known genes: median percentile %.1f (range %.1f-%.1f); receptor genes: median %.1f (range %.1f-%.1f)\n", median(ref$PPRE_score_percentile[ref$grp == "Known PPRE-responsive"]), min(ref$PPRE_score_percentile[ref$grp == "Known PPRE-responsive"]), max(ref$PPRE_score_percentile[ref$grp == "Known PPRE-responsive"]),
            median(ref$PPRE_score_percentile[ref$grp == "PPAR / RXR gene"]), min(ref$PPRE_score_percentile[ref$grp == "PPAR / RXR gene"]), max(ref$PPRE_score_percentile[ref$grp == "PPAR / RXR gene"])))
cat("Spearman n_sites vs promoter length:", cor(s$n_sites, s$promoter_length, method = "spearman"), "\n")
# PNG versions of the rerun PDFs
for (f in c("overlap_heatmap_42", "jaccard_heatmap_42", "venn_diagram_42", "H3K27ac_overlap_stackplot_42", "ORA_dotplot_CM_expressed_default_universe", "ORA_dotplot_all_genes_scored_default_universe")) {
  pdf_convert(file.path(res, paste0(f, ".pdf")), format = "png", pages = 1, dpi = 300, filenames = file.path(out, paste0(f, ".png")), verbose = FALSE)
}
print(list.files(out, pattern = "fig8|fig3|heatmap|venn|stackplot|ORA"))
