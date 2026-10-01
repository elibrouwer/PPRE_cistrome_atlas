# Composite PPRE-score figure (panels A-G) for the 42 human-source motifs.
# A: gene-set Venn, B: PPRE score median vs GTEx cardiomyocyte expression, C-F: per-gene score distributions (sum, max, median, average),
# G: GO-BP over-representation of the top 5000 genes by median score.
# Input: results/human42_pipeline/fimo_analysis/Bedfiles_output/PPRE_gtex_ranked_scores_abs_3k_test.bed (one row per gene, written by 03b).
# Per-hit score = FIMO score * exp(-|distance to TSS| / 3000)  (same decay constant as 03b; the older figure label said 2000).
suppressMessages({library(dplyr); library(ggplot2); library(patchwork); library(ggvenn); library(clusterProfiler); library(org.Hs.eg.db)})
root <- "C:/Users/brouw/PPRE_cistrome_atlas/results/human42_pipeline/fimo_analysis"
bed <- read.delim(file.path(root, "Bedfiles_output/PPRE_gtex_ranked_scores_abs_3k_test.bed"), check.names = FALSE)
# the bed holds one row per gene: gene-level sum/max/average as columns, and the median in the BED score column (see 03b: score = PPRE_score_median)
g <- bed %>% transmute(ensembl, symbol, PPRE_score_sum, PPRE_score_max, PPRE_score_median = score, PPRE_score_average, expression = GTEx_expression)
cat("genes with a scored PPRE:", nrow(g), "\n")
write.csv(g, file.path(root, "Scoring_FIMO/PPRE_gene_scores_human42.csv"), row.names = FALSE)

mark <- c("CPT1A", "HADHA", "HADHB", "ACADVL", "PPARA", "PPARG", "RXRA", "RXRB", "RXRG")
th <- theme_bw(base_size = 9) + theme(panel.border = element_blank(), panel.grid = element_blank(),
                                      axis.line = element_line(colour = "black"))
hist_panel <- function(col, title) {
  v <- g[g$symbol %in% mark, c("symbol", col)]; names(v)[2] <- "value"
  ggplot(g, aes(.data[[col]])) + geom_histogram(bins = 200, fill = "grey25") +
    geom_vline(data = v, aes(xintercept = value), colour = "blue", linetype = "dashed", linewidth = 0.3) +
    ggrepel::geom_text_repel(data = v, aes(x = value, y = Inf, label = symbol), inherit.aes = FALSE, angle = 90, direction = "x", vjust = 1, size = 2, colour = "blue", segment.size = 0.2) +
    scale_x_continuous(limits = c(0, NA), expand = c(0, 0)) + scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(title = title, x = NULL, y = "Frequency") + coord_cartesian(clip = "off") + th
}
# panel A: counts taken from venn_diagram.pdf of the same run
n <- c(g = 1298, ge = 227, gep = 15601, gp = 563, p = 3307)
ids <- function(k) paste0(k, seq_len(n[k]))
G <- c(ids("g"), ids("ge"), ids("gep"), ids("gp")); E <- c(ids("ge"), ids("gep")); P <- c(ids("gep"), ids("gp"), ids("p"))
pA <- ggvenn(list("GTEx genes" = G, "Expressed
cardiomyocyte genes" = E, "All motif
PPRE genes" = P), fill_alpha = 0, show_percentage = FALSE, stroke_size = 0.6, text_size = 2.6, set_name_size = 2.6)
pA <- pA + scale_x_continuous(expand = expansion(mult = 0.25)) + scale_y_continuous(expand = expansion(mult = 0.2))
pB <- ggplot(g, aes(expression, PPRE_score_median)) + geom_point(size = 0.3, alpha = 0.4) +
  geom_point(data = g[g$symbol %in% mark, ], colour = "blue", size = 1) +
  ggrepel::geom_text_repel(data = g[g$symbol %in% mark, ], aes(label = symbol), colour = "blue", size = 2.5) +
  labs(title = NULL, x = "Expression in log(TP10k+1)", y = "PPRE score median") + th

top <- g %>% arrange(desc(PPRE_score_median)) %>% slice_head(n = 5000) %>% pull(symbol)
map <- bitr(top, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
ego <- enrichGO(map$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 5000, readable = TRUE)
write.csv(as.data.frame(ego), file.path(root, "Scoring_FIMO/ORA_GOBP_median_top5000_human42.csv"), row.names = FALSE)
pG <- dotplot(ego, showCategory = 10, x = "GeneRatio", orderBy = "x", title = "PPRE score median") + theme(text = element_text(size = 7), axis.text.y = element_text(size = 6))

fig <- (pA | pB) / (hist_panel("PPRE_score_sum", "PPRE score sum") | hist_panel("PPRE_score_max", "PPRE score max")) /
  (hist_panel("PPRE_score_median", "PPRE score median") | hist_panel("PPRE_score_average", "PPRE score average")) / pG +
  plot_annotation(tag_levels = "A", caption = "Score = FIMO score x exp(-|distance to TSS| / 3000)") + plot_layout(heights = c(1, 1, 1, 1.3))
for (ext in c("png", "pdf")) ggsave(file.path(root, paste0("PPRE_score_figure_human42.", ext)), fig, width = 8, height = 12, dpi = 200)
cat("saved\n")
