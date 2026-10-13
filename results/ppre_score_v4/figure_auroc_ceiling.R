# Figure: AUROC of the site-count score against empirical ceilings (ChIP-binding oracle, label reproducibility between the two hepatocyte datasets).
suppressPackageStartupMessages(library(ggplot2))
v4 <- "C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Minor internship/R_code/results/PPRE_score_v4"
out <- "C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Minor internship/R_code/Claude outputs/PPRE_boxplots_png"
x <- read.delim(file.path(v4, "v4ceiling_results.tsv"), stringsAsFactors = FALSE)
x$top <- ifelse(grepl("top 1%", x$test), "Top 1% responders", "Top 5% responders"); x$target <- sub(".*(GSE\\d+) top.*", "\\1", x$test)
x$group <- ifelse(grepl("^ChIP oracle", x$test), "ChIP binding in promoter\n(oracle)", "Other hepatocyte dataset\n(label reproducibility)")
x$predictor <- ifelse(grepl("^ChIP oracle", x$test), sub("ChIP oracle (\\S+) ->.*", "\\1", x$test), "other dataset")
sc <- read.csv(file.path(v4, "v4repeat_enriched_results.tsv"), sep = "\t", check.names = FALSE)
s <- data.frame(test = "score", AUROC = unlist(sc[1, c("hep GSE17251 top 5%", "hep GSE17251 top 1%", "hep GSE53399 top 5%", "hep GSE53399 top 1%")]), n_genes = NA,
                top = rep(c("Top 5% responders", "Top 1% responders"), 2), target = rep(c("GSE17251", "GSE53399"), each = 2), group = "PPRE site-count score\n(promoter sequence only)", predictor = "score")
d <- rbind(s, x); d$group <- factor(d$group, levels = rev(c("PPRE site-count score\n(promoter sequence only)", "ChIP binding in promoter\n(oracle)", "Other hepatocyte dataset\n(label reproducibility)")))
d$top <- factor(d$top, levels = c("Top 5% responders", "Top 1% responders"))
p <- ggplot(d, aes(AUROC, group)) + geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey50") +
  geom_point(aes(shape = target, colour = group), position = position_jitter(height = 0.18, width = 0, seed = 1), size = 2.6, alpha = 0.85) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, linewidth = 0.4, fatten = 0, aes(xmin = after_stat(x), xmax = after_stat(x)), colour = "black") +
  facet_wrap(~top, ncol = 1) + scale_colour_manual(values = c("PPRE site-count score\n(promoter sequence only)" = "#4F86D9", "ChIP binding in promoter\n(oracle)" = "#D96C6C", "Other hepatocyte dataset\n(label reproducibility)" = "#7A7A7A"), guide = "none") +
  scale_x_continuous(limits = c(0.45, 1), breaks = seq(0.5, 1, 0.1)) + labs(x = "AUROC (predicting hepatocyte PPARα-agonist responders)", y = NULL, shape = "Labels from") +
  theme_classic() + theme(strip.background = element_blank(), strip.text = element_text(face = "bold", hjust = 0), legend.position = "bottom")
ggsave(file.path(out, "figure_auroc_ceiling.png"), p, width = 5.6, height = 4.6, dpi = 300); ggsave(file.path(v4, "figure_auroc_ceiling.pdf"), p, width = 5.6, height = 4.6)
