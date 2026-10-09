# Histogram of motif widths (whole base pairs) for the 42 motifs, stacked by family; replaces the Q-Q plot of width, which steps because widths are discrete.
suppressPackageStartupMessages({library(ggplot2)})
base <- "C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Minor internship"
out <- file.path(base, "PPAR_RXR_motifs_Hocomoco_Jaspar", "Human42_analysis", "Objective1_results", "Statistics_plots", "single_plots")
d <- read.csv(file.path(base, "R_code", "results", "PPRE_score_v4", "motif_normality_input.csv"), stringsAsFactors = FALSE)
d$width <- round(d$mean_width); d$Family <- factor(d$type_simple, levels = c("PPAR", "RXR", "PPAR::RXR"))
family_colors <- c("PPAR" = "#4F86D9", "RXR" = "#D96C6C", "PPAR::RXR" = "#9A86CC")
p <- ggplot(d, aes(x = width, fill = Family)) + geom_bar(width = 0.8) + scale_fill_manual(values = family_colors) +
  scale_x_continuous(breaks = seq(8, 20, 1)) + scale_y_continuous(breaks = seq(0, 10, 2), expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Motif width (bp)", y = "Number of motifs") + theme_classic() + theme(panel.grid = element_blank(), axis.line = element_line(colour = "black"), legend.title = element_blank(), legend.position = c(0.88, 0.8))
ggsave(file.path(out, "histogram_motif_width.pdf"), p, width = 4.6, height = 3.6); ggsave(file.path(out, "histogram_motif_width.png"), p, width = 4.6, height = 3.6, dpi = 300)
file.copy(file.path(out, "histogram_motif_width.png"), file.path(base, "R_code", "Claude outputs", "PPRE_boxplots_png", "histogram_motif_width.png"), overwrite = TRUE)
cat("saved", file.path(out, "histogram_motif_width.png"), "\n")
