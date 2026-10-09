# Boxplots of FIMO hits per motif and motif width for the 42-motif list (same style as Objective1_chipseeker_fimo_analysis_with_CISBP.Rmd).
# Rank-based tests are used because Shapiro-Wilk rejected normality (Wilcoxon for two groups, Kruskal-Wallis for more).
suppressPackageStartupMessages({library(ggplot2); library(ggpubr); library(patchwork)})
base <- "C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Minor internship"
out  <- file.path(base, "PPAR_RXR_motifs_Hocomoco_Jaspar", "Human42_analysis", "Objective1_results", "Statistics_plots")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
d <- read.csv(file.path(base, "R_code", "results", "PPRE_score_v4", "motif_normality_input.csv"), stringsAsFactors = FALSE)
d$Motif_db <- factor(d$database, levels = c("HOCOMOCO", "JASPAR", "CIS-BP"))
d$Motif_type_simple <- factor(d$type_simple, levels = c("PPAR", "RXR", "PPAR::RXR"))
d$Motif_type <- ifelse(grepl("::", d$receptor), "PPAR::RXR", d$receptor)
d$Motif_type <- factor(d$Motif_type, levels = c("PPARA", "PPARD", "PPARG", "RXRA", "RXRB", "RXRG", "PPAR::RXR"))
d$N_amount <- d$n_hits_1e4; d$Mean_width <- d$mean_width
database_colors <- c("HOCOMOCO" = "#6F6F6F", "JASPAR" = "#B5B5B5", "CIS-BP" = "#E4A11B")
family_colors <- c("PPAR" = "#4F86D9", "RXR" = "#D96C6C", "PPAR::RXR" = "#9A86CC")
type_colors <- c("PPARA" = "#73A2E6", "PPARD" = "#8CB4EE", "PPARG" = "#B7D0F5", "RXRA" = "#E58C8C", "RXRB" = "#EDAAAA", "RXRG" = "#F4C9C9", "PPAR::RXR" = "#9A86CC")
plot_theme <- theme_classic() + theme(panel.border = element_blank(), panel.grid = element_blank(), axis.line = element_line(colour = "black"), legend.title = element_blank(), legend.position = "none")
bp <- function(x, y, pal, xlab, ylab, ylog = FALSE) {
  p <- ggboxplot(d, x = x, y = y, fill = x, color = x, palette = pal, add = "jitter", width = 0.5, add.params = list(size = 2, alpha = 0.8)) +
    stat_compare_means(method = "kruskal.test", label.y = max(d[[y]], na.rm = TRUE) * 1.08) + labs(x = xlab, y = ylab) + plot_theme
  # linear axis, as in the original notebook (a log axis with the p-value label squashed the boxes)
  p
}
db_n <- bp("Motif_db", "N_amount", database_colors, "Motif database", "Number of FIMO hits (p <= 1e-4)", TRUE)
db_w <- bp("Motif_db", "Mean_width", database_colors, "Motif database", "Motif width (bp)")
fam_n <- bp("Motif_type_simple", "N_amount", family_colors, "Motif family", "Number of FIMO hits (p <= 1e-4)", TRUE)
fam_w <- bp("Motif_type_simple", "Mean_width", family_colors, "Motif family", "Motif width (bp)")
typ_n <- bp("Motif_type", "N_amount", type_colors, "Motif type", "Number of FIMO hits (p <= 1e-4)", TRUE)
typ_w <- bp("Motif_type", "Mean_width", type_colors, "Motif type", "Motif width (bp)")
ggsave(file.path(out, "even_boxplot_panel.pdf"), db_n + db_w + fam_n + fam_w + plot_layout(ncol = 4), width = 16, height = 4)
ggsave(file.path(out, "even_boxplot_panel_3.pdf"), typ_n + typ_w + plot_layout(ncol = 2), width = 16, height = 4)
ggsave(file.path(out, "FIMO_db_plts.pdf"), db_n + db_w + plot_layout(ncol = 2), width = 6, height = 3)
ggsave(file.path(out, "FIMO_mottyp_smp_plts.pdf"), fam_n + fam_w + plot_layout(ncol = 2), width = 6, height = 3)
sc <- ggplot(d, aes(x = Mean_width, y = N_amount, colour = Motif_type)) + geom_point(size = 3) +
  geom_smooth(aes(group = 1), method = "lm", se = FALSE, colour = "black") +                       # straight line of log10(hits) on width (the y axis is log10)
  stat_cor(aes(group = 1), method = "spearman", label.x.npc = "left", label.y.npc = "top", label.sep = ", ", colour = "black") +   # rank correlation: the data are not normal
  scale_colour_manual(values = type_colors) + scale_y_log10() +
  labs(x = "Motif width (bp)", y = "Number of FIMO hits (p <= 1e-4)") + plot_theme + theme(legend.position = "right")
ggsave(file.path(out, "FIMO_scat_plts.pdf"), sc, width = 8, height = 3, useDingbats = FALSE)
ggsave(file.path(out, "FIMO_scat_plts.png"), sc, width = 8, height = 3, dpi = 200)
cat("saved to", out, "\n"); print(list.files(out))
