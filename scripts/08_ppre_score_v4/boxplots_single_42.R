# Every boxplot of the 42-motif analysis as a separate file (PDF and PNG), so the panels can be assembled by hand.
suppressPackageStartupMessages({library(ggplot2); library(ggpubr)})
base <- "C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Minor internship"
root <- file.path(base, "PPAR_RXR_motifs_Hocomoco_Jaspar", "Human42_analysis", "Objective1_results")
out <- file.path(root, "Statistics_plots", "single_plots"); dir.create(out, recursive = TRUE, showWarnings = FALSE)
d <- read.csv(file.path(base, "R_code", "results", "PPRE_score_v4", "motif_normality_input.csv"), stringsAsFactors = FALSE)
rp <- read.csv(file.path(root, "hits_in_repeats_per_motif_42.csv"), stringsAsFactors = FALSE)
d <- merge(d, rp[, c("motif", "n_in_repeat", "pct_in_repeat")], by = "motif")
d$Motif_db <- factor(d$database, levels = c("HOCOMOCO", "JASPAR", "CIS-BP")); d$Motif_type_simple <- factor(d$type_simple, levels = c("PPAR", "RXR", "PPAR::RXR"))
d$Motif_type <- factor(ifelse(grepl("::", d$receptor), "PPAR::RXR", d$receptor), levels = c("PPARA", "PPARD", "PPARG", "RXRA", "RXRB", "RXRG", "PPAR::RXR"))
d$N_amount <- d$n_hits_1e4; d$Mean_width <- d$mean_width; d$log10_hits <- log10(d$n_hits_1e4)
pal <- list(Motif_db = c("HOCOMOCO" = "#6F6F6F", "JASPAR" = "#B5B5B5", "CIS-BP" = "#E4A11B"), Motif_type_simple = c("PPAR" = "#4F86D9", "RXR" = "#D96C6C", "PPAR::RXR" = "#9A86CC"),
            Motif_type = c("PPARA" = "#73A2E6", "PPARD" = "#8CB4EE", "PPARG" = "#B7D0F5", "RXRA" = "#E58C8C", "RXRB" = "#EDAAAA", "RXRG" = "#F4C9C9", "PPAR::RXR" = "#9A86CC"))
xlab <- c(Motif_db = "Motif database", Motif_type_simple = "Motif family", Motif_type = "Motif type"); xtag <- c(Motif_db = "database", Motif_type_simple = "family", Motif_type = "motif_type")
ylab <- c(N_amount = "Number of FIMO hits (p <= 1e-4)", Mean_width = "Motif width (bp)", n_in_repeat = "FIMO hits in repeats (n)", pct_in_repeat = "Hits in repeats (%)")
ytag <- c(N_amount = "hits_per_motif", Mean_width = "motif_width", n_in_repeat = "hits_in_repeats_n", pct_in_repeat = "hits_in_repeats_pct"); ylog <- c("N_amount", "n_in_repeat")
plot_theme <- theme_classic() + theme(panel.border = element_blank(), panel.grid = element_blank(), axis.line = element_line(colour = "black"), legend.title = element_blank(), legend.position = "none")
save2 <- function(p, name, w, h) { ggsave(file.path(out, paste0(name, ".pdf")), p, width = w, height = h); ggsave(file.path(out, paste0(name, ".png")), p, width = w, height = h, dpi = 300) }
for (x in names(pal)) for (y in names(ylab)) {
  p <- ggboxplot(d, x = x, y = y, fill = x, color = x, palette = pal[[x]], add = "jitter", width = 0.5, add.params = list(size = 2, alpha = 0.8)) +
    stat_compare_means(method = "kruskal.test", label.y = max(d[[y]], na.rm = TRUE) * 1.08) + labs(x = xlab[[x]], y = ylab[[y]]) + plot_theme
  save2(p, paste0("boxplot_", ytag[[y]], "_by_", xtag[[x]]), if (x == "Motif_type") 4.6 else 3.2, 3.6)
}
# Q-Q plots (normality), one file each
sw <- function(x) signif(shapiro.test(x)$p.value, 3)
for (v in c("Mean_width", "N_amount", "log10_hits")) {
  lab <- c(Mean_width = "Mean motif width", N_amount = "FIMO hits per motif", log10_hits = "log10 FIMO hits per motif")[[v]]
  q <- ggplot(d, aes(sample = .data[[v]])) + stat_qq() + stat_qq_line(color = "red") + labs(title = paste0("Q-Q plot: ", lab), subtitle = paste0("Shapiro-Wilk p = ", sw(d[[v]]), ", n = ", nrow(d)), x = "Theoretical quantiles", y = "Sample quantiles") + theme_classic()
  save2(q, paste0("qqplot_", c(Mean_width = "motif_width", N_amount = "hits_per_motif", log10_hits = "log10_hits_per_motif")[[v]]), 3.6, 3.6)
}
cat(length(list.files(out, pattern = "[.]pdf$")), "plots saved to", out, "\n"); print(sort(list.files(out, pattern = "[.]pdf$")))
