# Hits in repeat elements per motif (42-motif list): counts, boxplots by database / family / motif type (same style as boxplots_42.R).
suppressPackageStartupMessages({library(GenomicRanges); library(ggplot2); library(ggpubr); library(patchwork); library(data.table)})
base <- "C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Minor internship"
res_root <- file.path(base, "PPAR_RXR_motifs_Hocomoco_Jaspar", "Human42_analysis"); rng_dir <- file.path(res_root, "hit_ranges_rds")
out <- file.path(res_root, "Objective1_results", "Statistics_plots"); dir.create(out, showWarnings = FALSE, recursive = TRUE)
# repeat mask: UCSC RepeatMasker hg38 (bed already on disk), merged into one mask
rm_ <- fread("C:/Users/brouw/rmsk_hg38.bed", header = FALSE, select = 1:3); setnames(rm_, c("chr", "start", "end"))
rm_ <- rm_[chr %in% paste0("chr", c(1:22, "X", "Y"))]; rep_gr <- reduce(GRanges(rm_$chr, IRanges(rm_$start + 1L, rm_$end)), ignore.strand = TRUE)
cat("repeat intervals (merged):", length(rep_gr), "\n")
files <- list.files(rng_dir, pattern = "[.]rds$", full.names = TRUE); rows <- list()
for (f in files) {
  m <- sub("[.]rds$", "", basename(f)); gr <- readRDS(f); inrep <- overlapsAny(gr, rep_gr, ignore.strand = TRUE)
  rows[[m]] <- data.frame(motif = m, n_hits = length(gr), n_in_repeat = sum(inrep), pct_in_repeat = mean(inrep) * 100, n_not_in_repeat = sum(!inrep)); cat(sprintf("%-24s %9d hits, %5.1f%% in repeats\n", m, length(gr), mean(inrep) * 100))
}
R <- do.call(rbind, rows); M <- read.csv(file.path(base, "R_code", "results", "PPRE_score_v4", "motif_summary_42.csv"))
d <- merge(R, M[, c("motif", "receptor", "type", "database")], by = "motif"); write.csv(d, file.path(res_root, "Objective1_results", "hits_in_repeats_per_motif_42.csv"), row.names = FALSE)
cat(sprintf("all motifs: %.1f%% of hits in repeats (mean per motif %.1f%%, SD %.1f, range %.1f-%.1f)\n", sum(d$n_in_repeat) / sum(d$n_hits) * 100, mean(d$pct_in_repeat), sd(d$pct_in_repeat), min(d$pct_in_repeat), max(d$pct_in_repeat)))
d$Motif_db <- factor(d$database, levels = c("HOCOMOCO", "JASPAR", "CIS-BP")); d$Motif_type_simple <- factor(ifelse(d$type == "heterodimer", "PPAR::RXR", d$type), levels = c("PPAR", "RXR", "PPAR::RXR"))
d$Motif_type <- factor(ifelse(grepl("::", d$receptor), "PPAR::RXR", d$receptor), levels = c("PPARA", "PPARD", "PPARG", "RXRA", "RXRB", "RXRG", "PPAR::RXR"))
database_colors <- c("HOCOMOCO" = "#6F6F6F", "JASPAR" = "#B5B5B5", "CIS-BP" = "#E4A11B"); family_colors <- c("PPAR" = "#4F86D9", "RXR" = "#D96C6C", "PPAR::RXR" = "#9A86CC")
type_colors <- c("PPARA" = "#73A2E6", "PPARD" = "#8CB4EE", "PPARG" = "#B7D0F5", "RXRA" = "#E58C8C", "RXRB" = "#EDAAAA", "RXRG" = "#F4C9C9", "PPAR::RXR" = "#9A86CC")
plot_theme <- theme_classic() + theme(panel.border = element_blank(), panel.grid = element_blank(), axis.line = element_line(colour = "black"), legend.title = element_blank(), legend.position = "none")
bp <- function(x, y, pal, xlab, ylab, ylog = FALSE) { p <- ggboxplot(d, x = x, y = y, fill = x, color = x, palette = pal, add = "jitter", width = 0.5, add.params = list(size = 2, alpha = 0.8)) +
  stat_compare_means(method = "kruskal.test", label.y = max(d[[y]], na.rm = TRUE) * 1.08) + labs(x = xlab, y = ylab) + plot_theme; p }
a <- bp("Motif_db", "n_in_repeat", database_colors, "Motif database", "FIMO hits in repeats (n)", TRUE); b <- bp("Motif_db", "pct_in_repeat", database_colors, "Motif database", "Hits in repeats (%)")
c_ <- bp("Motif_type_simple", "n_in_repeat", family_colors, "Motif family", "FIMO hits in repeats (n)", TRUE); e <- bp("Motif_type_simple", "pct_in_repeat", family_colors, "Motif family", "Hits in repeats (%)")
f1 <- bp("Motif_type", "n_in_repeat", type_colors, "Motif type", "FIMO hits in repeats (n)", TRUE); g1 <- bp("Motif_type", "pct_in_repeat", type_colors, "Motif type", "Hits in repeats (%)")
ggsave(file.path(out, "repeat_boxplot_panel.pdf"), a + b + c_ + e + plot_layout(ncol = 4), width = 16, height = 4)
ggsave(file.path(out, "repeat_boxplot_panel_motif_type.pdf"), f1 + g1 + plot_layout(ncol = 2), width = 16, height = 4)
# per-motif bar chart of the repeat share, sorted
d$motif_f <- factor(d$motif, levels = d$motif[order(d$pct_in_repeat)])
bar <- ggplot(d, aes(motif_f, pct_in_repeat, fill = Motif_db)) + geom_col() + scale_fill_manual(values = database_colors) + coord_flip() + theme_classic() + labs(x = NULL, y = "Hits in repeats (%)", fill = NULL)
ggsave(file.path(out, "repeat_share_per_motif.pdf"), bar, width = 6, height = 8)
cat("Kruskal-Wallis, share in repeats: by database p =", signif(kruskal.test(pct_in_repeat ~ Motif_db, d)$p.value, 3), "| by family p =", signif(kruskal.test(pct_in_repeat ~ Motif_type_simple, d)$p.value, 3), "\n")
cat("saved:", paste(list.files(out, pattern = "repeat"), collapse = ", "), "\n")
