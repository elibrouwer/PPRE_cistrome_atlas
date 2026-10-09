# (1) Are motif hits enriched in repeats compared with the repeat share of the scanned genome? (2) Statistics table for the 12 boxplots (Kruskal-Wallis, epsilon-squared, BH, Dunn-Holm).
suppressPackageStartupMessages({library(GenomicRanges); library(data.table); library(BSgenome.Hsapiens.UCSC.hg38); library(rstatix)})
base <- "C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Minor internship"
root <- file.path(base, "PPAR_RXR_motifs_Hocomoco_Jaspar", "Human42_analysis", "Objective1_results")
chrs <- paste0("chr", c(1:22, "X", "Y")); bs <- BSgenome.Hsapiens.UCSC.hg38
# ---- expected repeat share: repeat bp / non-N bp over the scanned chromosomes
nonN <- GRanges(); for (ch in chrs) { s <- bs[[ch]]; nI <- as(masks(Biostrings::maskMotif(s, "N"))[[1]], "IRanges"); nonN <- c(nonN, GRanges(ch, gaps(nI, start = 1, end = length(s)))) }
nonN <- reduce(nonN)
rm_ <- fread("C:/Users/brouw/rmsk_hg38.bed", header = FALSE, select = 1:3); setnames(rm_, c("chr", "start", "end")); rm_ <- rm_[chr %in% chrs]
rep_gr <- reduce(GRanges(rm_$chr, IRanges(rm_$start + 1L, rm_$end)), ignore.strand = TRUE)
tot <- sum(width(nonN)); rep_bp <- sum(width(intersect(rep_gr, nonN, ignore.strand = TRUE)))
p0 <- rep_bp / tot; cat(sprintf("scanned (non-N) bp %.0f; repeat bp %.0f; expected repeat share %.2f%%\n", tot, rep_bp, p0 * 100))
# ---- per motif: observed share vs expected (exact binomial, BH across 42 motifs)
d <- read.csv(file.path(root, "hits_in_repeats_per_motif_42.csv"), stringsAsFactors = FALSE)
d$expected_pct <- p0 * 100; d$fold <- (d$pct_in_repeat / 100) / p0
d$p_binom <- mapply(function(k, n) binom.test(k, n, p0)$p.value, d$n_in_repeat, d$n_hits)
d$padj_BH <- p.adjust(d$p_binom, "BH")
d$direction <- ifelse(d$padj_BH < 0.05, ifelse(d$fold > 1, "enriched", "depleted"), "n.s.")
write.csv(d, file.path(root, "hits_in_repeats_vs_expected_42.csv"), row.names = FALSE)
w <- wilcox.test(d$pct_in_repeat, mu = p0 * 100)
cat(sprintf("motifs: mean %.1f%% (range %.1f-%.1f) vs expected %.1f%%; fold mean %.2f (range %.2f-%.2f); Wilcoxon signed-rank vs expected p = %.3g\n", mean(d$pct_in_repeat), min(d$pct_in_repeat), max(d$pct_in_repeat), p0 * 100, mean(d$fold), min(d$fold), max(d$fold), w$p.value))
print(table(d$direction)); print(aggregate(fold ~ type, d, function(x) round(mean(x), 3)))
# ---- plot: observed share per motif with the expected line
suppressPackageStartupMessages({library(ggplot2)})
d$Motif_type_simple <- factor(ifelse(d$type == "heterodimer", "PPAR::RXR", d$type), levels = c("PPAR", "RXR", "PPAR::RXR"))
d$motif_f <- factor(d$motif, levels = d$motif[order(d$pct_in_repeat)])
pl <- ggplot(d, aes(motif_f, pct_in_repeat, fill = Motif_type_simple)) + geom_col() + geom_hline(yintercept = p0 * 100, linetype = "dashed") +
  scale_fill_manual(values = c("PPAR" = "#4F86D9", "RXR" = "#D96C6C", "PPAR::RXR" = "#9A86CC")) + coord_flip() + theme_classic() + labs(x = NULL, y = "Hits in repeats (%)", fill = NULL, caption = sprintf("Dashed line: repeat share of the scanned genome (%.1f%%)", p0 * 100))
sdir <- file.path(root, "Statistics_plots", "single_plots"); ggsave(file.path(sdir, "repeat_share_vs_expected.png"), pl, width = 6, height = 8, dpi = 300); ggsave(file.path(sdir, "repeat_share_vs_expected.pdf"), pl, width = 6, height = 8)
# ---- statistics table for the 12 boxplots
m <- read.csv(file.path(base, "R_code", "results", "PPRE_score_v4", "motif_normality_input.csv"), stringsAsFactors = FALSE)
m <- merge(m, d[, c("motif", "n_in_repeat", "pct_in_repeat")], by = "motif")
m$Motif_db <- factor(m$database, levels = c("HOCOMOCO", "JASPAR", "CIS-BP")); m$Motif_family <- factor(m$type_simple, levels = c("PPAR", "RXR", "PPAR::RXR"))
m$Motif_type <- factor(ifelse(grepl("::", m$receptor), "PPAR::RXR", m$receptor), levels = c("PPARA", "PPARD", "PPARG", "RXRA", "RXRB", "RXRG", "PPAR::RXR"))
ys <- c(Hits = "n_hits_1e4", Width = "mean_width", Hits_in_repeats_n = "n_in_repeat", Hits_in_repeats_pct = "pct_in_repeat"); xs <- c("Motif_family", "Motif_db", "Motif_type")
rows <- list(); post <- list()
for (x in xs) for (yn in names(ys)) {
  df <- data.frame(y = m[[ys[[yn]]]], g = m[[x]]); k <- kruskal.test(y ~ g, df); n <- nrow(df); kk <- nlevels(droplevels(df$g))
  rows[[length(rows) + 1]] <- data.frame(grouping = x, measure = yn, n_groups = kk, min_group_n = min(table(droplevels(df$g))), H = unname(k$statistic), p = k$p.value, epsilon2 = unname((k$statistic - kk + 1) / (n - kk)), plot_role = ifelse(x == "Motif_type", "supplement", "main"))
  df$g <- droplevels(df$g)
}
S <- do.call(rbind, rows); S$p_BH <- p.adjust(S$p, "BH"); S$significant_BH <- S$p_BH < 0.05
for (i in which(S$significant_BH)) {
  df <- data.frame(y = m[[ys[[S$measure[i]]]]], g = droplevels(m[[S$grouping[i]]])); dt <- as.data.frame(dunn_test(df, y ~ g, p.adjust.method = "holm")); dt$grouping <- S$grouping[i]; dt$measure <- S$measure[i]; post[[length(post) + 1]] <- dt[, c("grouping", "measure", "group1", "group2", "n1", "n2", "statistic", "p", "p.adj")]
}
P <- do.call(rbind, post); write.csv(S, file.path(root, "Statistics_plots", "boxplot_statistics_42.csv"), row.names = FALSE); write.csv(P, file.path(root, "Statistics_plots", "boxplot_posthoc_dunn_holm_42.csv"), row.names = FALSE)
print(S, digits = 3); print(P, digits = 3)
