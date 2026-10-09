# Repeat share of motif hits vs a GC-matched expectation: genome tiled in 500 bp windows, repeat share per window-GC bin (2% bins),
# expected share for a motif = mean over its hits of the repeat share of the hit's window-GC bin.
suppressPackageStartupMessages({library(GenomicRanges); library(data.table); library(BSgenome.Hsapiens.UCSC.hg38); library(ggplot2)})
base <- "C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Minor internship"
res_root <- file.path(base, "PPAR_RXR_motifs_Hocomoco_Jaspar", "Human42_analysis"); root <- file.path(res_root, "Objective1_results")
chrs <- paste0("chr", c(1:22, "X", "Y")); bs <- BSgenome.Hsapiens.UCSC.hg38; W <- 500L
rm_ <- fread("C:/Users/brouw/rmsk_hg38.bed", header = FALSE, select = 1:3); setnames(rm_, c("chr", "start", "end")); rm_ <- rm_[chr %in% chrs]
rep_gr <- reduce(GRanges(rm_$chr, IRanges(rm_$start + 1L, rm_$end)), ignore.strand = TRUE)
tabs <- list()
for (ch in chrs) {
  s <- bs[[ch]]; L <- length(s); st <- seq(1L, L - W + 1L, by = W); v <- Views(s, IRanges(st, width = W)); lf <- letterFrequency(v, c("G", "C", "N"), as.prob = FALSE)
  nn <- lf[, 3]; ok <- nn <= 0.1 * W; gc <- (lf[, 1] + lf[, 2]) / pmax(W - nn, 1)
  cv <- coverage(rep_gr[seqnames(rep_gr) == ch], width = L)[[ch]]; rb <- as.numeric(sum(Views(cv, IRanges(st, width = W))))
  tabs[[ch]] <- data.table(chr = ch, widx = seq_along(st), gcbin = round(gc / 0.02) * 0.02, rep_bp = rb, bp = W - nn, ok = ok)
}
T <- rbindlist(tabs); Tok <- T[ok == TRUE]; B <- Tok[, .(share = sum(rep_bp) / sum(bp), n = .N), by = gcbin][order(gcbin)]; print(B[n > 1000])
shareby <- setNames(B$share, as.character(B$gcbin))
files <- list.files(file.path(res_root, "hit_ranges_rds"), pattern = "[.]rds$", full.names = TRUE); rows <- list()
for (f in files) {
  m <- sub("[.]rds$", "", basename(f)); gr <- readRDS(f); gr <- gr[as.character(seqnames(gr)) %in% chrs]
  mid <- (start(gr) + end(gr)) %/% 2L; wi <- (mid - 1L) %/% W + 1L; key <- data.table(chr = as.character(seqnames(gr)), widx = wi)
  j <- merge(key[, .(chr, widx, i = .I)], T[, .(chr, widx, gcbin, ok)], by = c("chr", "widx"), all.x = TRUE); j <- j[!is.na(gcbin) & ok == TRUE]
  pexp <- shareby[as.character(j$gcbin)]; obs <- mean(overlapsAny(gr[j$i], rep_gr, ignore.strand = TRUE))
  rows[[m]] <- data.frame(motif = m, n_used = nrow(j), obs_pct = obs * 100, exp_gc_pct = mean(pexp, na.rm = TRUE) * 100, mean_hit_gc_window = mean(j$gcbin)); cat(sprintf("%-24s obs %.1f%% exp(GC) %.1f%%\n", m, obs * 100, mean(pexp, na.rm = TRUE) * 100))
}
R <- do.call(rbind, rows); R$fold_gc <- R$obs_pct / R$exp_gc_pct
M <- read.csv(file.path(root, "hits_in_repeats_vs_expected_42.csv")); R <- merge(M[, c("motif", "type", "receptor", "database", "fold")], R, by = "motif"); names(R)[names(R) == "fold"] <- "fold_genome"
write.csv(R, file.path(root, "hits_in_repeats_gc_matched_42.csv"), row.names = FALSE)
cat(sprintf("fold vs genome: mean %.2f (%.2f-%.2f) | fold vs GC-matched: mean %.2f (%.2f-%.2f) | Wilcoxon vs 1: p = %.3g\n", mean(R$fold_genome), min(R$fold_genome), max(R$fold_genome), mean(R$fold_gc), min(R$fold_gc), max(R$fold_gc), wilcox.test(R$fold_gc, mu = 1)$p.value))
print(aggregate(cbind(fold_genome, fold_gc, mean_hit_gc_window) ~ type, R, mean)); print(R[order(-R$fold_gc), c("motif", "type", "mean_hit_gc_window", "fold_genome", "fold_gc")][c(1:5, 38:42), ], digits = 3)
cat("Spearman fold_genome vs window GC:", cor(R$fold_genome, R$mean_hit_gc_window, method = "spearman"), "\n")
R$Motif_type_simple <- factor(ifelse(R$type == "heterodimer", "PPAR::RXR", R$type), levels = c("PPAR", "RXR", "PPAR::RXR"))
p <- ggplot(R, aes(fold_genome, fold_gc, colour = Motif_type_simple)) + geom_abline(linetype = "dashed") + geom_hline(yintercept = 1, colour = "grey60") + geom_point(size = 3) +
  scale_colour_manual(values = c("PPAR" = "#4F86D9", "RXR" = "#D96C6C", "PPAR::RXR" = "#9A86CC")) + theme_classic() + theme(legend.title = element_blank()) + labs(x = "Fold in repeats vs whole genome", y = "Fold in repeats vs GC-matched expectation")
sdir <- file.path(root, "Statistics_plots", "single_plots"); ggsave(file.path(sdir, "repeat_fold_genome_vs_gc_matched.png"), p, width = 4.6, height = 3.6, dpi = 300); ggsave(file.path(sdir, "repeat_fold_genome_vs_gc_matched.pdf"), p, width = 4.6, height = 3.6)
