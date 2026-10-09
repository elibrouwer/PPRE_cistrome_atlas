# Pipeline for the 42-motif list, part 3: pairwise motif overlap (percent of hits, minoverlap 5) and base-pair Jaccard heatmaps (notebook chunks overlap-heatmap and motif-jaccard-heatmap).
suppressPackageStartupMessages({library(GenomicRanges); library(pheatmap); library(dplyr)})
base <- "C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Minor internship"
res_root <- file.path(base, "PPAR_RXR_motifs_Hocomoco_Jaspar", "Human42_analysis"); rng_dir <- file.path(res_root, "hit_ranges_rds"); res_dir <- file.path(res_root, "Objective1_results")
M <- read.csv(file.path(base, "R_code", "results", "PPRE_score_v4", "motif_summary_42.csv"))
files <- list.files(rng_dir, pattern = "\\.rds$", full.names = TRUE); nm <- sub("\\.rds$", "", basename(files)); gl <- lapply(files, readRDS); names(gl) <- nm
n <- length(gl); cat("motifs", n, "\n"); t0 <- Sys.time()
red <- lapply(gl, reduce)
ov <- matrix(100, n, n, dimnames = list(nm, nm)); jac <- matrix(100, n, n, dimnames = list(nm, nm))
for (i in seq_len(n)) {
  for (j in seq_len(n)) if (i != j) ov[i, j] <- mean(overlapsAny(gl[[i]], gl[[j]], minoverlap = 5)) * 100
  for (j in seq_len(n)) if (j > i) { inter <- sum(width(intersect(red[[i]], red[[j]]))); uni <- sum(width(reduce(c(red[[i]], red[[j]])))); jac[i, j] <- jac[j, i] <- inter / uni * 100 }
  cat(sprintf("%s done (%d/%d) %.1f min\n", nm[i], i, n, as.numeric(difftime(Sys.time(), t0, units = "mins")))); flush.console()
}
write.csv(round(ov, 2), file.path(res_dir, "motif_overlap_percent_42.csv")); write.csv(round(jac, 2), file.path(res_dir, "motif_jaccard_percent_42.csv"))
ann <- data.frame(Family = M$type, Database = M$database, row.names = M$motif); ann <- ann[nm, , drop = FALSE]
pdf(file.path(res_dir, "overlap_heatmap_42.pdf"), width = 12.3, height = 11.7); pheatmap(round(ov, 1), cluster_rows = TRUE, cluster_cols = TRUE, annotation_row = ann, annotation_col = ann, fontsize = 7); dev.off()
pdf(file.path(res_dir, "jaccard_heatmap_42.pdf"), width = 12.3, height = 11.7); pheatmap(round(jac, 1), cluster_rows = TRUE, cluster_cols = TRUE, annotation_row = ann, annotation_col = ann, fontsize = 7); dev.off()
off <- ov[upper.tri(ov) | lower.tri(ov)]; offj <- jac[upper.tri(jac)]
cat(sprintf("overlap percent (off-diagonal): mean %.1f, median %.1f, max %.1f | Jaccard: mean %.1f, median %.1f, max %.1f\n", mean(off), median(off), max(off), mean(offj), median(offj), max(offj)))
