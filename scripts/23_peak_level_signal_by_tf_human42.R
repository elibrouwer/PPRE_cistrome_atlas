# Peak-level motif signal for PPARG, PPARA and RXRA ChIP-seq (ChIP-Atlas), per cell group: do FIMO hits of the motif subsets separate peak windows from GC-matched controls?
# (derived from 19_ppar_chip_calibration_human42.R, 10,000 peaks max per set, no model fitting)
# Windows: 500 bp centred on peaks (MACS2 score >= 250, <= 20,000 peaks per group) vs 3 GC-matched control windows per peak (1 kb GC class of 2.5 %, N-free, not overlapping a peak).
# Features: per-motif maximum relative FIMO score (score / motif max) in the window, from the unmasked p <= 1e-4 hits of the 42 motifs.
# Models: covariate only (1 kb GC), all-motif max, heterodimer-motif max, elastic net on the 42 per-motif features (+ GC). Chromosome-blocked evaluation (train odd / test even chromosomes and vice versa) in adipocyte, then transfer to the other groups.
BENCH <- "C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Minor internship/R_code/benchmark"
suppressPackageStartupMessages({ library(data.table); library(GenomicRanges); library(glmnet) })
set.seed(20251001)
cache_dir <- file.path(BENCH, "cache"); STD <- paste0("chr", c(1:22, "X"))
out <- "C:/Users/brouw/PPRE_cistrome_atlas/results/human42_pipeline/chip_calibration"; dir.create(out, showWarnings = FALSE)
fimo_root <- "C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Minor internship/PPAR_RXR_motifs_Hocomoco_Jaspar/FIMO_results_human42_unmasked"
inv <- fread("C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Minor internship/PPAR_RXR_motifs_Hocomoco_Jaspar/Motifs_update_2026-09-30/motif_inventory_gates.csv", encoding = "UTF-8")
het <- inv[category == "PPAR:RXR heterodimer", motif_id]
gc <- fread("C:/Users/brouw/PPRE_cistrome_atlas/results/human42_pipeline/clustering/gc_bins_1kb.tsv"); gc[, `:=`(gcf = gc / (length - n), cls = floor(100 * gc / (length - n) / 2.5))]
gcok <- gc[n == 0 & length == 1000]

pk <- rbindlist(lapply(c("PPARG", "PPARA", "RXRA"), function(tf) { d <- fread(file.path(cache_dir, paste0("chipatlas_", tf, "_min.tsv")), header = FALSE); setnames(d, c("chr", "start", "end", "srx", "group", "cell", "score")); d[, tf := tf]; d }))
pk <- pk[chr %in% STD & score >= 250]; pk[, group := paste(tf, gsub("%20", " ", group))]
pk <- pk[group %in% pk[, .N, by = group][N >= 300, group]]
cat("peak sets (score>=250, >=300 peaks):
"); print(pk[, .N, by = group])
W <- list()
for (g in unique(pk$group)) {
  p <- pk[group == g]; p <- p[sample(.N, min(.N, 10000))]
  p[, mid := (start + end) %/% 2L]
  gr_p <- GRanges(pk[group == g]$chr, IRanges(pk[group == g]$start + 1L, pk[group == g]$end))
  pw <- p[, .(chr, pos = mid, y = 1L)]
  pw[, bin := pos %/% 1000L]; pw <- merge(pw, gc[, .(chr = chrom, bin, gcf, cls)], by = c("chr", "bin"), all.x = TRUE); pw <- pw[!is.na(cls)]
  cw <- rbindlist(lapply(split(pw, pw$cls), function(d) { pool <- gcok[cls == d$cls[1]]; if (!nrow(pool)) return(NULL)
    s <- pool[sample(.N, 3 * nrow(d), replace = TRUE)]; data.table(chr = s$chrom, pos = s$bin * 1000L + sample(250:749, nrow(s), TRUE), y = 0L, gcf = s$gcf, cls = s$cls) }))
  ok <- !overlapsAny(GRanges(cw$chr, IRanges(cw$pos - 250L, cw$pos + 250L)), gr_p); cw <- cw[ok]
  w <- rbind(pw[, .(chr, pos, y, gcf)], cw[, .(chr, pos, y, gcf)]); w[, group := g]; W[[g]] <- w
}
W <- rbindlist(W); W[, id := .I]
gw <- GRanges(W$chr, IRanges(W$pos - 250L, W$pos + 250L))

# per-motif max relative score per window
feat <- list()
for (f in list.files(fimo_root, "fimo\\.tsv$", recursive = TRUE, full.names = TRUE)) {
  d <- fread(f, sep = "\t", select = c("motif_id", "sequence_name", "start", "stop", "score"), colClasses = list(character = c("motif_id", "sequence_name")))
  d <- d[sequence_name %in% as.character(c(1:22, "X")) & !startsWith(motif_id, "#")]
  d[, `:=`(start = as.integer(start), stop = as.integer(stop), score = as.numeric(score))]
  d[, rel := score / max(score), by = motif_id]
  ov <- findOverlaps(GRanges(paste0("chr", d$sequence_name), IRanges(d$start, d$stop)), gw)
  a <- data.table(id = subjectHits(ov), motif = d$motif_id[queryHits(ov)], rel = d$rel[queryHits(ov)])[, .(rel = max(rel)), by = .(id, motif)]
  feat[[f]] <- a; cat(basename(dirname(f)), "done\n")
}
F <- rbindlist(feat); motifs <- sort(unique(inv$motif_id[inv$motif_id %in% F$motif]))
M <- matrix(0, nrow(W), length(motifs), dimnames = list(NULL, motifs)); M[cbind(F$id, match(F$motif, motifs))] <- F$rel
auc <- function(x, y) { r <- rank(x); n1 <- sum(y == 1); n0 <- sum(y == 0); (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
cg <- setNames(inv$category, inv$motif_id)[motifs]
sub <- list(`all 42` = rep(TRUE, length(motifs)), `PPAR-containing 14` = cg %in% c("PPAR", "PPAR:RXR heterodimer"), `heterodimer 3` = cg == "PPAR:RXR heterodimer", `RXR-only 28` = cg == "RXR")
R <- rbindlist(lapply(unique(W$group), function(g) { i <- which(W$group == g); y <- W$y[i]
  data.table(peak_set = g, n_peaks = sum(y == 1), GC = auc(W$gcf[i], y), t(sapply(sub, function(s) auc(apply(M[i, s, drop = FALSE], 1, max), y))),
             pct_peaks_het_hit_rel0.7 = 100 * mean(apply(M[i[y == 1], cg == "PPAR:RXR heterodimer", drop = FALSE], 1, max) >= 0.7),
             pct_controls_het_hit_rel0.7 = 100 * mean(apply(M[i[y == 0], cg == "PPAR:RXR heterodimer", drop = FALSE], 1, max) >= 0.7)) }))
fwrite(R, file.path(out, "peak_level_motif_signal_by_tf_and_group.csv"))
print(R[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)], nrows = 50)
