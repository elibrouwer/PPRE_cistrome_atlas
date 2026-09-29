## Threshold-free validation of PPARD.H13CORE.0.PSM.A against real ChIP-seq,
## the way motif-quality benchmarking papers do it (e.g. Weirauch et al. 2013,
## PMC4821295): instead of picking one p-value cutoff and asking "how many
## hits survive", score real ChIP peak sequences (positives) against matched
## background sequences (negatives) with the same PWM, then report a
## threshold-free discrimination metric -- AUC-ROC and AUPRC -- computed
## across every possible score cutoff at once.
##
## Ground truth: same 8 HUVEC PPARD ChIP-seq experiments used previously
## (recurrent peak set, >=2/8 experiments agree; see Objective1_fimo_chip_validation.R).
##
## IMPORTANT CORRECTION: an earlier version of this script scored regions using
## only the already-thresholded FIMO output (p<=1e-4). That left 99.5% of BOTH
## positive and negative regions tied at score=0 (most real ChIP peaks are too
## short/weak to contain a hit that strong), which makes ROC/AUC meaningless --
## nearly the entire ranking is an unbroken tie. Fixed here by rescanning every
## region's actual genomic DNA sequence directly against the real PWM (log-odds
## vs. uniform background), via Biostrings::PWMscoreStartingAt on both strands,
## with NO threshold at all -- a true continuous score for every region.

library(dplyr)
library(readr)
library(GenomicRanges)
library(GenomeInfoDb)
library(ggplot2)
library(scales)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)

base_dir <- "data/PPAR_RXR_motifs_Hocomoco_Jaspar"
hocomoco_meme_dir <- file.path(base_dir, "HOCOMOCO")
chip_bed_dir <- "data/pparD_huvec_chip"
motif_id <- "PPARD.H13CORE.0.PSM.A"
set.seed(1)

## ---- 0. Build a log-odds PWM directly from HOCOMOCO's probability matrix ----
parse_meme_pfm <- function(path) {
  lines <- readLines(path)
  header_idx <- grep("^letter-probability matrix", lines)
  w <- as.integer(sub(".*w=\\s*(\\d+).*", "\\1", lines[header_idx]))
  mat_lines <- lines[(header_idx + 1):(header_idx + w)]
  mat <- t(do.call(rbind, lapply(mat_lines, function(l) as.numeric(strsplit(trimws(l), "\\s+")[[1]]))))
  rownames(mat) <- c("A", "C", "G", "T")
  mat
}
prob_mat <- parse_meme_pfm(file.path(hocomoco_meme_dir, paste0(motif_id, "_meme_format.meme")))
prob_mat_smoothed <- pmax(prob_mat, 0.005)
prob_mat_smoothed <- sweep(prob_mat_smoothed, 2, colSums(prob_mat_smoothed), "/")
logodds_pwm <- log2(prob_mat_smoothed / 0.25)
motif_width <- ncol(logodds_pwm)
cat(sprintf("PWM for %s: width %d, max possible score %.2f\n",
            motif_id, motif_width, sum(apply(logodds_pwm, 2, max))))

score_region <- function(gr) {
  seqlengths(gr) <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)[seqlevels(gr)]
  seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38, gr)
  vapply(seqs, function(sq) {
    if (length(sq) < motif_width) return(NA_real_)
    starts <- seq_len(length(sq) - motif_width + 1)
    fwd <- max(PWMscoreStartingAt(logodds_pwm, sq, starting.at = starts))
    rev <- max(PWMscoreStartingAt(logodds_pwm, reverseComplement(sq), starting.at = starts))
    max(fwd, rev)
  }, numeric(1))
}

## ---- 1. Rebuild the real PPARD ChIP peak sets (union + recurrent) ----
chip_beds <- list.files(chip_bed_dir, pattern = "\\.bed$", full.names = TRUE)
stopifnot(length(chip_beds) == 8)
chip_peaks_list <- lapply(chip_beds, function(f) {
  df <- suppressMessages(read_tsv(f, col_names = c("chrom", "start", "end", "name", "score", "strand",
                                   "signalValue", "pValue", "qValue", "peak"),
                  col_types = cols(.default = "c")))
  GRanges(seqnames = df$chrom, ranges = IRanges(start = as.integer(df$start) + 1, end = as.integer(df$end)))
})
chip_gr_all <- suppressWarnings(do.call(c, unname(chip_peaks_list)))
chip_peaks_union <- suppressWarnings(reduce(chip_gr_all))
cov <- coverage(chip_gr_all)
chip_peaks_recurrent <- as(cov >= 2, "GRanges")
chip_peaks_recurrent <- chip_peaks_recurrent[chip_peaks_recurrent$score == TRUE]
seqlevelsStyle(chip_peaks_union) <- "UCSC"
seqlevelsStyle(chip_peaks_recurrent) <- "UCSC"

positives <- chip_peaks_recurrent
cat(sprintf("Positive set (real PPARD ChIP peaks, >=2/8 experiments): %d regions\n", length(positives)))

## ---- 2. Build matched negative/background set: same-width windows shifted
## away from every real peak (1:1 with positives), discarded if the shifted
## window itself lands on a real peak (prefer +5kb, fall back to -5kb, else drop) ----
offset <- 5000L
neg_plus  <- suppressWarnings(GenomicRanges::shift(positives, offset))
neg_minus <- suppressWarnings(GenomicRanges::shift(positives, -offset))
seqlevelsStyle(neg_plus) <- "UCSC"
seqlevelsStyle(neg_minus) <- "UCSC"

valid_plus  <- !suppressWarnings(overlapsAny(neg_plus, chip_peaks_union, minoverlap = 1))
valid_minus <- !suppressWarnings(overlapsAny(neg_minus, chip_peaks_union, minoverlap = 1))

negatives <- neg_plus
negatives[!valid_plus] <- neg_minus[!valid_plus]
keep <- valid_plus | valid_minus   # drop pairs where neither shifted window is clean

positives_matched <- positives[keep]
negatives_matched <- negatives[keep]

# keep only pairs on standard chromosomes BSgenome can actually serve
std_chroms <- paste0("chr", c(1:22, "X", "Y", "M"))
std_ok <- as.character(seqnames(positives_matched)) %in% std_chroms &
          as.character(seqnames(negatives_matched)) %in% std_chroms
positives_matched <- positives_matched[std_ok]
negatives_matched <- negatives_matched[std_ok]

cat(sprintf("Matched negative background windows (+/-%d bp from each peak, real peaks excluded, standard chromosomes only): %d\n",
            offset, length(negatives_matched)))

## ---- 3. Score every region's real DNA sequence directly against the PWM
## (continuous log-odds score, no FIMO / no threshold involved) ----
pos_scores <- score_region(positives_matched)
neg_scores <- score_region(negatives_matched)

# drop any region BSgenome couldn't fetch (off-scaffold, too short, etc.)
ok <- !is.na(pos_scores) & !is.na(neg_scores)
pos_scores <- pos_scores[ok]; neg_scores <- neg_scores[ok]

cat(sprintf("\nScored %d matched positive/negative pairs (dropped %d unscoreable)\n",
            length(pos_scores), sum(!ok)))
cat(sprintf("Positive region scores: median %.2f, range %.2f to %.2f\n",
            median(pos_scores), min(pos_scores), max(pos_scores)))
cat(sprintf("Negative region scores: median %.2f, range %.2f to %.2f\n",
            median(neg_scores), min(neg_scores), max(neg_scores)))

## ---- 4. Manual ROC / PR curve + AUC (no extra package dependency) ----
scored <- tibble(score = c(pos_scores, neg_scores),
                  label = c(rep(1L, length(pos_scores)), rep(0L, length(neg_scores))))

thresholds <- sort(unique(scored$score), decreasing = TRUE)
n_pos <- sum(scored$label == 1)
n_neg <- sum(scored$label == 0)

roc_pr <- lapply(thresholds, function(th) {
  tp <- sum(scored$score >= th & scored$label == 1)
  fp <- sum(scored$score >= th & scored$label == 0)
  fn <- n_pos - tp
  tibble(threshold = th, tpr = tp / n_pos, fpr = fp / n_neg,
         precision = tp / (tp + fp), recall = tp / n_pos)
}) %>% bind_rows() %>%
  arrange(fpr, tpr) %>%
  add_row(threshold = Inf, tpr = 0, fpr = 0, precision = 1, recall = 0, .before = 1) %>%
  add_row(threshold = -Inf, tpr = 1, fpr = 1, precision = n_pos / (n_pos + n_neg), recall = 1)

auc_roc <- sum(diff(roc_pr$fpr) * (head(roc_pr$tpr, -1) + tail(roc_pr$tpr, -1)) / 2)
pr_sorted <- roc_pr %>% arrange(recall)
auprc <- sum(diff(pr_sorted$recall) * (head(pr_sorted$precision, -1) + tail(pr_sorted$precision, -1)) / 2)
baseline_auprc <- n_pos / (n_pos + n_neg)

cat(sprintf("\nAUC-ROC = %.3f (0.5 = random)\nAUPRC = %.3f (baseline/random = %.3f)\n",
            auc_roc, auprc, baseline_auprc))

## ---- 5. Plot ROC + PR curves side by side ----
p_roc <- ggplot(roc_pr, aes(x = fpr, y = tpr)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
  geom_step(color = "#377EB8", linewidth = 1) +
  labs(x = "False positive rate", y = "True positive rate",
       title = sprintf("ROC (AUC = %.3f)", auc_roc)) +
  coord_equal() +
  theme_minimal(base_size = 11)

p_pr <- ggplot(pr_sorted, aes(x = recall, y = precision)) +
  geom_hline(yintercept = baseline_auprc, linetype = "dashed", color = "grey60") +
  geom_step(color = "#D73027", linewidth = 1) +
  labs(x = "Recall", y = "Precision",
       title = sprintf("Precision-Recall (AUPRC = %.3f, baseline = %.3f)", auprc, baseline_auprc)) +
  coord_equal() +
  theme_minimal(base_size = 11)

combined <- gridExtra::arrangeGrob(p_roc, p_pr, ncol = 2,
  top = grid::textGrob(paste0(motif_id, " -- threshold-free validation vs. real HUVEC PPARD ChIP-seq\n",
                               sprintf("%d positive (real peak) vs %d matched negative (flanking, +/-%dbp) regions, scored by direct PWM rescan",
                                       length(pos_scores), length(neg_scores), offset)),
                        gp = grid::gpar(fontsize = 12)))

out_path <- file.path(base_dir, "Objective1_results", paste0(motif_id, "_auc_validation.png"))
ggsave(out_path, combined, width = 11, height = 5.5, dpi = 150, bg = "white")
message("Saved: ", out_path)

write_csv(scored, file.path(base_dir, "Objective1_results", paste0(motif_id, "_auc_validation_scores.csv")))
