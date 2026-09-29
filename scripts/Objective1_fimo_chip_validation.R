## Validate FIMO threshold selection against REAL ChIP-seq peaks, instead of
## the HCM/PLN differential-accessibility proxy used previously.
##
## ChIP-Atlas lists 8 "PPARG" ChIP-seq experiments in Cardiovascular tissue
## (HUVEC), but inspecting each one individually shows they are all actually
## PPARD (PPARbeta/delta) ChIP -- mislabeled by ChIP-Atlas's automatic antigen
## curation (title/chip-antibody field says "PPARbeta/delta", curated
## "Antigen" field wrongly says "PPARG"). Since this project's FIMO output
## already includes PPARD.H13CORE.0.PSM.A, these 8 experiments are used as
## real, tissue-matched (cardiovascular) ground truth for THAT motif instead.
##
## GEO series SRA098847 (GSM1214678-81, GSM1226372-75): HUVEC, PPARbeta/delta
## ChIP, DMSO vs. agonist x normoxia vs. hypoxia x 2 replicates.

library(dplyr)
library(readr)
library(GenomicRanges)
library(GenomeInfoDb)
library(rtracklayer)
library(ggplot2)
library(scales)

base_dir <- "data/PPAR_RXR_motifs_Hocomoco_Jaspar"
fimo_dir <- file.path(base_dir, "FIMO")
chip_bed_dir <- "data/pparD_huvec_chip"
motif_id <- "PPARD.H13CORE.0.PSM.A"

## ---- 1. Load and merge the 8 real PPARD ChIP-seq peak sets (hg38, q<1e-5) ----
chip_beds <- list.files(chip_bed_dir, pattern = "\\.bed$", full.names = TRUE)
stopifnot(length(chip_beds) == 8)

chip_peaks_list <- lapply(chip_beds, function(f) {
  df <- read_tsv(f, col_names = c("chrom", "start", "end", "name", "score", "strand",
                                   "signalValue", "pValue", "qValue", "peak"),
                  col_types = cols(.default = "c"))
  GRanges(seqnames = df$chrom,
          ranges = IRanges(start = as.integer(df$start) + 1, end = as.integer(df$end)))
})
names(chip_peaks_list) <- basename(chip_beds)
cat("Peaks per experiment:\n")
print(sapply(chip_peaks_list, length))

chip_peaks_merged <- reduce(do.call(c, unname(chip_peaks_list)))
cat(sprintf("\nMerged PPARD ChIP peak set (union across 8 experiments): %d peaks, %s bp total\n",
            length(chip_peaks_merged), format(sum(width(chip_peaks_merged)), big.mark = ",")))

## also build a "reproducible" set: peaks supported by >=2 of 8 experiments,
## a stricter and more defensible ground truth than the raw union
chip_gr_all <- do.call(c, unname(chip_peaks_list))
cov <- coverage(chip_gr_all)
chip_peaks_recurrent <- as(cov >= 2, "GRanges")
chip_peaks_recurrent <- chip_peaks_recurrent[chip_peaks_recurrent$score == TRUE]
cat(sprintf("Recurrent PPARD ChIP peak set (>=2/8 experiments agree): %d peaks, %s bp total\n",
            length(chip_peaks_recurrent), format(sum(width(chip_peaks_recurrent)), big.mark = ",")))

## ---- 2. Load FIMO hits for PPARD.H13CORE.0.PSM.A ----
fimo_col_types <- cols(
  motif_id = col_character(), motif_alt_id = col_character(),
  sequence_name = col_character(), start = col_double(), stop = col_double(),
  strand = col_character(), score = col_double(), `p-value` = col_double(),
  `q-value` = col_double(), matched_sequence = col_character()
)
fimo_files <- list.files(fimo_dir, pattern = "\\.tsv$", full.names = TRUE)
hits <- lapply(fimo_files, function(f) {
  read_tsv(f, comment = "#", col_types = fimo_col_types) %>% filter(motif_id == !!motif_id)
}) %>% bind_rows()
stopifnot(nrow(hits) > 0)
cat(sprintf("\nTotal %s FIMO hits genome-wide: %s\n", motif_id, format(nrow(hits), big.mark = ",")))

hits_gr <- GRanges(seqnames = hits$sequence_name,
                    ranges = IRanges(start = hits$start, end = hits$stop),
                    strand = hits$strand, score = hits$score,
                    pvalue = hits$`p-value`, qvalue = hits$`q-value`)
seqlevelsStyle(hits_gr) <- "UCSC"
seqlevelsStyle(chip_peaks_merged) <- "UCSC"
seqlevelsStyle(chip_peaks_recurrent) <- "UCSC"

## ---- 3. Background: bp covered by ChIP peaks vs. whole genome ----
genome_bp <- 3.05e9

run_enrichment <- function(ground_truth_gr, label) {
  expected_fraction <- sum(width(ground_truth_gr)) / genome_bp
  p_thresholds <- c(1e-4, 1e-5, 1e-6, 1e-7, 1e-8)
  tbl <- lapply(p_thresholds, function(pt) {
    gr_sub <- hits_gr[hits_gr$pvalue <= pt]
    n_total <- length(gr_sub)
    if (n_total == 0) return(NULL)
    n_in <- length(unique(queryHits(findOverlaps(gr_sub, ground_truth_gr, minoverlap = 1))))
    obs_frac <- n_in / n_total
    fold <- obs_frac / expected_fraction
    p <- binom.test(n_in, n_total, p = expected_fraction, alternative = "greater")$p.value
    tibble(ground_truth = label, p_threshold = pt, n_total_hits = n_total,
           n_in_chip_peak = n_in, observed_fraction = obs_frac,
           fold_enrichment = fold, binom_p = p)
  }) %>% bind_rows()
  tbl
}

enrichment_union <- run_enrichment(chip_peaks_merged, "union (>=1/8 experiments)")
enrichment_recurrent <- run_enrichment(chip_peaks_recurrent, "recurrent (>=2/8 experiments)")
enrichment_tbl <- bind_rows(enrichment_union, enrichment_recurrent)

print(enrichment_tbl, n = Inf)

## ---- 4. Plot ----
p_enrich <- ggplot(enrichment_tbl, aes(x = factor(p_threshold, levels = rev(unique(p_threshold))),
                                        y = fold_enrichment, fill = ground_truth)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  geom_text(aes(label = comma(n_total_hits)), position = position_dodge(width = 0.8),
            vjust = -0.4, size = 2.8) +
  scale_fill_manual(values = c("union (>=1/8 experiments)" = "#8FAF3C", "recurrent (>=2/8 experiments)" = "#377EB8"),
                     name = "PPARD ChIP\nground truth") +
  labs(x = "p-value threshold (hits kept: p \u2264 threshold)",
       y = "fold enrichment in real PPARD ChIP peaks\n(HUVEC, vs. genome-wide chance)",
       title = paste0(motif_id, " -- FIMO hits validated against real HUVEC PPARD ChIP-seq"),
       subtitle = "8 GEO experiments (SRA098847), hg38, q<1e-5 peak calls; labels = FIMO hits kept at that threshold") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

out_path <- file.path(base_dir, "Objective1_results", paste0(motif_id, "_chip_validation.png"))
ggsave(out_path, p_enrich, width = 9.5, height = 6.5, dpi = 150, bg = "white")
message("Saved: ", out_path)

write_csv(enrichment_tbl, file.path(base_dir, "Objective1_results", paste0(motif_id, "_chip_validation.csv")))

## ---- 5. Same validation, but thresholding on q-value instead of p-value ----
# FIMO's q-value for this motif only spans ~0.079-0.252 (FDR correction across
# the whole genome-wide multi-motif scan compresses it) -- fixed round-number
# cutoffs like q<=0.05 or q<=0.01 return ZERO hits immediately, since nothing
# in the data reaches that low. Quantile-based thresholds (top X% of hits
# ranked by q-value) give a usable graduated series from this same
# distribution instead of an unusable fixed scale.
cat(sprintf("\nq-value range for %s: %.4f to %.4f\n",
            motif_id, min(hits_gr$qvalue), max(hits_gr$qvalue)))

q_fracs <- c(0.01, 0.05, 0.10, 0.25, 0.50, 1.00)   # top X% of hits by q-value (ascending = most confident)
q_thresholds_val <- quantile(hits_gr$qvalue, probs = q_fracs)

run_enrichment_q <- function(ground_truth_gr, label) {
  expected_fraction <- sum(width(ground_truth_gr)) / genome_bp
  tbl <- lapply(seq_along(q_fracs), function(i) {
    qt <- q_thresholds_val[i]
    gr_sub <- hits_gr[hits_gr$qvalue <= qt]
    n_total <- length(gr_sub)
    if (n_total == 0) return(NULL)
    n_in <- length(unique(queryHits(findOverlaps(gr_sub, ground_truth_gr, minoverlap = 1))))
    obs_frac <- n_in / n_total
    fold <- obs_frac / expected_fraction
    p <- binom.test(n_in, n_total, p = expected_fraction, alternative = "greater")$p.value
    tibble(ground_truth = label, top_fraction = q_fracs[i], q_threshold = qt,
           n_total_hits = n_total, n_in_chip_peak = n_in,
           observed_fraction = obs_frac, fold_enrichment = fold, binom_p = p)
  }) %>% bind_rows()
  tbl
}

enrichment_union_q <- run_enrichment_q(chip_peaks_merged, "union (>=1/8 experiments)")
enrichment_recurrent_q <- run_enrichment_q(chip_peaks_recurrent, "recurrent (>=2/8 experiments)")
enrichment_tbl_q <- bind_rows(enrichment_union_q, enrichment_recurrent_q)

print(enrichment_tbl_q, n = Inf)

q_label <- function(frac, qt) sprintf("top %g%%\n(q≤%.3f)", frac * 100, qt)
enrichment_tbl_q <- enrichment_tbl_q %>%
  mutate(x_label = q_label(top_fraction, q_threshold))
x_order <- enrichment_tbl_q %>% distinct(top_fraction, x_label) %>% arrange(top_fraction) %>% pull(x_label)

p_enrich_q <- ggplot(enrichment_tbl_q, aes(x = factor(x_label, levels = x_order),
                                            y = fold_enrichment, fill = ground_truth)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  geom_text(aes(label = comma(n_total_hits)), position = position_dodge(width = 0.8),
            vjust = -0.4, size = 2.8) +
  scale_fill_manual(values = c("union (>=1/8 experiments)" = "#8FAF3C", "recurrent (>=2/8 experiments)" = "#377EB8"),
                     name = "PPARD ChIP\nground truth") +
  labs(x = "hits kept (top X% by q-value, most confident first)",
       y = "fold enrichment in real PPARD ChIP peaks\n(HUVEC, vs. genome-wide chance)",
       title = paste0(motif_id, " -- same validation, thresholded by q-value quantile"),
       subtitle = sprintf("q-value only spans %.3f-%.3f for this motif, so fixed cutoffs (q≤0.05 etc.) are unusable -- quantiles used instead; labels = FIMO hits kept",
                           min(hits_gr$qvalue), max(hits_gr$qvalue))) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

out_path_q <- file.path(base_dir, "Objective1_results", paste0(motif_id, "_chip_validation_qvalue.png"))
ggsave(out_path_q, p_enrich_q, width = 9.5, height = 6.5, dpi = 150, bg = "white")
message("Saved: ", out_path_q)

write_csv(enrichment_tbl_q, file.path(base_dir, "Objective1_results", paste0(motif_id, "_chip_validation_qvalue.csv")))
