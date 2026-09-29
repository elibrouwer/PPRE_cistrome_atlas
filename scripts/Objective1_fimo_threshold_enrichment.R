## Data-driven FIMO threshold selection: instead of picking a score/q-value
## cutoff by convention, test a range of thresholds against the real
## HCM/PLN differentially-accessible regions already used elsewhere in this
## project (Objective1_chipseeker_fimo_analysis_cleaned.Rmd, "overlap-regions"
## chunk) and see which threshold actually enriches for hits that fall inside
## real regulatory regions, versus hits expected by chance at that region
## density.

library(dplyr)
library(readr)
library(readxl)
library(GenomicRanges)
library(GenomeInfoDb)
library(rtracklayer)
library(ggplot2)
library(scales)

base_dir <- "data/PPAR_RXR_motifs_Hocomoco_Jaspar"
fimo_dir <- file.path(base_dir, "FIMO")
motif_id <- "PPARG.H13CORE.0.P.B"

## ---- 1. Load FIMO hits for this motif (same as the ladder-plot script) ----
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

hits_gr <- GRanges(
  seqnames = hits$sequence_name,
  ranges = IRanges(start = hits$start, end = hits$stop),
  strand = hits$strand,
  score = hits$score,
  pvalue = hits$`p-value`,
  qvalue = hits$`q-value`
)
seqlevelsStyle(hits_gr) <- "UCSC"

## ---- 2. Load HCM + PLN differential regions, lift to hg38 (same pipeline
## as the existing overlap-regions chunk) ----
HCM_differentially_regions <- read_xlsx(file.path(base_dir, "HCM_differentially_regions.xlsx"),
                                         sheet = "Supplementary Table 2A")
grange_HCM <- GRanges(seqnames = HCM_differentially_regions$`Chromosome number`,
                       ranges = IRanges(start = HCM_differentially_regions$`Region start`,
                                        end = HCM_differentially_regions$`Region end`))

PLN_differentially_regions <- read_xlsx(file.path(base_dir, "PLN_differentially_regions.xlsx"), skip = 1) %>%
  filter(!is.na(`Region start†`))
grange_PLN <- GRanges(seqnames = PLN_differentially_regions$Chromosome,
                       ranges = IRanges(start = PLN_differentially_regions$`Region start†`,
                                        end = PLN_differentially_regions$`Region end†`))

chain <- import.chain(file.path(base_dir, "hg19ToHg38.over.chain"))
grange_HCM_hg38 <- unlist(liftOver(grange_HCM, chain))
grange_PLN_hg38 <- unlist(liftOver(grange_PLN, chain))

reg_regions <- reduce(c(grange_HCM_hg38, grange_PLN_hg38))   # merge HCM+PLN, collapse overlaps
seqlevelsStyle(reg_regions) <- "UCSC"

## ---- 3. Genome-wide background: bp covered by regulatory regions vs whole genome ----
# Same genome FIMO was scanned against (Ensembl homo_sapiens_110 primary assembly, ~3.05 Gb non-N)
genome_bp <- 3.05e9
reg_bp <- sum(width(reg_regions))
expected_fraction <- reg_bp / genome_bp

cat(sprintf("Regulatory region coverage: %.0f bp across %d regions (%.4f%% of genome)\n",
            reg_bp, length(reg_regions), 100 * expected_fraction))

## ---- 4. For a range of p-value thresholds, compute observed enrichment ----
# (q-value was tried first but FIMO's q-values here only span ~0.016-0.31 --
# the FDR correction across the whole genome-wide multi-motif scan compresses
# them too tightly to give a useful graduated series. p-value (and score,
# which is monotonic with it for a fixed-width motif) spans many more orders
# of magnitude and matches the threshold table already used in the ladder plot.)
p_thresholds <- c(1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10)

enrichment_tbl <- lapply(p_thresholds, function(pt) {
  gr_sub <- hits_gr[hits_gr$pvalue <= pt]
  n_total <- length(gr_sub)
  if (n_total == 0) return(NULL)
  n_in_region <- length(unique(queryHits(findOverlaps(gr_sub, reg_regions, minoverlap = 1))))
  observed_fraction <- n_in_region / n_total
  fold_enrichment <- observed_fraction / expected_fraction
  # binomial test: is observed_fraction significantly above the by-chance expectation?
  p <- binom.test(n_in_region, n_total, p = expected_fraction, alternative = "greater")$p.value
  tibble(p_threshold = pt, n_total_hits = n_total, n_in_regulatory_region = n_in_region,
         observed_fraction = observed_fraction, fold_enrichment = fold_enrichment, binom_p = p)
}) %>% bind_rows()

print(enrichment_tbl, n = Inf)

## ---- 5. Plot enrichment vs threshold ----
p_enrich <- ggplot(enrichment_tbl, aes(x = factor(p_threshold, levels = rev(p_thresholds)), y = fold_enrichment)) +
  geom_col(fill = "#3B6FB6") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  geom_text(aes(label = comma(n_total_hits)), vjust = -0.4, size = 3) +
  labs(x = "p-value threshold (hits kept: p \u2264 threshold)",
       y = "fold enrichment in HCM/PLN regulatory regions\n(vs. genome-wide chance)",
       title = paste0(motif_id, " -- does a stricter FIMO threshold actually find more real regulatory sites?"),
       subtitle = "bars = fold enrichment over chance; labels = number of hits kept at that threshold; dashed line = no enrichment (fold = 1)") +
  theme_minimal(base_size = 11)

out_path <- file.path(base_dir, "Objective1_results", paste0(motif_id, "_threshold_enrichment.png"))
ggsave(out_path, p_enrich, width = 9, height = 6, dpi = 150, bg = "white")
message("Saved: ", out_path)

write_csv(enrichment_tbl, file.path(base_dir, "Objective1_results", paste0(motif_id, "_threshold_enrichment.csv")))
