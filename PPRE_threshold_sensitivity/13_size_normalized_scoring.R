# Fixes a confound discovered while validating the scores from 11: promoter
# window total size (union of -3000/+1000 across all UCSC knownGene isoforms)
# ranges from 4,000bp to 109,608bp across genes (27-fold), because genes with
# many widely-dispersed alternative TSSs (e.g. MACF1: 21 windows, 109,608bp)
# get proportionally more chances to accumulate motif hits by pure chance --
# nothing to do with real PPAR biology. This is exactly why the first
# RobustRankAggreg run (12) surfaced MACF1/TBCD/OBSCN/DMD-type genes (huge,
# complex loci) at the top instead of known PPAR targets: every raw-count and
# naive Fisher/Stouffer score in 11 was silently rewarding window size.
#
# Fix: the same size-normalized enrichment test this project already uses for
# exactly this kind of question (binomial/Poisson rate test against a
# genome-wide chance rate -- see Objective1_fimo_threshold_enrichment.R and
# 07_chip_validation.R), applied per gene instead of per threshold/family
# globally. For each gene and motif family: expected_hits = genome-wide rate
# (intervals per bp, from 07's n_total_intervals at p<=1e-4) x this gene's
# total promoter bp; one-sided Poisson test of observed vs. expected. A gene
# with a huge promoter is only rewarded if it has MORE hits than its own size
# alone would predict.

suppressPackageStartupMessages(library(data.table))
source("./PPRE_threshold_sensitivity/00_config.R")

genome_bp <- 3.05e9  # same constant as 07_chip_validation.R / Objective1_fimo_threshold_enrichment.R

## ---- per-gene total promoter width (bp), and per-gene/family observed loci counts ----
bed <- fread(promoter_windows_nonzero_bed, header = FALSE,
             col.names = c("chrom", "start0", "end", "gene", "score", "strand"))
gene_width <- bed[, .(promoter_bp = sum(end - start0), n_windows = .N), by = gene]

fisher_input <- fread(file.path(out_results, "gene_scoring_methods.csv"))[, .(
  gene, fisher_n_loci_combined, fisher_n_loci_heterodimer, fisher_n_loci_PPAR_half, fisher_n_loci_RXR_half
)]
gene_tbl <- merge(gene_width, fisher_input, by = "gene", all.x = TRUE)
for (col in c("fisher_n_loci_combined","fisher_n_loci_heterodimer","fisher_n_loci_PPAR_half","fisher_n_loci_RXR_half"))
  gene_tbl[is.na(get(col)), (col) := 0]

## ---- genome-wide per-bp rate, per family, at the recommended p<=1e-4 tier
## (n_total_intervals column already computed in 07_chip_validation.R) ----
chip_val <- fread(file.path(out_results, "chip_validation_PPARA_liver_by_threshold.csv"))
rate_tbl <- chip_val[threshold == baseline_threshold, .(view, n_total_intervals)]
rate <- setNames(rate_tbl$n_total_intervals / genome_bp, rate_tbl$view)
cat("Genome-wide per-bp interval rate (p<=1e-4), used as the Poisson null for every gene's promoter:\n")
print(rate)

poisson_enrich <- function(observed, width_bp, lambda_per_bp) {
  expected <- width_bp * lambda_per_bp
  # one-sided: is observed higher than expected given this gene's own promoter size?
  neglog10p <- -ppois(observed - 1L, expected, lower.tail = FALSE, log.p = TRUE) / log(10)
  fold <- observed / pmax(expected, 1e-9)
  list(neglog10p = neglog10p, expected = expected, fold = fold)
}

r <- poisson_enrich(gene_tbl$fisher_n_loci_combined,    gene_tbl$promoter_bp, rate["combined"]);    gene_tbl[, `:=`(poisson_neglog10p_combined = r$neglog10p, poisson_fold_combined = r$fold, expected_combined = r$expected)]
r <- poisson_enrich(gene_tbl$fisher_n_loci_heterodimer, gene_tbl$promoter_bp, rate["heterodimer"]); gene_tbl[, `:=`(poisson_neglog10p_heterodimer = r$neglog10p, poisson_fold_heterodimer = r$fold, expected_heterodimer = r$expected)]
r <- poisson_enrich(gene_tbl$fisher_n_loci_PPAR_half,   gene_tbl$promoter_bp, rate["PPAR_half"]);   gene_tbl[, `:=`(poisson_neglog10p_PPAR_half = r$neglog10p, poisson_fold_PPAR_half = r$fold, expected_PPAR_half = r$expected)]
r <- poisson_enrich(gene_tbl$fisher_n_loci_RXR_half,    gene_tbl$promoter_bp, rate["RXR_half"]);    gene_tbl[, `:=`(poisson_neglog10p_RXR_half = r$neglog10p, poisson_fold_RXR_half = r$fold, expected_RXR_half = r$expected)]

fwrite(gene_tbl, file.path(out_results, "gene_scoring_size_normalized.csv"))

## ---- sanity check: does the previously-anomalous size confound go away? ----
big_genes <- c("MACF1", "TBCD", "OBSCN", "DMD", "FOXP1")
known_genes <- c("CPT1A", "CPT1B", "HADHA", "HADHB", "ACADVL")
check <- gene_tbl[gene %in% c(big_genes, known_genes),
                   .(gene, promoter_bp, fisher_n_loci_combined, poisson_fold_combined, poisson_neglog10p_combined)]
cat("\n=== Size-confound check: large-promoter genes should no longer dominate once corrected ===\n")
print(check[order(-poisson_neglog10p_combined)])

cat(sprintf("\nCorrelation of promoter_bp with RAW loci count (combined): %.3f\n",
            cor(gene_tbl$promoter_bp, gene_tbl$fisher_n_loci_combined, method = "spearman")))
cat(sprintf("Correlation of promoter_bp with SIZE-NORMALIZED score (combined): %.3f (should be near 0)\n",
            cor(gene_tbl$promoter_bp, gene_tbl$poisson_neglog10p_combined, method = "spearman")))

message("\nDone. Wrote gene_scoring_size_normalized.csv")
