# Shared paths and settings for the PPRE threshold-sensitivity analysis.
# Sourced by every other script in this folder -- do not run directly.
#
# Scope decision (confirmed with the project owner, 2026-09-25): this analysis
# covers p<=1e-4, 1e-5, 1e-6 only. The 5e-4 and 1e-3 tiers requested in the
# original brief are NOT computable from existing data -- every complete FIMO
# run in this project (FIMO_results_no_max, FIMO_results_genome_masked,
# FIMO_results_promoter_ownbg) was itself run with FIMO's own --thresh 1e-4,
# so hits weaker than p=1e-4 were never scored or stored by FIMO and cannot be
# recovered after the fact. See README.md for the exact commands needed to
# extend this analysis to 5e-4/1e-3 if a rerun becomes available.

base_dir   <- "data/PPAR_RXR_motifs_Hocomoco_Jaspar"
repo_dir   <- "."
gtex_dir   <- "data/GTEX_single_nuclear_data"

# FIMO_results_no_max is the run used here: unmasked GRCh38 primary assembly,
# genome-wide order-0 background, --max-stored-scores 20000000 (no truncation),
# all 30 motifs across 8 shards. This is the same run PPARD_threshold_validation
# identified as authoritative after comparing all four available runs -- see
# the audit notes in README.md for why the original FIMO/ directory (truncated,
# 22 motifs only) and the masked/promoter-own-background runs are NOT used here.
fimo_dir   <- file.path(base_dir, "FIMO_results_no_max")
fimo_shards <- c("1_4", "5_8", "9_12", "13_16", "17_20", "21_22", "23_26", "27_30")

out_dir     <- file.path(base_dir, "PPRE_threshold_sensitivity")
out_data    <- file.path(out_dir, "data")
out_results <- file.path(out_dir, "results")
out_figures <- file.path(out_dir, "figures")
for (d in c(out_data, out_results, out_figures)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

motif_family_csv <- file.path(out_results, "motif_family_table.csv")

# Nested thresholds: p<=1e-4 is a superset of p<=1e-5 is a superset of p<=1e-6,
# so each tier can be obtained by further filtering the same loaded data
# (FIMO itself was run at output threshold 1e-4, i.e. 1e-4 is everything FIMO kept).
p_thresholds <- c(1e-4, 1e-5, 1e-6)
baseline_threshold <- 1e-4   # "current" analysis threshold, used for rank-stability comparisons

# Cardiomyocyte expression gene sets already exported by
# Objective1_Expressed_genes_GTEx_TSS.Rmd (GTEx v9 snRNA-seq, Myocyte cells,
# Heart tissue). Re-used as-is, not re-derived.
expr_nonzero_csv <- file.path(gtex_dir, "GTEx_non_zeros_ensembl.csv")
expr_top20_csv   <- file.path(gtex_dir, "GTEx_top20%_ensembl.csv")

promoter_windows_nonzero_bed <- file.path(out_data, "promoter_windows_nonzero_expressed.bed")
promoter_windows_top20_bed   <- file.path(out_data, "promoter_windows_top20_expressed.bed")

message("Config loaded. Output directory: ", out_dir)
