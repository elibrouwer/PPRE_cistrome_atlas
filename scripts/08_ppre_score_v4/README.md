# Score validation and final score (v4 work)

Scripts for the evaluation of PPRE gene scores on human data, the pre-registered tests of score variants, and the rerun of the motif pipeline for the 42 motifs. Design, pre-set rules and results of every test are in `results/ppre_score_v4/DESIGN_V4.md` (Addenda 1 to 14) and `RESULT_GATE1.md`.

| Script | Content |
|---|---|
| `gate1_baselines.py` | baselines on human hepatocyte agonist responses (GSE17251, GSE53399) |
| `v4_windows.py`, `v4robust.py` | pairs of hits and genes (10 and 50 kb windows), liver ATAC gate, robustness grid, per-gene intervals, leftover confounds |
| `v4_expr_adjust.py`, `v4e_expressed.py`, `v4tf.py` | expression adjustment, expression-aware correction, training-free score |
| `v4cd.py`, `v4k27.py`, `v4k27_posthoc.py`, `v4cardio.py` | conservation, dimer support, H3K27ac layers and the cardiomyocyte test |
| `v4motif.py`, `v4motif2.py`, `v4repeats.py` | motif reliability and subsets (leave-one-dataset-out), repeat-embedded hits |
| `v4simple.py`, `v4easy.py`, `v4easy2.py`, `v4tepic.py` | simpler score definitions and head-to-head comparison with the TEPIC-style sum |
| `export_sitecount_scores.py` | final score for all genes (number of sites, percentile within promoter-length and GC groups) |
| `cardio_labels.R`, `gse262419_labels.R`, `atac_peaks.R` | labels and checks for the cardiomyocyte datasets GSE160987, GSE262419 and GSE178984 |
| `pipeline42_part1.R`, `part1b.R`, `part2.R`, `part3.R`, `part4.R` | per-motif ChIPseeker annotation (one motif at a time), tables, gene lists, heatmaps, genome coverage, over-representation analysis for the 42 motifs |
| `boxplots_42.R`, `boxplots_repeats_42.R`, `boxplots_single_42.R` | boxplots of hits and width per motif, hits in repeats per motif, and every plot and Q-Q plot as a separate file (`Statistics_plots/single_plots/`) |
| `v4_standard.py`, `v4motif3.py`, `v4motif_oi.py`, `v4motif_s0.py`, `v4motif_sizecontrol.py` | motif prioritisation follow-ups: field-standard benchmarking, outcome-independent sets, RXR half-site removal, size-matched control, cardiomyocyte confirmation (Addenda 11 to 15) |
| `v4final.py`, `v4stab.py` | recomputation of robustness, intervals and the ranked list for the site-count percentile score |
| `build_motif_page.py`, `update_pages_addendum15.py` | build the HTML summary pages (not needed to reproduce results) |
| `v3_human.py`, `common.py` | helpers imported by the v4 scripts |

The scripts contain absolute paths of the local analysis folders and expect the cache of hits (`PPRE_CACHE_DIR`) and ChIP-seq peak sets (`HPEAKS`) described in `DESIGN_V4.md`; adjust the paths before running. The per-motif ChIPseeker tables (about 800 MB) and the hit ranges are not stored here.
