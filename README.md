# PPRE cistrome atlas

Code and results for the manuscript *A comprehensive atlas of the predicted human PPRE cistrome decodes transcription factor target heterogeneity in cardiomyocytes*.

The single-nucleus RNA-seq analysis of the cardiomyocyte subpopulations is in [PPRE_cardiomyocytes_manuscript](https://github.com/elibrouwer/PPRE_cardiomyocytes_manuscript).

## Approach

1. 30 PPAR/RXR motifs (HOCOMOCO, JASPAR and CIS-BP) are aligned and merged into consensus logos.
2. The motifs are scanned genome-wide with FIMO (MEME Suite 5.5.7, GRCh38, p <= 1e-4, order-0 genome background).
3. The hits are annotated with ChIPseeker (3 kb upstream to 1 kb downstream of the TSS) and overlapped with genes expressed in GTEx cardiomyocytes.
4. Every gene gets a PPRE score, which is the FIMO score weighted by an exponential decay of the distance to the TSS. The sum, mean, median and maximum per gene are calculated.

## Update: 42 human-source motifs (branch `human42-motif-update`)

The 30-motif set described above was replaced by a systematically assembled set. All PPAR, RXR and PPAR:RXR position weight matrices were retrieved from JASPAR 2026 (CORE, all versions), HOCOMOCO (CORE, v13, identical to v14 for these TFs) and CIS-BP 3.10 (direct evidence, human TFs; three Transfac matrices are not published). This gave 71 motifs (`motifs/PPAR_RXR_all_motifs_combined.meme`, inventory in `results/motif_selection/motif_inventory.csv`).

Motifs whose data come from mouse, or whose species is unclear, were excluded (`results/motif_selection/motif_inventory_excluded_mouse_unclear.csv`). The species of each motif comes from JASPAR and HOCOMOCO annotation; for CIS-BP it is inferred from the source study and flagged as such in `species_basis`. **The final set is 42 human-source motifs** (11 JASPAR, 6 HOCOMOCO, 25 CIS-BP; `motifs/human42/`, `results/motif_selection/motif_inventory_human_only.csv`).

Quality and redundancy of the set were characterised before scanning (`scripts/00_motif_selection/`): fraction of exact zeros (sparse matrices are flagged), half-site architecture (DR1, IR0, IR3, single half-site), and pairwise similarity with clustering (`motif_inventory_gates.csv`, `motif_pairwise_similarity.csv`, `gate2_redundancy_heatmap_human42.png`). All 42 motifs are scanned and every hit is kept; cluster assignments are provided so that redundancy can be handled when scoring (16 clusters at r >= 0.98, `scoring_set_r098_representatives.csv`).

| Step | File |
|---|---|
| Motif retrieval and selection | `scripts/00_motif_selection/` (`jaspar_query.py`, `build_motif_set.py`, `gates.py`, `plot_gate2.py`) |
| Alignment and consensus. The rule-based ("manual") placement reproduces the hand-curated start columns and orientations for the 15 motifs shared with the 30-motif set (15/15). The unsupervised `DNAmotifAlignment()` result depended on the order in which motifs were passed in when it was allowed to reverse-complement (1 to 10 flipped motifs, consensus 3.1 to 6.3 bits), so motifs are oriented by the half-site rule first and the aligner runs with `revcomp = FALSE`; thresholds 0.2 to 0.6 then give identical results and 0.5 is kept (`alignment_threshold_and_order_test.txt`) | `scripts/01b_motif_alignment_human42.Rmd`, `results/motif_alignment_human42/` (`manual_alignment_human42_stack.pdf`, `aligned_motifs_human42_threshold0.5_preoriented.pdf`, `PPRE_architecture_validation_human42_preoriented.pdf`, `manual_alignment_human42.csv`) |
| FIMO | `motifs/fimo_batches/run_fimo_human42_unmasked.sh`: unmasked GRCh38 primary assembly, `--thresh 1e-4`, `--max-stored-scores 20000000`, order-0 genome background, 11 batches of 4 motifs (MEME Suite 5.5.7). 45,704,060 hits in total |
| FIMO summary per motif | `scripts/06_fimo_human42_summary.py`, `results/fimo_human42/fimo_hits_per_motif.csv` (the raw `fimo.tsv` files are about 3.4 GB and are not stored here) |
| ChIPseeker annotation and PPRE score | `scripts/03b_fimo_chipseeker_PPRE_score_human42.Rmd` (generated from script 03 by `scripts/00_motif_selection/make_03b_human42.py`; the only changes are the motif set and `sameStrand = FALSE`, so hits on either strand are assigned to the nearest TSS) |
| Final motif table | `results/motif_selection/Table_S1_final_42_human_source_motifs.csv` |

Output of `scripts/03b_...` is in `results/human42_pipeline/fimo_analysis/` (per-motif table `motif_results_table.xlsx`, motif overlap heatmaps, statistics plots, GTEx overlap `Motifs_GTEX_cm_expressed.csv`: 15,713 cardiomyocyte-expressed genes with at least one promoter hit, PPRE score histograms, ORA). The per-motif ChIPseeker tables (1.5 GB) and the bed file with all 45.7 million hits (2.3 GB) are not in the repository. The per-motif hit and gene counts are not comparable with the earlier 30-motif tables: this run uses the unmasked genome, 42 motifs and `sameStrand = FALSE`.

Note that very short motifs cannot reach small p-values: for example the 8 bp motif `RXRA_M08962_3.10` has no hit at p <= 1e-5, so per-motif p-value tiers are not comparable across widths.

The sections below describe the earlier 30-motif run; scripts 02 to 05 have not been rerun for the new motif set unless stated in the results folders.

## Scripts

Run the scripts in this order, from the root of the repository.

| Script | Content |
|---|---|
| `scripts/01_motif_alignment.Rmd` | motif import, unsupervised and manual alignment, consensus logos, MEME export for FIMO |
| `scripts/02_expressed_genes_GTEx_cardiomyocytes.Rmd` | genes expressed in GTEx heart cardiomyocytes (non-zero and top 20%) |
| `scripts/03_fimo_chipseeker_PPRE_score.Rmd` | ChIPseeker annotation of the FIMO hits, motif overlap (Jaccard), motif statistics, overlap with GTEx and the HCM/PLN H3K27ac regions, GO over-representation, PPRE score |
| `scripts/04_UCSC_tracks_and_TSS_metaplot.Rmd` | UCSC custom tracks, FAM53A multi-track figure, TSS density with a shuffled-motif control |
| `scripts/05_PPARA_chipseq_overlap/` | overlap of the promoter hits of CPT1A, CPT1B, HADHA, HADHB and ACADVL with PPARA ChIP-seq peaks |

## Results

Results are in `results/`: `pan_ppre/` (consensus motif and motif similarity), `tables/`, `figures/`, `beds/` and `PPARA_chipseq_overlap/`.

## Input data

The scripts read the input data from a `data/` folder that is not in this repository:

- `data/PPAR_RXR_motifs_Hocomoco_Jaspar/`: motif files and the FIMO output (`FIMO_results_no_max/<batch>/fimo.tsv`)
- `data/GTEX_single_nuclear_data/`: GTEx snRNA-seq atlas (`GTEx_8_tissues_snRNAseq_atlas_071421.public_obs.h5ad`)
- `data/downloads/`: ChIP-seq peaks (GSM4748812 and ChIP-Atlas `Oth.ALL.05.PPARA.AllCell.bed`)

## Software

R with ChIPseeker, GenomicRanges, rtracklayer, TxDb.Hsapiens.UCSC.hg38.knownGene, org.Hs.eg.db, clusterProfiler, motifStack, universalmotif, DiffLogo, pheatmap, ggvenn, ggplot2, Gviz and ComplexHeatmap. FIMO comes from the MEME Suite.

## Cluster propensity (exploratory, after Wang, Wang & Zang, NAR 2025)

`scripts/11a*.py` (gaps, GC windows), `11b_motif_cluster_propensity.py` (signed KS statistic of nearest-downstream-site distances vs. controls) and `11c_plot_cluster_propensity.py` (figure), run on the 42-motif FIMO hits at p ≤ 1e-6 and p ≤ 1e-4. Controls: uniform outside N-gaps and the ENCODE hg38 blacklist v2, repeat-matched, and GC-matched (1 kb windows). Results in `results/human42_pipeline/clustering/`. After GC matching, most of the apparent clustering disappears (all |CP| ≤ 0.03 at p ≤ 1e-4 except the IR0/palindromic motifs, whose negative CP is an artefact of self-overlapping hits). 10 iterations at p ≤ 1e-6, 5 at p ≤ 1e-4 (paper: 100). `gc_bins_1kb.tsv` (69 MB) is regenerated by `11a2_gc_bins.py` and not tracked.

## Update: score validation, final score and rerun of the pipeline for the 42 motifs

Scripts in `scripts/08_ppre_score_v4/` (see the README there), design records and tables in `results/ppre_score_v4/`, pipeline outputs in `results/human42_pipeline_rerun/`.

* **The original score** (FIMO score weighted by exp(-distance/3000), summed per gene) mostly measures promoter size: its rank correlation with promoter length is 0.64 and with the number of TSSs 0.72.
* **Final score:** the number of PPRE sites in the promoter (hits p <= 1e-5 within -3 kb to +1 kb of a TSS, hits within 10 bp merged), as a percentile among genes of similar promoter length (deciles) and GC (tertiles). It is fixed-weight, needs no ChIP-seq or tissue data and removes the size artefact (correlation 0.02 with length, 0.09 with TSS count, 0.14 with GC). All variants tried (wider windows, chromatin, H3K27ac, conservation, dimer support, motif subsets) were tested against rules written before the run; the results, including negative ones, are in `results/ppre_score_v4/DESIGN_V4.md`.
* **Reading:** gene-level prediction of ligand response is weak (AUROC about 0.53 to 0.6 in expressed hepatocyte genes) and single-gene ranks are uncertain (median 95% interval about 72 percentile points), so the score is meant for groups of genes; the clearest support is at hit level for the retained motifs.
* **Rerun for 42 motifs:** 45,555,607 hits (p <= 1e-4), merged into 9,801,562 regions covering 6.10% of the genome; 3,540,037 promoter hits (ChIPseeker, -3 kb to +1 kb, same strand); 19,418 protein-coding genes with a promoter PPRE, 15,734 of them expressed in cardiomyocytes. Boxplots, Q-Q plots and Shapiro-Wilk tests (the data are not normal, so rank-based tests are used) and hits in repeats per motif (54% of all hits) are in `results/human42_pipeline_rerun/`.
* The per-motif ChIPseeker tables (about 800 MB) and the hit ranges are not stored in the repository.
