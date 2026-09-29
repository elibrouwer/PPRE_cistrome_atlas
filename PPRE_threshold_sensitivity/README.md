# PPRE threshold-sensitivity analysis

Reproducible threshold-sensitivity analysis for the genome-wide PPAR::RXR
response element (PPRE) FIMO scan, covering full heterodimer, PPAR half-site,
and RXR half-site motifs, and their overlap with cardiomyocyte-expressed
gene promoters.

## Scope actually covered, and why

**Threshold tiers analyzed: p<=1e-4, 1e-5, 1e-6 only.** The originally
requested 5e-4 and 1e-3 tiers are **not computable from any FIMO output that
currently exists in this project**. Every complete genome-wide FIMO run
(`FIMO_results_no_max`, `FIMO_results_genome_masked`,
`FIMO_results_promoter_ownbg`) was itself run with FIMO's own `--thresh 1e-4`
— i.e. FIMO discarded any candidate weaker than p=1e-4 *before* writing
output, so those hits were never scored or stored and cannot be recovered
after the fact from existing files. This was confirmed against the project
owner before proceeding (2026-09-25): rather than launch a new genome-wide
FIMO run at p<=1e-3 now (estimated ~10x the hits of the current run, i.e.
tens of millions of additional rows per motif shard and a long runtime), the
decision was to report this as an explicit limitation and hand over the
exact commands to extend the analysis later (below), rather than spend that
compute during this pass.

### To extend to p<=1e-3 later

Run (via the same environment/host that produced `FIMO_results_no_max` —
confirmed as MEME Suite 5.5.7, reachable via `bash` in a separate session)
one command per motif shard, identical to the existing `FIMO_results_no_max`
commands except for `--thresh` and the output directory:

```bash
cd "data/PPAR_RXR_motifs_Hocomoco_Jaspar"
for shard in 1_4 5_8 9_12 13_16 17_20 21_22 23_26 27_30; do
  fimo --thresh 1e-3 --max-stored-scores 20000000 \
       --bgfile genome_bg_order0.txt \
       --oc "FIMO_results_thresh1e-3/${shard}" \
       "RXR_PPAR_combined_files_${shard}.meme" \
       Homo_sapiens.GRCh38.dna.primary_assembly.fa
done
```

Keep `--max-stored-scores 20000000` (or raise it further and check the FIMO
log for a "scores discarded" warning) — a looser `--thresh` produces
substantially more hits per motif, so silently truncating at the default
100,000-score buffer would reproduce exactly the problem this analysis
exists to avoid (see "Why the original FIMO/ run is excluded" below). This
will need noticeably more disk (`FIMO_results_no_max`, at p<=1e-4, is
already ~2.5GB across 8 shards) and runtime than the existing run. Once
done, point `fimo_dir` in `00_config.R` at the new directory, add `1e-3` and
`5e-4` (5e-4 requires its own separate FIMO run — a `--thresh` looser than a
tier does not retroactively produce that tier without rerunning; if both are
wanted, run at `--thresh 1e-3` once and both 1e-3 and 5e-4 fall out of that
same output, exactly like 1e-4/1e-5/1e-6 fall out of the current run) to
`p_thresholds`, and re-run `03_load_fimo_hits.R` onward (delete the cached
`data/all_hits_p1e-4.rds` first, and rename appropriately).

## Why `FIMO_results_no_max` (not the other 3 runs)

Four genome-wide/region FIMO runs exist in this project; only one is used
here:

| Run | Used here? | Why / why not |
|---|---|---|
| `FIMO/` (original, 2025-03) | No | **Truncated.** ~18x fewer hits than the untruncated rerun for the same motif shard (confirmed by direct `wc -l` comparison), almost certainly from FIMO's default `--max-stored-scores 100000` never having been overridden — no log/`fimo.xml` survives to confirm the exact command. Also only covers 22 of the current 30 motifs. |
| `FIMO_results_no_max/` (2026-09-16) | **Yes** | Unmasked GRCh38 primary assembly, genome-wide order-0 background, `--max-stored-scores 20000000` (no truncation), all 30 motifs. Same run `PPARD_threshold_validation/` independently identified as authoritative after comparing all four. |
| `FIMO_results_genome_masked/` (2026-09-16) | No | RepeatMasked genome — different sequence composition (repeat regions removed) makes hit counts and p-values not directly comparable to the unmasked run without re-deriving everything on the masked genome too. Kept as a separate, uncombined dataset per the project's own convention. |
| `FIMO_results_promoter_ownbg/` (2026-09-23) | No | Scanned only a promoter FASTA directly, using a **promoter-composition background** (`promoter_bg_order0.txt`) instead of the genome-wide background the other runs use. p-values from a different background model are not comparable to genome-wide p-values at the same nominal cutoff, even though the nominal threshold (1e-4) is the same number. |

**Two existing validation scripts still point at the truncated original
`FIMO/` run** — `Objective1_fimo_chip_validation.R` and
`Objective1_fimo_auc_validation.R` (both dated 2026-09-22/23) hard-code
`fimo_dir <- file.path(base_dir, "FIMO")`. Their PPARD.H13CORE.0.PSM.A
AUC/enrichment numbers should be treated as unreliable until pointed at
`FIMO_results_no_max` instead — **not changed by this analysis**, per the
instruction to preserve existing outputs; flagged here for a future fix.

## Pipeline

Run in order (`Rscript.exe` used directly; R 4.5.2, packages confirmed
installed: data.table, dplyr, readr, tidyr, tibble, GenomicRanges,
GenomeInfoDb, IRanges, S4Vectors, TxDb.Hsapiens.UCSC.hg38.knownGene,
org.Hs.eg.db, rtracklayer, ggplot2, scales):

```bash
Rscript 01_motif_family_table.R
Rscript 02_promoter_windows_expressed_genes.R
Rscript 03_load_fimo_hits.R
Rscript 04_dedup_intervals_by_threshold.R
Rscript 05_promoter_overlap_gene_ranking.R
Rscript 06_family_contribution.R
Rscript 07_chip_validation.R
Rscript 08_summary_tables_and_plots.R
```

`00_config.R` is sourced by every script (not run directly) and holds all
paths/settings in one place.

1. **`01_motif_family_table.R`** — exports the motif-family classification
   (heterodimer / PPAR half-site / RXR half-site) as a standalone,
   machine-checked CSV. This classification previously existed only inline
   in one Rmd chunk (`Objective1_motif_alignment_CISBP_additions.Rmd`,
   `manual_alignment` data frame) — transcribed here verbatim and
   cross-checked against the actual `MOTIF` lines in the 8 `.meme` shard
   files, so a future motif-set change that isn't reflected here fails
   loudly (`stop()`) instead of silently going stale.
   Counts: **3 heterodimer, 15 PPAR half-site, 12 RXR half-site** (30 total).
2. **`02_promoter_windows_expressed_genes.R`** — builds -3000/+1000-of-TSS
   promoter windows (union across UCSC knownGene transcript isoforms per
   gene; same rule as `PPARA_chipseq_overlap/01_get_promoter_windows.R` and
   the ChIPseeker `tssRegion=c(-3000,1000)` calls elsewhere in this
   project), for the two cardiomyocyte-expressed gene sets already exported
   by `Objective1_Expressed_genes_GTEx_TSS.Rmd` (GTEx v9 snRNA-seq, Myocyte
   cells, Heart tissue; re-used as-is, not re-derived): **nonzero
   expression** (any detectable expression) and **top20** (top 20% by mean
   expression).
3. **`03_load_fimo_hits.R`** — loads all `FIMO_results_no_max` hits
   (32,295,620 rows across 30 motifs and 8 shards, everything FIMO kept at
   its own `--thresh 1e-4`) into one compact table. Because the 1e-5 and
   1e-6 tiers are strict subsets of what FIMO already computed at 1e-4, a
   single load supports all three tiers — no re-reading of the raw TSVs per
   threshold.
4. **`04_dedup_intervals_by_threshold.R`** — for each threshold: (a) raw
   motif-call counts per motif/family (JASPAR/HOCOMOCO/CIS-BP overlapping
   calls at the same locus **not yet** collapsed), and (b) deduplicated
   genomic intervals (`GenomicRanges::reduce`, strand-ignored — a PPRE-type
   element is the same site whichever strand FIMO's best-scoring
   orientation landed on), retaining which motif(s) and family(ies) support
   each merged interval. Heterodimer / PPAR-half / RXR-half hits are kept
   distinguishable throughout, never silently pooled.
5. **`05_promoter_overlap_gene_ranking.R`** — overlaps deduplicated
   intervals with both expressed-gene promoter sets, separately for 4
   motif "views" (heterodimer-only, PPAR-half-only, RXR-half-only,
   combined/any). Computes a per-gene score (`sum(-log10(best_pvalue))`
   across its promoter-overlapping intervals) and rank at each threshold,
   then Spearman rank correlation and top-100/top-500 gene-list overlap
   against the **p<=1e-4 baseline** ("the current analysis" in the original
   brief) — per view and per gene universe, kept separate rather than
   pooled.
6. **`06_family_contribution.R`** — gene-level heterodimer-vs-half-site
   breakdown: of the genes assigned a candidate PPRE (combined view), how
   many have that support from a full heterodimer motif at all, vs. only
   from a PPAR or RXR half-site with no heterodimer hit anywhere in the same
   promoter.
7. **`07_chip_validation.R`** — genome-wide fold-enrichment validation
   against real ChIP-seq (`PPARA_chipseq_overlap/chipatlas_liver_peaks_raw.bed`,
   ChIP-Atlas PPARA peaks aggregated across 7 liver-context experiments,
   reduced to non-redundant loci — same step
   `06_make_liver_centrimo_fasta.R` already used). Uses the same
   fold-enrichment-vs-genome-background + one-sided binomial test method
   already established in this project's
   `Objective1_fimo_threshold_enrichment.R` — **deliberately not** the more
   elaborate permutation/blacklist-masking framework built for
   `PPARD_threshold_validation/`, per an explicit decision not to extend
   that pipeline in this pass (see "Related prior validation work" below).
   **This is INDIRECT validation**: PPARA antigen, liver tissue — not
   cardiomyocyte, and not PPARD/PPARG/RXR ChIP. No RXR ChIP-seq and no
   cardiac-tissue PPAR ChIP-seq exists anywhere in this project (confirmed
   during the audit); RXR half-site hits are scored against the same PPARA
   peaks only on the rationale that they can sit in the same DR1 element a
   PPARA:RXR heterodimer ChIP would capture, not because RXR was the ChIP
   antigen. Absence from a ChIP peak is not evidence of no binding — different
   tissue, different condition, and ChIP peak-calling itself has a
   sensitivity floor.
8. **`08_summary_tables_and_plots.R`** — assembles
   `results/MASTER_results_by_cutoff_and_family.csv` and 4 summary figures.

## Related prior validation work (cited, not extended here)

`PPARD_threshold_validation/` (in the data directory, not this repo) already
ran a rigorous, pre-registered 11-point threshold sweep for
`PPARD.H13CORE.0.PSM.A` specifically, against 8 real HUVEC PPARD ChIP-seq
experiments (GEO SRA098847), with permutation nulls and ENCODE
blacklist/assembly-gap masking. **Its headline finding is directly relevant
here and should inform interpretation of this analysis's threshold
recommendation**: at the project's standard p<=1e-4, FIMO hits were
*significantly depleted* (not enriched) in real ChIP peaks pooled across
conditions (fold enrichment 0.25-0.29x, empirical P=0.001), and **no single
genome-wide threshold met its pre-registered "good enough" bar** — though 2
of 8 individual experiments showed strong, real, monotonic enrichment (up to
25-35x). Per an explicit decision made with the project owner, this
project's specific pipeline/methodology was **cited, not re-run or extended**
here; see `PPARD_threshold_validation/FINAL_INTERPRETATION.md` for the full
result.

## Known gaps / traceability notes (not fixed here, flagged for awareness)

- `promoter_regions_tss3000_1000.fa` / `_repeatmasked.fa` (used by
  `FIMO_results_promoter_ownbg`, not by this analysis) have no generating
  script found anywhere in either directory tree — naming is consistent with
  the `-3000/+1000` rule but this was not confirmed by inspecting a script.
- No genome-provenance record (exact Ensembl GRCh38 release/date) exists
  anywhere in the project for the `Homo_sapiens.GRCh38.dna.primary_assembly.fa`
  file FIMO was run against.
- `R_code/README.md` says "21 motifs"; the actual counts are 22 (original
  truncated run) and 30 (current, used here) — never reconciled.

## Part 2: gene-scoring / enriched-list methods (scripts 09-14)

Added after the threshold-sensitivity analysis (above) to answer a follow-up
question: given that promoter overlap is nearly saturating and p-value
threshold doesn't track real occupancy (both established in Part 1), how
should candidate genes actually be scored and ranked? Five established
methods were implemented and validated -- not invented ad hoc -- against
independent, community-curated PPAR gene sets from MSigDB (`msigdbr`) via
`fgsea`, the standard GSEA implementation (same framework this project
already uses, see `Objective1_results/GSEA_*_GSEA_results.csv`).

- **`09_dr1_hexamer_geometry.R`**: PWM-independent scan of promoter DNA for
  the classical AGGTCA-like half-site hexamer (Umesono et al. 1991 direct-
  repeat rule), looking for two hexamers on the same strand with the
  canonical 1bp DR1 spacer. Honest negative-ish result: the gap-distance
  distribution shows no clean peak at 1bp (a 6bp pattern with 1 mismatch
  allowed recurs too often by chance), and heterodimer-PWM-supported genes at
  p<=1e-4 show a DR1-hexamer rate (49.7%) statistically indistinguishable
  from background (48.9%). Only the strictest heterodimer tier (p<=1e-6, the
  rare/strong hits) shows real enrichment (79.8%) -- a useful cross-check,
  not a strong standalone score.
- **`10_reference_gene_sets.R`**: caches 4 external MSigDB gene sets for
  validation -- `PPAR_DR1_Q2` (TRANSFAC, genes with an independently
  predicted PPAR DR1 element -- the closest methodological analogue to this
  project's own approach), `KEGG_PPAR_SIGNALING_PATHWAY`,
  `REACTOME_REGULATION_OF_LIPID_METABOLISM_BY_PPARALPHA`,
  `SANDERSON_PPARA_TARGETS` -- plus this project's own 5 ChIP-checked PPARA
  fatty-acid-oxidation genes (CPT1A, CPT1B, HADHA, HADHB, ACADVL).
- **`11_gene_scoring_methods.R`**: Fisher's combined probability test
  (Fisher 1925) and Stouffer's weighted Z-method (Stouffer et al. 1949,
  heterodimer loci weighted 3x half-site loci), per gene per motif family.
  **Fixed a real numerical bug**: with a median of ~47 supporting loci per
  gene, linear-scale combined p-values underflow to exact 0 for most genes,
  destroying rank information -- recomputed on the log scale throughout.
  Also fixed `qnorm(1-p)` in the Stouffer step, which silently returns `Inf`
  for any p below ~1e-16 (`1 - p` rounds to exactly `1.0` in double
  precision) -- replaced with `qnorm(p, lower.tail=FALSE)`.
- **`13_size_normalized_scoring.R`**: **the important fix.** The naive scores
  from 11 turned out to be dominated by promoter-window size (Spearman
  rho=0.862 between promoter_bp and raw loci count) -- genes with many
  alternative TSSs (e.g. `MACF1`: 21 windows, 109,608bp unioned promoter)
  accumulate more hits by pure chance, nothing to do with real biology. Fixed
  with a size-normalized Poisson enrichment test (observed vs. expected loci
  given that gene's own promoter size and the genome-wide per-bp rate for
  that motif family) -- the same logic this project already uses for its
  ChIP fold-enrichment tests, applied per-gene instead of per-threshold.
  Post-fix correlation with promoter size: 0.037.
- **`12_method_validation_and_aggregation.R`** / **`14_final_validation_and_aggregation.R`**:
  validate every method (naive and size-corrected) via `fgsea` against the 4
  external gene sets, then combine methods that clear an explicit bar
  (padj<0.05 for >=1 of the 3 substantial gene sets) via RobustRankAggreg
  (Kolde et al. 2012). **Finding: the RRA consensus performed WORSE than the
  single best method** (NES dropped from 1.35/padj=3.8e-7 to
  non-significant) -- the 4 family-specific Poisson scores are too
  redundant to combine profitably, so the consensus blend was NOT used as
  the final answer; reported as an honest negative result instead.

### Recommended method and final list

**`poisson_neglog10p_heterodimer`** (size-normalized Poisson enrichment of
heterodimer-motif loci) is the single best-validated method: NES=1.35,
padj=3.8e-7 against `PPAR_DR1_Q2` -- the one external gene set that is
methodologically comparable to this project's own approach (an independent,
TRANSFAC-based DR1-motif prediction, not our own FIMO output). It does not
reach significance against `KEGG_PPAR_SIGNALING_PATHWAY` or
`SANDERSON_PPARA_TARGETS` (broader pathway-membership sets, less specific to
promoter motif architecture -- an honest limitation). It also doesn't reach
significance against the REACTOME lipid-metabolism set; only the naive,
size-biased method does there, but REACTOME-set genes have significantly
larger promoters than average (Wilcoxon p=0.02) -- that "win" for the naive
method is itself suspected to be the same size artifact reappearing, not a
genuine advantage, and should not be read as evidence for using the naive
method instead.

The 5 project-confirmed PPARA fatty-acid-oxidation genes land at ranks
1,119 (HADHB), 2,098 (CPT1A), 2,280 (ACADVL), 2,629 (HADHA), and 3,076
(CPT1B) out of 15,703 -- consistently in the top 7-20%, a real and
consistent improvement over chance without claiming these are "the" top
hits (sequence motif density is only one contributor to real PPARA
responsiveness; co-activator recruitment and chromatin state are not
captured here).

Final output: `results/FINAL_recommended_enriched_gene_list.csv` (all
15,703 genes ranked by `poisson_neglog10p_heterodimer`), and the comparison
figure `figures/fig6_final_method_validation_NES.png`.

## Key results (2026-09-25 run)

- **Deduplication collapse**: 32,295,620 raw hits (p<=1e-4) -> 10,067,664 merged intervals (3.20x); 4,011,452 -> 1,873,171 (2.13x) at p<=1e-5; 619,098 -> 393,003 (1.57x) at p<=1e-6.
- **Half-site motifs dominate raw hit counts** at every threshold by roughly an order of magnitude over heterodimer motifs (e.g. p<=1e-4: 15.09M PPAR-half + 14.38M RXR-half vs. 2.83M heterodimer raw hits) — expected, since a heterodimer call requires both half-sites AND the correct spacer simultaneously.
- **Promoter overlap is near-saturating and threshold-insensitive as a yes/no signal**: 15,703/15,703 (100%) of non-zero-expressed genes have >=1 candidate-PPRE-containing promoter interval at p<=1e-4, still 12,659/15,703 (80.6%) at p<=1e-6. A bare presence/absence call over this threshold range carries little discriminating information on its own — the continuous per-gene score and which motif family supports it matter far more than whether a gene clears the cutoff at all.
- **Heterodimer support erodes far faster than half-site support as the threshold tightens**: of genes assigned a candidate PPRE, 98.1% have direct heterodimer motif evidence at p<=1e-4, but only 15.0% still do at p<=1e-6 — i.e. 85.0% of the "high-confidence" p<=1e-6 gene set is actually supported ONLY by an isolated half-site, not a full heterodimer call. Tightening the threshold does not make the heterodimer evidence more selective; it just removes it, leaving half-site-only calls as the overwhelming majority of what remains.
- **Gene-ranking stability** (Spearman rho vs. the p<=1e-4 baseline, nonzero-expressed genes, combined motif view): 0.92 at p<=1e-5, 0.71 at p<=1e-6. Heterodimer-only ranking is consistently the least stable view (0.62, then 0.29); RXR-half is the most stable (0.92, then 0.71).
- **ChIP validation (indirect — liver PPARA, not cardiac)**: modest enrichment at p<=1e-4 (fold 1.07-1.15x across motif views, statistically significant only because of very large N) that does NOT improve, and in fact drops below chance (fold 0.75-0.93x), at both stricter thresholds for every motif family. This is the same qualitative pattern `PPARD_threshold_validation/` found independently for PPARD specifically: a stricter FIMO p-value does not track with better real ChIP-seq recovery in this dataset.
- **No single threshold in {1e-4, 1e-5, 1e-6} can be empirically justified as biologically optimal from this analysis** — consistent with the cited PPARD-specific prior finding. See the recommendation given to the user for the suggested broad-candidate (p<=1e-4) vs. stricter sequence-priority (require heterodimer-motif support specifically, not just a lower p-value) framing.

Full numbers: `results/MASTER_results_by_cutoff_and_family.csv`, `results/rank_stability_vs_baseline.csv`, `results/family_contribution_by_threshold.csv`, `results/chip_validation_PPARA_liver_by_threshold.csv`.

**Disk note**: `results/deduplicated_intervals_p1e-04.csv` is ~1.3GB (10M rows, one per merged interval); the 1e-5/1e-6 versions are 213MB/41MB. Delete these three per-interval CSVs if disk space is tight — everything else in `results/` is small and derived from them plus the cached `data/all_hits_p1e-4.rds` and `data/deduplicated_intervals_by_threshold.rds`, which can regenerate them.

## Outputs

- `data/` — `all_hits_p1e-4.rds` (cached, all FIMO_results_no_max hits),
  `promoter_windows_{nonzero,top20}_expressed.bed`,
  `deduplicated_intervals_by_threshold.rds`.
- `results/` — per-threshold and summary CSVs, including
  **`MASTER_results_by_cutoff_and_family.csv`** (the requested table of
  results by cutoff and motif family) and `motif_family_table.csv`.
- `figures/` — 4 summary plots (raw hit counts, genes assigned, rank
  stability, ChIP fold-enrichment tradeoff).

**Nothing in this folder modifies or overwrites any pre-existing analysis
output elsewhere in the project** — all inputs (`FIMO_results_no_max`,
the GTEx expressed-gene CSVs, the ChIP-Atlas liver peaks BED) are read
read-only; all new output goes under this dedicated
`PPRE_threshold_sensitivity/` folder.
