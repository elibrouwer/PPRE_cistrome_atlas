# Addendum to PPRE_score_validation_summary.md (written 2026-10-05, session "PPAR/RXR motif files collection")

Purpose: (1) correct numbers in the summary that were affected by a coding bug, (2) add experiments run after that summary was written (idea 2 motif + expression, idea 3 distal calibrated model, idea 6 6-mer sequence model, corrected heart-ATAC and distal-kernel tests). To be merged into the summary once the other sessions have finished.

## 1. The bug and what it affects
In scripts 16b, 22, 24 and 28 the distance weight `exp(-|d|/tau)` was passed as a global-length vector inside a per-gene `max()`, so it was recycled over the whole table and the decay was not applied per hit (the window cut-offs were applied correctly). Fixed by storing the weight as a column (`gagg()` helper). All dependent analyses were rerun (16b, 17, 18, 20, 21, 22, 24, 25, 26, 27b, 28, 29b); CSVs in `results/human42_pipeline/` are the corrected ones (commits d1357ae, ae13764).
Not affected: ChIPseeker hit-level scores (13 to 15), window-level peak/calibration work (19, 23), site-level literature test (25 site level), the first ChIP-calibrated gene score after its own earlier fix (20), the 6-mer model (29).
Other sessions: scripts 28a/28b (p-cutoff and distal calibrated scores) were written in parallel; if they use a global-length weight vector inside a grouped max they carry the same bug. Their output folder `results/human42_pipeline/score_improvement/` was empty when checked.

## 2. Corrections to the summary
| Summary statement | Corrected |
|---|---|
| "Wide windows (+-100 kb, raw hits) saturated, AUROC 0.50" | With the decay applied: PPARgene 0.52 to 0.53, heart 0.56 to 0.57, overall 0.54 to 0.55 (all-motif promoter 0.54, heterodimer promoter 0.55). Not at chance, not better than the heterodimer promoter score |
| "Heart-H3K27ac-linked distal weighting: small, inconsistent gains" | Best overall of the hand-built variants: ac_tau20k_sum 0.58 (PPARgene 0.558, heart 0.602), ac_tau20k_max 0.57; rank average with the promoter score +0.022 on PPARgene (0.006 to 0.039) and +0.043 on heart sets (0.003 to 0.089) vs the all-motif promoter score. Other differences have intervals including 0 |
| "Heterodimer-only: PPARgene 0.57, heart 0.55" | PPARgene 0.551, heart 0.553 (overall 0.552); promoter-peak genes 0.566 |
| Section 5 subset means | PPARgene+heart: all 42 0.537, PPAR-containing 14 0.549, heterodimer 3 0.552, PPAR-only 11 0.548, RXR-only 0.536. Promoter-peak genes: 0.528 / 0.561 / 0.566 / 0.559 / 0.527. Paired (PPARG promoter-peak genes): PPAR14 minus all +0.080 (0.060 to 0.099); heterodimer minus all +0.101 (0.081 to 0.122); PPAR14 minus heterodimer -0.021 (-0.032 to -0.011). PPARgene+heart: all differences within +-0.06, intervals include 0 |
| Section 6 literature gene level (tier 1 / tier 1+2) | all 42 0.563 / 0.568; PPAR14 0.556 / 0.586; heterodimer 0.604 / 0.626; ChIP-calibrated top-3 0.630 / 0.636; expression 0.597 / 0.599. PPARA/PPARD genes (21): all 42 0.476, PPAR14 0.654, heterodimer 0.704, calibrated 0.644, expression 0.566. PPARG genes (19): all 42 0.662, PPAR14 0.529, heterodimer 0.557, expression 0.681. Intervals about +-0.1. Site-level numbers unchanged |
| PPARD myofibroblast, gene level | promoter-peak genes (800): all 42 0.505, PPAR14 0.529, heterodimer 0.523, calibrated 0.521, expression 0.570. Peak <= 10 kb and GW501516-induced (31): 0.529 / 0.529 / 0.520 / 0.536, expression 0.631. Peak-level numbers unchanged |
| ChIP-calibrated model vs hand-built | vs heterodimer promoter: PPARgene +0.023 (-0.003 to +0.049), heart +0.027 (-0.029 to +0.085); vs all-motif promoter: PPARgene +0.048 (0.012 to 0.083) |
| Promoter-peak gold (script 21), heterodimer score | PPARG blood 0.73, others 0.67, liver 0.60, digestive 0.59, adipocyte 0.57; PPARA liver 0.50; RXRA 0.50 to 0.57 |
| Heart ATAC (ENCFF174HPZ) | restricting the heterodimer score to ATAC bins: -0.046 on PPARgene+heart (-0.085 to -0.009), -0.066 on PPARG promoter peaks; all-motif score improves +0.026 / +0.036 on PPARG promoter peaks; rank average with heterodimer ties |

## 3. New experiments
Gold registry: `results/human42_pipeline/score_exploration/gold_registry.rds` (script 27a): families A literature tier 1/1+2, B PPARgene + mouse heart, C promoter-peak genes (PPARG groups, PPARA/RXRA liver, PPARD myofibroblast), D distal-peak-only genes (3 to 10 kb, 10 to 50 kb), E PPARD peak + GW501516-induced. Matched-background AUROC, family means.

**Idea 2, motif + expression (27b)**: rank average of motif score and cardiomyocyte expression. Overall mean: expression 0.569; PPAR14 + expression 0.570; heterodimer + expression 0.576; ChIP-model + expression 0.582; motif alone 0.545 to 0.564. Gain over expression alone: promoter-peak genes +0.030 (0.022 to 0.037); PPARgene+heart +0.029 for ChIP-model + expression (0.006 to 0.051), PPAR14 +0.009 (n.s.); literature +0.012 to +0.052 (intervals include 0); distal-peak -0.003 to -0.009; PPARD induced -0.04 (n.s.).

**Idea 3, distal calibrated model (28)**: adipocyte-trained elastic net applied to 500 bp bins genome-wide, max over bins within +-100 kb with decay. Overall: wide, no GC term, tau 50 kb 0.580 (literature 0.608, PPARgene+heart 0.574, promoter peak 0.560, distal peak 0.509, PPARD induced 0.649); tau 20 kb 0.577; top-3 mean tau 20 kb 0.573; promoter-only calibrated 0.557; heterodimer promoter 0.554; all-motif promoter 0.532. Distal-peak-only genes stay at 0.49 to 0.51 for every score. No paired intervals computed for these differences.

**Idea 6, 6-mer sequence model (29, 29b)**: canonical 6-mer counts of 500 bp windows, elastic net, PPARG peaks vs GC-matched controls. Peak level (AUROC): adipocyte held-out chromosomes 0.85 (motif 0.76, motif + 6-mer 0.89); liver 0.78; digestive 0.77; blood 0.84; HUVEC 0.86 (motif 0.48); PPARD myofibroblast 0.67. Strong positive 6-mers include CpG-rich words and bZIP-like TTGCAA besides AGGTCA/GGGTCA. The high HUVEC transfer shows the model mostly learns generic open-chromatin/regulatory sequence. Gene level (promoter tiles): overall 0.539 vs motif 0.545 to 0.553 and expression 0.570.

## 4. Where this leaves the claim
Unchanged in substance, numbers slightly different: PPAR-containing (or heterodimer) promoter motif evidence, AUROC about 0.55 to 0.65 depending on gold standard, with distal calibrated and expression-combined variants adding about 0.01 to 0.03. No tested idea raised any gold standard above about 0.65. Cardiac validation remains absent. Evaluation caveats as in section 9 of the summary (many variants on the same gold sets, expression bias, study bias).

## 5. To do when the other sessions have concluded
1. Merge section 2 into the summary tables (sections 3, 5, 6) and add section 3 of this addendum.
2. Check 28a/28b and the ceiling-test scripts for the same grouped-weight bug; rerun if affected.
3. Decide the primary score (PPAR-containing 14 vs heterodimer 3 vs ChIP-calibrated) and update the figure (script 12 still shows the median).
