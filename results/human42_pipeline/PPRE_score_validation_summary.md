# PPRE score: what we tested and what it shows (42 human-source motifs)

Status: all numbers below are from the scripts in `scripts/` (11 to 26) and tables in `results/human42_pipeline/`. AUROC = area under the ROC curve; 0.5 = chance.
"Matched" = against genes matched on promoter length and GC (benchmark evaluator), unless stated.

## 1. Where we are

| Part | Status |
|---|---|
| Motif set | 42 human-source motifs (11 JASPAR, 6 HOCOMOCO, 25 CIS-BP): 28 RXR, 11 PPAR, 3 PPAR:RXR heterodimer (MA1148.1, MA1148.2, MA0065.1). Mouse and unclear-species motifs excluded; 16 clusters at r >= 0.98 |
| FIMO | unmasked GRCh38, p <= 1e-4, order-0 background, 45.7 M hits (45.56 M on standard chromosomes), all released |
| Atlas outputs | footprint tables, hits-per-gene histograms (fig7-fig12), alignments, ChIPseeker/GTEx pipeline, all on GitHub main |
| Cluster propensity (NAR 2025 method) | no real clustering once GC is matched; comparable to the published distribution |
| Gene-level PPRE score | weak signal only (see section 3). Manuscript score (median) is at chance |
| Cardiac validation | none available. Only cardiac-context ChIP found (iPSC-cardiac PPARA 317 peaks; HUVEC PPARG) does not work |

## 2. Cluster propensity (Wang, Wang & Zang, NAR 2025 approach)
- Paper (571 JASPAR motifs, Table S1): mean CP 0.015, median 0.023, 68% positive, range -0.237 to +0.669.
- Ours, paper-style control (p <= 1e-4): 74% positive, CP -0.11 to +0.12; GC-matched control: 62% positive but all |CP| <= 0.03 except the IR0/palindromic motifs (artefactual negative CP).
- Four motifs shared with the paper: our CP is close to theirs (correlation 0.99 for GC-matched).
- Conclusion: the PPAR/RXR motifs show no meaningful clustering beyond base composition; not a finding.

## 3. Gene-level score tests (chronological)

| Test | Result |
|---|---|
| GO fatty-acid gene sets, original definitions | median score about 0.5; max-based 0.59 (beta-oxidation); decay 3000 best; raw FIMO score >= -log10 p |
| Confounder control (GC, length, expression) | max score keeps AUROC about 0.58 (OR 1.30/SD, p = 0.03 for beta-oxidation); median has none |
| 16 alternative scores (relative, percentile, cluster-sum, noisy-OR, TRAP-like, top-3 loci, kernels, GC-adjusted, ensemble) | none robustly better; honest best-of selection gain +0.001 |
| Earlier 30-motif benchmark (`benchmark/`, PPARgene 329 targets + mouse-heart PPARA-response sets) | score ceiling about 0.55-0.60; hit count not worse than score |
| Wide windows (+-100 kb, raw hits) | saturated, AUROC 0.50 |
| Heart-H3K27ac-linked distal weighting | small, inconsistent gains |
| Plain hit counts per gene (fig7-fig12 definitions) | 0.51-0.52, worse than the heterodimer FIMO score (-0.051, interval -0.076 to -0.025) |
| Heterodimer-only (3 motifs) near-TSS score | best hand-built: PPARgene 0.57, heart sets 0.55 |
| ChIP-calibrated model (adipocyte PPARG), gene score on promoter tiles | 0.576-0.581 overall = tie with heterodimer-only; +0.062 over all-motif score on PPARgene |
| ENCODE heart ATAC (ENCFF174HPZ; 9,717 merged peaks) weighting | no help (heterodimer score -0.056) |

## 4. ChIP-seq based tests

**Where PPAR peaks are**: 25-38% of PPARG/RXRA/PPARA peaks lie within 3 kb of a TSS (random 13-14%); median 6-14 kb; most binding is distal.

**Peak-level (do hits separate peaks from GC-matched windows?)**
| Peak set | strong heterodimer hit in peaks | controls | AUROC heterodimer motifs |
|---|---|---|---|
| PPARG adipocyte / liver / digestive / blood / others | 5-12% | 0.5% | 0.65-0.75 |
| RXRA liver / others | 2.6% / 5.1% | 0.5% | 0.61 / 0.65 |
| PPARA liver | 0.4% | 0.5% | 0.49 (no motif enrichment) |
| PPARD myofibroblast (E-MTAB-371, 4,531 lifted peaks) | 5.4% | 0.8% | 0.56 (PPAR-containing 14: 0.567) |
| PPARG HUVEC; PPARA iPSC-cardiac | about 0 | 0.45% | no signal; peaks not near TSS either |
- Elastic net on 42 motifs (adipocyte, held-out chromosomes): 0.76; heterodimer max 0.72; all-motif max 0.59; GC alone 0.50. Transfers to liver (0.78), others (0.75), digestive tract (0.73), blood (0.68); not to HUVEC (0.48).

**Gene-level vs genes with a promoter-proximal peak (AUROC)**
- PPARG groups: heterodimer 0.55-0.70, all-motif about 0.51-0.53.
- PPARA liver 0.52, RXRA about 0.50-0.55: at chance. PPARD promoter peak genes: 0.52-0.53 (PPAR-containing), 0.49 (all 42).

## 5. Which motifs should the score use?
- All 42 behaves like the RXR-only subset (identical AUROC), because 28 of 42 motifs are RXR. RXR half-sites are not PPAR-specific (RXR partners with many receptors) and RXR-only motifs predict nothing, even for RXRA peaks (0.51).
- PPAR-containing 14 (11 PPAR + 3 heterodimer): PPARG promoter-peak genes +0.075 over all 42 (interval 0.058 to 0.094); PPARgene+heart sets +0.014 (not significant).
- Heterodimer 3 vs PPAR-containing 14: +0.025 on PPARG peaks only; selection after seeing results, so treat as sensitivity.
- Recommended primary definition: PPAR-containing 14 (chosen on biology). Heterodimer-only as sensitivity result.

## 6. Literature-curated human PPREs (`literature_PPRE_human_prioritised.csv`)
- Gene level (tier 1: 25 genes; tier 1+2: 40 genes): heterodimer 0.62/0.61, PPAR-containing 0.59/0.59, all 42 0.61/0.59, ChIP-calibrated model 0.63/0.64, cardiomyocyte expression alone 0.60/0.60. No paired difference separable from zero (interval about +-0.1). PPARA/PPARD genes (21): 0.63-0.64 vs all-42 0.54; PPARG genes (19): all-42 0.63, het 0.55.
- Site level (motif evidence at the exact validated position vs random genes, same TSS-relative coordinate): tier 1: mean percentile 0.63, 50% in top quarter (expected 0.50 and 25%); tier 1+2: 0.60 and 40%. PPARA/PPARD sites better than PPARG sites.
- Caveats: study bias (well-studied metabolic genes), partly circular (sites found by prediction), expression not matched.

## 7. Why PPARG enriches and PPARA/RXRA do not (interpretation, partly tested)
- Tested: PPARG peaks are strongly motif-driven (10-26x heterodimer hit enrichment); PPARA liver peaks are not; RXRA in between.
- Hypotheses (not tested): PPARA peaks reflect tethering or co-binding (e.g. HNF4A); PPARG data come from strongly ligand-activated, high-expression settings; the PPARG::RXRA matrix may be the best match to its direct DR1 binding.
- Note: cardiomyocytes mainly express PPARA and PPARD, so the PPARG-validated signal may be least informative for heart.


## 7b. Ceiling test (script 27): how well does real ChIP-seq binding predict the same gold sets?
Gene score = strongest ChIP peak in the promoter window or within +-10 kb; same gold sets (4 PPARgene + 4 mouse-heart) and matched-background AUROC.
- Best real-binding scores: RXRA liver +-10 kb 0.635, PPARG adipocyte +-10 kb 0.627, RXRA liver promoter 0.612, PPARG liver +-10 kb 0.607, PPARG adipocyte promoter 0.604.
- PPARA liver ChIP (the isoform of the heart gold sets) only 0.51-0.52; PPARA neural, RXRA digestive/blood, PPARG HUVEC at 0.50.
- Motif scores for comparison: heterodimer 3: 0.562; PPAR-containing 14: 0.552; all 42: 0.540. The heterodimer motif score ranks above 18 of the 28 ChIP-based scores.
- Conclusion: these gold standards cap even real binding data at about 0.6-0.64, so the motif score (0.55-0.56) is not far below what binding data itself delivers; the limit is the endpoint, not only the score.

## 8. Corrections made along the way
- GSE53399 is a microarray (primary human hepatocytes, GW7647), not ChIP-seq. Earlier search summary was wrong. Not downloaded.
- ENCODE heart ATAC IDR file has duplicated lines: 26,918 lines = 9,717 merged peaks, not 30 Mb.
- A bug in the first ChIP-calibrated gene score (max over tiles took the global max) was fixed before the reported numbers.

## 9. Limits of the evidence
- No cardiomyocyte gold standard; PPARgene is liver/adipose dominated; heart gold sets are small and include secondary responders.
- AUROC against "other genes" is capped because most PPAR binding is distal and many negatives are unvalidated true targets.
- Many variants were tried on the same gold sets; intervals are not corrected for this.
- Promoter-peak and literature gold sets are biased towards open, active, well-studied promoters; expression alone gives 0.55-0.63 on several of them.

## 10. Defensible claim for the manuscript
"The 42-motif atlas and all FIMO hits are released. A promoter score based on PPAR-containing motifs enriches modestly (AUROC about 0.55-0.65) for validated PPARG-type promoter binding, PPARD-bound promoters and literature PPRE genes in non-cardiac cells, with no score separable from the others at current sample sizes. It is a prioritisation heuristic and is not validated in cardiomyocytes."

## 11. Open decisions
1. Switch the figure and tables to the PPAR-containing 14-motif score as primary (the figure still shows the median).
2. Add CIs to the promoter-peak table; write the claim into the README.
3. Ideas to raise AUROC (ranked): (a) filter hits by shuffled-motif FDR or p <= 1e-5/1e-6; (b) distal ChIP-calibrated model with distance decay; (c) combined motif + expression score with cross-validation; (d) conservation and strong-locus counts; (e) train on PPARD/PPARA human ChIP; (f) sequence model (CNN/gkm-SVM); (g) cardiac co-factor motifs. Evaluation: also match negatives on expression.
4. Cardiac anchors still missing: Hocker snATAC cell-type files (GEO GSE165837, 403 from this machine, download in browser); more inclusive heart ATAC peak set (ENCODE pseudoreplicated file ENCFF544RZR).
5. Table S2-S4 and the PPRE_full tab in the Google Sheet still show the old 22-motif run; alignment anchoring for RXR IR0 motifs undecided.

## Key scripts
`13` score variants and GO validation; `14` confounders; `15` 16 alternative scores; `16a/16b` binning and distal kernels; `17` hit counts vs score; `18` cardiac-context ChIP; `19` PPARG calibration; `20` ChIP-calibrated gene score; `21` promoter-peak gold; `22` motif subsets; `23` peak-level by TF; `24` heart ATAC; `25` literature PPREs; `26` PPARD myofibroblast.
