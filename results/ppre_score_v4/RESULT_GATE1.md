# Gate 1 first result (2026-10-08): baselines on human hepatocyte agonist responses

Code `gate1_baselines.py`; tables `gate1_baselines_frac0.05.tsv` (top 5% up) and `gate1_baselines_frac0.01.tsv` (top 1% up). Labels: paired t per gene (agonist vs vehicle, same donor), up = top fraction. GSE17251: Wy14643 vs DMSO, 24 h, 6 donors. GSE53399: GW7647 10 uM vs vehicle, 24 h, 3 donors. 16,604 genes. Exploratory development sets.

## Results (AUROC, 95% CI)
Top 5%: every score is 0.49-0.535; v3 0.522 / 0.521 (CI lower bounds 0.502 / 0.500); original sum 0.514 / 0.522; GC only 0.490 / 0.519; heart ATAC only 0.530 / 0.513; v3 + ATAC 0.535 / 0.524.
Top 1%: v3 0.536 (0.49-0.58) / 0.566 (0.52-0.61); original sum 0.560 / 0.529; hit count 0.553 / 0.522; promoter length 0.550 / 0.514; v3 + ATAC 0.533 / 0.588.

## Reading
* Promoter-window motif scores barely predict ligand response in hepatocytes. No score is consistently best; the ordering of v3 and the original sum flips between the two datasets, so v3 does not win here.
* Heart ATAC is the wrong tissue for liver cells, so the "v3 + ATAC" test is not a fair test of the accessibility idea; liver ATAC/DNase would be needed.
* Labels from 3-6 donors are noisy; AUCs of 0.52-0.59 cannot separate scores.
* GSE115827 (SGBS adipocytes) is not usable: its control is undifferentiated D0, so rosiglitazone vs control is confounded with differentiation.
* Consequence for design: only two usable development systems (both hepatocyte), so leave-one-system-out is not possible; the success rule has to be rewritten before more modelling.

## Step 2 (2026-10-08): wider window and tissue-matched ATAC (`v4_windows.py`, `v4_windows_results.tsv`)
Variants fixed in the script header before scoring. Liver ATAC ENCFF488BRH (adult liver, 136,959 merged peaks); gate = hit weight 1 inside a peak, 0.25 outside.
AUROC (GSE17251 top 5% / GSE53399 top 5% / GSE17251 top 1% / GSE53399 top 1% / literature genes):
* v3 promoter: 0.522 / 0.520 / 0.540 / 0.571 / 0.598
* **promoter + liver ATAC gate: 0.547 / 0.528 / 0.572 / 0.595 / 0.582**
* 10 kb: 0.512 / 0.522 / 0.528 / 0.567 / 0.580; 10 kb + ATAC: 0.540 / 0.530 / 0.565 / 0.592 / 0.577
* 50 kb: 0.488 / 0.520 / 0.504 / 0.546 / 0.534; 50 kb + ATAC: 0.520 / 0.531 / 0.546 / 0.575 / 0.553
* liver ATAC at the promoter alone: 0.518 / 0.516 / 0.546 / 0.528 / 0.537; original sum 0.514 / 0.522 / 0.560 / 0.529 / 0.566
Reading: the ATAC gate improves on v3 in 4 of 4 response contrasts (+0.008 to +0.032), literature genes unchanged within noise (0.582 vs 0.598); a wider window never helps and 50 kb is worse. Gains are small and the CIs overlap; not adjusted for expression, so part of the gain may be "open promoter = expressed gene". Tissue-matched ATAC is required (heart ATAC on liver cells gave no gain).
Provisional v4 = v3 site model + tissue-matched ATAC gate, promoter window only.

## Step 3 (2026-10-08): expression adjustment (`v4_expr_adjust.py`, `v4_expr_adjust.tsv`)
**The ATAC-gate gain does not survive.** Adjusting for promoter length, GC, #TSS and mean expression, adjusted OR per SD: GSE17251 top 5% v3 1.09 vs gate 1.12; top 1% 1.16 vs 1.12; GSE53399 top 5% 1.13 vs 1.12; top 1% 1.35 vs 1.36 (all CIs overlap). In expressed genes only, the paired AUC difference (gate - v3) is +0.006 (-0.007 to 0.020), -0.025 (-0.050 to -0.004), -0.011 (-0.027 to 0.003), +0.004 (-0.022 to 0.029). The step-2 gain was an expression effect (open promoter = expressed gene), not extra PPAR information.
v3 itself keeps adjusted OR 1.09-1.35 (lower CI > 1 in 3 of 4 contrasts) after expression adjustment. **Conclusion: no tested v4 change (extended window, ATAC gate) improves on v3; heart v4 not frozen.**
