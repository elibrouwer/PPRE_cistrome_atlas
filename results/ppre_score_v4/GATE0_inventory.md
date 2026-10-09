# Gate 0: data inventory for v4 (2026-10-08)

Metadata only (NCBI GEO via E-utilities). **No data files downloaded yet.** Search: human, expression profiling, PPAR agonist terms (133 GEO records matched; the candidates below were picked for vehicle control + replicates + a clean agonist contrast).

## Functional (gene-level ligand response) candidates
| Role | Accession | System | Ligand | Platform | Samples | Notes / concerns |
|---|---|---|---|---|---|---|
| **Locked test** | GSE160987 | hPSC-cardiomyocytes | GW0742, WY14643, rosiglitazone | RNA-seq | in project (DESeq2 labels built) | evaluated once, at the end |
| Dev 1 | GSE53399 | primary human hepatocytes | GW7647 (PPARA) | array GPL13158 | 120 (dose x time x donor) | paired by donor; choose one fixed contrast (e.g. 24 h, highest dose vs vehicle) before looking at scores |
| Dev 2 | GSE17251 | human hepatocytes / liver | Wy14643 (PPARA) | array GPL570 | 24 (donor x time x DMSO/Wy) | same tissue as Dev 1, different ligand; counts as a partly dependent system |
| Dev 3 | GSE115827 | SGBS adipocytes | rosiglitazone (PPARG) | RNA-seq | 12 (incl. dual agonists, D0 control) | only rosiglitazone vs control arms used; check replicates per arm |
| Optional | GSE33152 | primary hepatocytes | dual PPARa/g agonists | array GPL6884 | 54 | alternative to Dev 2 |
| Too small | GSE50378 | HUVEC | GW501516 | array GPL570 | 4 total | no replicates per arm, cannot define response labels; dropped from the functional set (HUVEC stays a binding set) |

Honest limits: only one system is PPARG-driven (adipocyte), two are PPARA in liver, and the locked set is a cardiomyocyte mix of ligands. "Three independent cell systems" is met only as hepatocyte / adipocyte / (locked) cardiomyocyte, with hepatocyte as one system, so the leave-one-system-out rule becomes **2 development systems + 1 locked**. DESIGN_V4 section 6 criterion 1 must be amended to "2 of 2 development systems and the locked set" (deviation to be logged).

## Binding sets: already in project
PPARG adipocyte, digestive tract (HT-29), blood (THP-1), HUVEC PPARD (GSE50144), HepG2 PPARG (ENCODE), RXRA sets.

## Chromatin layers
* ENCODE heart left-ventricle ATAC ENCFF174HPZ: in project.
* phastCons 100-way: in project.
* Cardiac H3K27ac peaks: **not in project**; needs download (ENCODE heart LV H3K27ac narrowPeak, size to be confirmed, expected a few MB).

## Proposed downloads (need your yes)
1. GSE53399 series matrix (small, <10 MB expected, 120 arrays).
2. GSE17251 series matrix (<5 MB).
3. GSE115827 count/supplementary table (<5 MB).
4. ENCODE heart LV H3K27ac peak file (a few MB).
Sizes are expected, not confirmed; I will print the real sizes before fetching.
