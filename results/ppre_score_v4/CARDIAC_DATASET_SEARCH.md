# Search for new human cardiac PPAR-response data (2026-10-08)

Method: NCBI GEO (E-utilities, db=gds) with 5 query families: (1) PPAR agonist names x cardiomyocyte/cardiac/heart (40 records), (2) iPSC/hPSC cardiomyocytes x agonist/fatty acid/PPAR (32), (3) PPARA/D/G, RXRA, PGC-1a knockdown/knockout/overexpression x cardiac (3), (4) iPSC-cardiomyocyte drug and chemical panels (91), (5) PPAR in fibroblast/endothelial/myocardium/heart failure (33); about 70 record summaries read, plus supplement lists and ENA assay types for the leads. Web search for iPSC-CM PPAR agonist transcriptome. Not searched: ArrayExpress, SRA free text, LINCS, CMap. Nothing was downloaded except metadata.

## Candidates (human only)
| Rank | Accession | What it is | Use | Concern |
|---|---|---|---|---|
| 1 | GSE295444, GSE295447 | 12 iPSC-derived cell types x 12 compounds, 7 h, 3 donors, 3 replicates, DMSO controls; includes **cardiomyocytes with rosiglitazone 1 uM** | independent cardiac PPARG-agonist contrast | rosiglitazone gave almost no response in GSE160987 cardiomyocytes (1 gene at padj < 0.05); likely weak here too; 1.7 GB tars, but per-sample Salmon files are about 2 MB each, so only the needed samples would be fetched |
| 2 | GSE262419 | iPSC-cardiomyocytes, 464 chemicals, concentration-response, 16 plates (TempO-seq-style), 1,152 samples | independent cardiac set if PPAR agonists (fenofibrate, glitazones, PFAS) are among the 464 | chemical list not visible in GEO metadata; the 47 KB `hash.csv.gz` likely maps wells to chemicals |
| 3 | GSE125862 | hiPSC-cardiomyocytes treated with FM19G11 + **WY14643** + T3/IGF-1/dexamethasone vs DMSO (3 v 3), plus 3 adult left ventricle samples | PPARA-agonist-containing cardiac contrast | combination treatment, so response mixes maturation signalling with PPARA; not a clean PPAR contrast |
| 4 | GSE178984 | **ATAC-seq of the same hPSC-cardiomyocytes as GSE160987** (Control, LCFA, GW0742, GSK0660 +/- LCFA, 3 replicates each; 2.9 MB normalised counts) | cardiomyocyte-specific accessibility layer; GW0742-induced regions as a cardiac binding proxy | same cells as GSE160987, so not independent of it |
| 5 | GSE218902 | cosmetic-ingredient panel in MCF7, HepG2, A549 and cardiomyocytes, 420 samples | possible if PPAR-active compounds (parabens) are included | compounds weakly PPAR-active; not checked |
| - | GSE83668 | iPSC-CM rosiglitazone vs untreated, AmpliSeq | none | one sample per donor, no replicates |
| - | GSE124057, GSE298270 | fatty acids (4 samples); CRISPRa PPARGC1A in AC16 (4 samples) | none | too small / coactivator not receptor |

Not relevant: ERRg (GSE135319, GSE113760), RAR agonist Am80 (GSE176303), statins (GSE113546), KDM5 (GSE250210), PPARG knockdown in valve endothelial cells (GSE206927, not cardiomyocyte), non-human sets.

## Reading
No clean, independent human cardiomyocyte PPAR-response set with a strong effect was found beyond the GW0742 contrast of GSE160987 (already used). The best chances are GSE262419 (if PPAR agonists are in the chemical list) and the cardiomyocyte rosiglitazone contrast in GSE295444/GSE295447 (probably weak). A sanity gate must be fixed before scoring any of them: known PPAR targets (PDK4, ANGPTL4, CPT1A, CPT1B) induced and a minimum number of differentially expressed genes; otherwise the set is reported as uninformative.

## Update 2026-10-08: GSE262419 downloaded and checked
Plates 8, 9, 10, 12, 14 (about 27 MB). Pirinixic acid (WY-14643) 10 uM: 10 genes at padj < 0.05, no PPAR target induced; troglitazone: 0-1 genes; PFOA, PFNA, PFDA, PFOS: 6-117 genes, no target induced. Sanity gate failed in all 8 contrasts (see DESIGN_V4.md Addendum 8), so GSE262419 is uninformative and was not scored. Remaining leads: GSE295444/GSE295447 (cardiomyocyte rosiglitazone, expected weak), GSE125862 (combination with WY-14643), GSE178984 (cardiomyocyte ATAC, same cells as GSE160987).

## Update 2026-10-08 (2): GSE178984 ATAC checked
GW0742 effect on accessibility: 1 of 35,383 peaks at FDR < 0.05 (gate: 100), so it cannot serve as a binding proxy. The PDK4 enhancer opens (log2FC +2.9, FDR 0.17). Not scored (DESIGN_V4.md Addendum 9).
## Conclusion of the search
No dataset with an informative human cardiomyocyte PPAR response beyond the GW0742 RNA-seq contrast of GSE160987 (already used) was found. Untested leads: GSE295444/GSE295447 (cardiomyocyte rosiglitazone; rosiglitazone gave no response in GSE160987 or GSE262419-type conditions) and GSE125862 (WY-14643 combined with other factors).
