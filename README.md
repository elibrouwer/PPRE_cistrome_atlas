# PPRE cistrome atlas

Code and results for the manuscript *A comprehensive atlas of the predicted human PPRE cistrome decodes transcription factor target heterogeneity in cardiomyocytes*.

The single-nucleus RNA-seq analysis of the cardiomyocyte subpopulations is in [PPRE_cardiomyocytes_manuscript](https://github.com/elibrouwer/PPRE_cardiomyocytes_manuscript).

## Approach

1. 30 PPAR/RXR motifs (HOCOMOCO, JASPAR and CIS-BP) are aligned and merged into consensus logos.
2. The motifs are scanned genome-wide with FIMO (MEME Suite 5.5.7, GRCh38, p <= 1e-4, order-0 genome background).
3. The hits are annotated with ChIPseeker (3 kb upstream to 1 kb downstream of the TSS) and overlapped with genes expressed in GTEx cardiomyocytes.
4. Every gene gets a PPRE score, which is the FIMO score weighted by an exponential decay of the distance to the TSS. The sum, mean, median and maximum per gene are calculated.

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
