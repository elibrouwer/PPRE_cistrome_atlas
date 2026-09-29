# PPRE cistrome atlas

Genome-wide identification and scoring of PPAR::RXR response elements (PPREs)
and their overlap with cardiomyocyte-expressed genes, supporting the PPRE
cardiomyocytes manuscript. Companion repository (single-cell / MERSCOPE
analysis, Python): [PPRE_cardiomyocytes_manuscript](https://github.com/elibrouwer/PPRE_cardiomyocytes_manuscript).

## Approach

1. Publicly available PPAR and RXR motifs (HOCOMOCO, JASPAR, plus CIS-BP
   additions) are aligned/merged into a pan-PPRE model (`results/pan_ppre/`).
2. Motifs are scanned genome-wide with FIMO (MEME Suite 5.5.7, GRCh38,
   order-0 genomic background, p <= 1e-4).
3. Hits are annotated to promoters (ChIPseeker) and intersected with genes
   expressed in GTEx cardiomyocyte snRNA-seq and with PamGene data.
4. Hits are validated against ChIP-seq (PPARA, PPARD) and by threshold
   sensitivity analysis.

## Layout

| Path | Contents |
|---|---|
| `scripts/` | R / Rmd analysis scripts (`Objective1_*` motif analysis, `Objective3_*` multi-omics comparison, GTEx/PamGene scripts) |
| `reports/` | Rendered HTML of the Rmd notebooks |
| `PPRE_threshold_sensitivity/` | Numbered pipeline (`00_config.R` ... `24_*`) for FIMO threshold sensitivity; see its own README |
| `PPARA_chipseq_overlap/` | PPARA ChIP-seq overlap and CentriMo analysis (`01_` ... `06_`) |
| `results/figures`, `results/tables`, `results/beds`, `results/pan_ppre` | Key outputs used in the manuscript |

## Reproducing

Scripts expect raw inputs under a `data/` directory (git-ignored), mirroring:

- `data/PPAR_RXR_motifs_Hocomoco_Jaspar/` : motif files, FIMO output
  (`FIMO_results_no_max`, `FIMO_results_genome_masked`), `CIS-BP/`, `Objective1_results/`
- `data/GTEX_single_nuclear_data/` : GTEx snRNA-seq expression tables
- `data/Pam_chip/`, `data/Objective3_MultiOmics_overlaps_PPRE/` : PamGene data
- `data/downloads/` : ChIP-seq peak files (GSM4748812; ChIP-Atlas `Oth.ALL.05.PPARA.AllCell.bed`)

Run scripts from the repository root. The 214 MB `ChIPAtlas_liver_peaks_centered500bp.fa`
and the raw ChIP-Atlas bed are not included; regenerate them with
`PPARA_chipseq_overlap/06_make_liver_centrimo_fasta.R`.

R packages include ChIPseeker, GenomicRanges, rtracklayer, ggplot2, ComplexHeatmap
and their dependencies; FIMO/CentriMo come from the MEME Suite.
