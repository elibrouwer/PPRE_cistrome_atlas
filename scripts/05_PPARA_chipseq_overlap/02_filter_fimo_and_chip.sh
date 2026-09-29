#!/bin/bash
# Step 2: select the FIMO hits and the ChIP-seq peaks that fall in the promoter windows.
# awk is used so the large fimo.tsv files do not have to be loaded into R.

set -e

OUT_DIR="results/PPARA_chipseq_overlap"
FIMO_DIR="data/PPAR_RXR_motifs_Hocomoco_Jaspar/FIMO_results_no_max"
REGIONS="$OUT_DIR/gene_promoter_windows.bed"

GSM_PEAKS="data/downloads/GSM4748812_1_S1_peaks.narrowPeak.gz"
CHIPATLAS_PEAKS="data/downloads/Oth.ALL.05.PPARA.AllCell.bed"

BATCHES="1_4 5_8 9_12 13_16 17_20 21_22 23_26 27_30"

echo "=== Filtering genome-wide FIMO batches to promoter windows ==="
FIMO_OUT="$OUT_DIR/promoter_fimo_hits.tsv"
: > "$FIMO_OUT"
for b in $BATCHES; do
  f="$FIMO_DIR/$b/fimo.tsv"
  echo "  scanning $b ..."
  awk -v regions="$REGIONS" '
    BEGIN{
      FS="\t"; OFS="\t"; n=0
      while ((getline line < regions) > 0) {
        split(line, a, "\t")
        chrom=a[1]; sub(/^chr/,"",chrom)
        n++; rchrom[n]=chrom; rstart[n]=a[2]+1; rend[n]=a[3]; rgene[n]=a[4]
      }
      close(regions)
    }
    NR==1{next}  # skip header
    {
      for (i=1;i<=n;i++) {
        if ($3==rchrom[i] && $4<=rend[i] && $5>=rstart[i]) {
          print $0, rgene[i]
        }
      }
    }
  ' "$f" >> "$FIMO_OUT"
done
echo "  -> $(wc -l < "$FIMO_OUT") promoter-window FIMO hits written to $FIMO_OUT"

echo "=== Filtering GSM4748812 (single-experiment HepG2 PPARA) narrowPeak ==="
GSM_OUT="$OUT_DIR/promoter_GSM4748812_peaks.bed"
zcat "$GSM_PEAKS" | awk -v regions="$REGIONS" '
  BEGIN{
    FS="\t"; OFS="\t"; n=0
    while ((getline line < regions) > 0) {
      split(line, a, "\t")
      n++; rchrom[n]=a[1]; rstart[n]=a[2]; rend[n]=a[3]; rgene[n]=a[4]
    }
    close(regions)
  }
  {
    for (i=1;i<=n;i++) {
      if ($1==rchrom[i] && $2<=rend[i] && $3>=rstart[i]) {
        print $0, rgene[i]
      }
    }
  }
' > "$GSM_OUT"
echo "  -> $(wc -l < "$GSM_OUT") peaks written to $GSM_OUT"

echo "=== Filtering ChIP-Atlas all-cell-type PPARA aggregate (739 MB, streaming) ==="
CA_OUT="$OUT_DIR/promoter_chipatlas_peaks.bed"
awk -v regions="$REGIONS" '
  BEGIN{
    FS="\t"; OFS="\t"; n=0
    while ((getline line < regions) > 0) {
      split(line, a, "\t")
      n++; rchrom[n]=a[1]; rstart[n]=a[2]; rend[n]=a[3]; rgene[n]=a[4]
    }
    close(regions)
  }
  NR==1{next}  # skip track header line
  {
    for (i=1;i<=n;i++) {
      if ($1==rchrom[i] && $2<=rend[i] && $3>=rstart[i]) {
        print $0, rgene[i]
      }
    }
  }
' "$CHIPATLAS_PEAKS" > "$CA_OUT"
echo "  -> $(wc -l < "$CA_OUT") peak records written to $CA_OUT"

echo "Done."
