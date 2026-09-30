#!/bin/bash
# FIMO on 42 human-source PPAR/RXR motifs, unmasked GRCh38, standard p-value threshold, no practical score cap.
# Run inside WSL from ~/PPRE_FIMO. Same settings as FIMO_results_no_max (unmasked):
#   --thresh 1e-4 (FIMO default) --max-stored-scores 20000000 --bgfile genome_bg_order0.txt
set -u
cd ~/PPRE_FIMO
OUT=FIMO_results_human42_unmasked
mkdir -p $OUT/logs
run_one () {
  n="$1"
  ~/meme/bin/fimo --thresh 1e-4 --max-stored-scores 20000000 --bgfile genome_bg_order0.txt \
    --oc "$OUT/${n#RXR_PPAR_human42_files_}" "human42_batches/$n.meme" \
    Homo_sapiens.GRCh38.dna.primary_assembly.fa > "$OUT/logs/$n.log" 2>&1
  echo "$n exit=$?" >> $OUT/logs/_done.txt
}
export -f run_one; export OUT
cat human42_batches/batches.txt | xargs -P 4 -I{} bash -c 'run_one {}'
echo ALL_DONE >> $OUT/logs/_done.txt
