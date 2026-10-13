#!/bin/bash
set -u
cd ~/ctrl42; mkdir -p out logs; rm -f logs/_done.txt
run_one(){ n="$1"; ~/meme/bin/fimo --thresh 1e-4 --max-stored-scores 50000000 --bgfile ~/PPRE_FIMO/genome_bg_order0.txt --oc out/$n memes/$n.meme promoter_union.fa > logs/$n.log 2>&1; echo "$n exit=$?" >> logs/_done.txt; }
export -f run_one
printf "real\nperm01\nperm02\nperm03\nperm04\nperm05\nrev\n" | xargs -P 7 -I{} bash -c 'run_one {}'
echo ALL_DONE >> logs/_done.txt
