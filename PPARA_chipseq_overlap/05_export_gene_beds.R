# Step 5: write one promoter-FIMO-hit BED per gene, same 6-column format as the existing
# CPT1B_promoter_motif_hits.bed, plus two extra columns flagging whether each hit overlaps
# the single-experiment (GSM4748812) or aggregate (ChIP-Atlas) PPARA ChIP-seq peak sets.

suppressPackageStartupMessages(library(data.table))

out_dir <- "PPARA_chipseq_overlap"
fimo <- fread(file.path(out_dir, "promoter_fimo_hits_annotated.csv"))

for (g in sort(unique(fimo$gene))) {
  sub <- fimo[gene == g, .(chrom, start, end = stop, motif_id, score,
                            strand, in_GSM4748812, in_ChIPAtlas_any, in_ChIPAtlas_ge3)]
  setorder(sub, start)
  fwrite(sub, file.path(out_dir, paste0(g, "_promoter_motif_hits_with_chipseq.bed")),
         sep = "\t", col.names = FALSE)
}

cat("Per-gene BED files written to", out_dir, "\n")
