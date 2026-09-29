# Step 6: build a CentriMo FASTA from the ChIP-Atlas peaks restricted to liver-context
# PPARA ChIP-seq experiments (7 SRX experiments, dominated by the HepG2 series that
# includes GSM4748812/816/817, plus one more liver SRX). Peaks from the different
# experiments overlap heavily, so they're reduce()'d to non-redundant loci first, each
# then centered on its midpoint +/-250bp to match the GSM4748812 CentriMo run.

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(BSgenome.Hsapiens.UCSC.hg38)
})

out_dir <- "PPARA_chipseq_overlap"
raw <- fread(file.path(out_dir, "chipatlas_liver_peaks_raw.bed"), header = FALSE,
             col.names = c("chrom", "start", "end"))
cat("Raw liver-tagged peak records:", nrow(raw), "\n")

genome <- BSgenome.Hsapiens.UCSC.hg38
raw <- raw[chrom %in% seqnames(genome)]

gr <- GRanges(raw$chrom, IRanges(raw$start + 1, raw$end))
gr_reduced <- reduce(gr)
cat("Non-redundant merged loci:", length(gr_reduced), "\n")

mid <- (start(gr_reduced) + end(gr_reduced)) %/% 2L
half_width <- 250L
win <- GRanges(seqnames(gr_reduced), IRanges(mid - half_width, mid + half_width))
win <- win[start(win) > 0 & end(win) <= seqlengths(genome)[as.character(seqnames(win))]]

seqs <- getSeq(genome, win)
names(seqs) <- sprintf("%s:%d-%d", seqnames(win), start(win), end(win))

fasta_path <- file.path(out_dir, "ChIPAtlas_liver_peaks_centered500bp.fa")
Biostrings::writeXStringSet(seqs, fasta_path)
cat("Wrote", length(seqs), "sequences to", fasta_path, "\n")
