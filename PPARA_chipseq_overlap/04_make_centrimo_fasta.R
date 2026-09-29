# Step 4: build the peak-centered FASTA CentriMo needs. Uses the FULL genome-wide
# GSM4748812 (HepG2, PPARA-only condition) narrowPeak set - not just the 5 gene
# promoters - because CentriMo needs many peaks to get a reliable positional-enrichment
# signal. Peaks are recentered on their MACS2 summit +/- 250bp so every sequence has the
# same width and CentriMo's central-enrichment test is meaningful.

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(BSgenome.Hsapiens.UCSC.hg38)
})

out_dir <- "PPARA_chipseq_overlap"
peaks_path <- "data/downloads/GSM4748812_1_S1_peaks.narrowPeak.gz"

np_cols <- c("chrom","start","end","name","score","strand","signalValue","pValue","qValue","summit")
peaks <- fread(peaks_path, header = FALSE, col.names = np_cols)
cat("Loaded", nrow(peaks), "peaks\n")

# narrowPeak is 0-based half-open; summit is an OFFSET from start (0-based)
peaks[, summit_pos := start + summit + 1L]  # convert to 1-based genomic coordinate

half_width <- 250L
peaks[, win_start := summit_pos - half_width]
peaks[, win_end   := summit_pos + half_width]

genome <- BSgenome.Hsapiens.UCSC.hg38
valid_chroms <- seqnames(genome)
peaks <- peaks[chrom %in% valid_chroms]
peaks <- peaks[win_start > 0]

gr <- GRanges(peaks$chrom, IRanges(peaks$win_start, peaks$win_end))
seqs <- getSeq(genome, gr)
names(seqs) <- sprintf("%s:%d-%d_%s", peaks$chrom, peaks$win_start, peaks$win_end, peaks$name)

fasta_path <- file.path(out_dir, "GSM4748812_peaks_centered500bp.fa")
Biostrings::writeXStringSet(seqs, fasta_path)
cat("Wrote", length(seqs), "sequences to", fasta_path, "\n")
