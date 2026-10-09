# Pipeline for the 42-motif list, part 1 (per-motif FIMO summary, ChIPseeker promoter annotation, HCM/PLN H3K27ac overlap).
# Equivalent to the chunks load-data, annotate-fimo, summary-table and overlap-regions of Objective1_chipseeker_fimo_analysis_with_CISBP.Rmd,
# but one motif at a time (45.5 million hits do not fit in one GRanges). ChIPseeker annotates each peak independently of the other peaks,
# so annotating only hits near a TSS gives the same "Promoter" rows as annotating all hits.
suppressPackageStartupMessages({library(data.table); library(GenomicRanges); library(ChIPseeker); library(TxDb.Hsapiens.UCSC.hg38.knownGene); library(org.Hs.eg.db); library(rtracklayer); library(readxl); library(GenomeInfoDb)})
args <- commandArgs(trailingOnly = TRUE); only_dir <- if (length(args)) args[1] else NA
base <- "C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Minor internship"
fimo_root <- file.path(base, "PPAR_RXR_motifs_Hocomoco_Jaspar", "FIMO_results_human42_unmasked")
res_root <- file.path(base, "PPAR_RXR_motifs_Hocomoco_Jaspar", "Human42_analysis")
chip_dir <- file.path(res_root, "ChIPseeker_output"); res_dir <- file.path(res_root, "Objective1_results"); rng_dir <- file.path(res_root, "hit_ranges_rds")
for (d in c(chip_dir, res_dir, rng_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
prom_zone <- reduce(suppressWarnings(trim(promoters(txdb, upstream = 3100, downstream = 1100))), ignore.strand = TRUE); prom_zone <- keepStandardChromosomes(prom_zone, pruning.mode = "coarse")
# H3K27ac regions (hg19) lifted to hg38, as in the notebook
HCM <- read_xlsx(file.path(base, "PPAR_RXR_motifs_Hocomoco_Jaspar", "HCM_differentially_regions.xlsx"), sheet = "Supplementary Table 2A")
PLN <- read_xlsx(file.path(base, "PPAR_RXR_motifs_Hocomoco_Jaspar", "PLN_differentially_regions.xlsx"), skip = 1)
ci <- function(d, pat) grep(pat, names(d))[1]                       # column lookup by prefix (names carry footnote symbols)
PLN <- PLN[!is.na(PLN[[ci(PLN, "^Region start")]]), ]
g_hcm <- GRanges(HCM[["Chromosome number"]], IRanges(as.integer(HCM[["Region start"]]), as.integer(HCM[["Region end"]])))
g_pln <- GRanges(PLN[[ci(PLN, "^Chromosome")]], IRanges(as.integer(PLN[[ci(PLN, "^Region start")]]), as.integer(PLN[[ci(PLN, "^Region end")]])))
chain <- import.chain(file.path(base, "PPAR_RXR_motifs_Hocomoco_Jaspar", "hg19ToHg38.over.chain"))
g_hcm <- unlist(liftOver(g_hcm, chain)); g_pln <- unlist(liftOver(g_pln, chain))
cat("H3K27ac regions after liftOver: HCM", length(g_hcm), "PLN", length(g_pln), "\n")
dirs <- list.dirs(fimo_root, recursive = FALSE); if (!is.na(only_dir)) dirs <- dirs[basename(dirs) == only_dir]
summ <- list(); ov <- list(); t0 <- Sys.time()
for (dd in dirs) {
  f <- file.path(dd, "fimo.tsv"); cat("\n==", basename(dd), format(Sys.time()), "\n")
  x <- fread(f, sep = "\t", header = TRUE, fill = TRUE, colClasses = "character", select = c("motif_id", "sequence_name", "start", "stop", "strand", "score", "p-value", "q-value", "matched_sequence"))
  x <- x[!is.na(start) & start != "" & !startsWith(motif_id, "#")]; x[, `:=`(start = as.integer(start), stop = as.integer(stop), score = as.numeric(score), pv = as.numeric(`p-value`), qv = as.numeric(`q-value`))]
  for (m in unique(x$motif_id)) {
    xm <- x[motif_id == m]; gr <- GRanges(xm$sequence_name, IRanges(xm$start, xm$stop), strand = xm$strand, score = xm$score, pvalue = xm$pv, qvalue = xm$qv, motif_id = xm$motif_id, motif_seq = xm$matched_sequence)
    suppressWarnings(seqlevelsStyle(gr) <- "UCSC"); gr <- keepStandardChromosomes(gr, pruning.mode = "coarse")
    saveRDS(granges(gr, use.mcols = FALSE), file.path(rng_dir, paste0(m, ".rds")))
    n_all <- nrow(xm); minp <- min(xm$pv); mw <- mean(xm$stop - xm$start + 1)
    o_h <- length(queryHits(findOverlaps(g_hcm, gr, ignore.strand = TRUE, minoverlap = 1))); o_p <- length(queryHits(findOverlaps(g_pln, gr, ignore.strand = TRUE, minoverlap = 1)))
    near <- gr[overlapsAny(gr, prom_zone, ignore.strand = TRUE)]
    pk <- annotatePeak(near, tssRegion = c(-3000, 1000), TxDb = txdb, annoDb = "org.Hs.eg.db", sameStrand = TRUE, verbose = FALSE)
    pdf_ <- as.data.frame(pk); prom <- pdf_[grepl("Promoter", pdf_$annotation), ]
    write.table(prom, file.path(chip_dir, paste0(m, ".txt")), sep = " ")
    summ[[m]] <- data.frame(N_amount = n_all, Min_pvalue = minp, Mean_width = mw, Motif_name = m, Motif_seq = xm$matched_sequence[1], N_amount_genes = nrow(prom), Mean_dist_tss = mean(prom$distanceToTSS), n_near_tss = length(near))
    ov[[m]] <- data.frame(motif_id = m, total_hits = n_all, overlap_hcm = o_h, overlap_pln = o_p)
    cat(sprintf("%-24s hits %9s | near TSS %8s | promoter rows %8s | HCM %d PLN %d | %.1f min\n", m, format(n_all, big.mark = ","), format(length(near), big.mark = ","), format(nrow(prom), big.mark = ","), o_h, o_p, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
  rm(x); gc(verbose = FALSE)
}
S <- do.call(rbind, summ); O <- do.call(rbind, ov)
suffix <- if (is.na(only_dir)) "" else paste0("_", only_dir)
write.csv(S, file.path(res_dir, paste0("motif_summary_chipseeker_42", suffix, ".csv")), row.names = FALSE); write.csv(O, file.path(res_dir, paste0("H3K27ac_overlap_42", suffix, ".csv")), row.names = FALSE)
cat("done; total", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
