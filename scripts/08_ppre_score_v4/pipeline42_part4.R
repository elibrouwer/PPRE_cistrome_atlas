# Part 4: genome coverage and number of merged PPRE regions for the 42 motifs (union of all hits, strand ignored).
suppressPackageStartupMessages({library(GenomicRanges); library(BSgenome.Hsapiens.UCSC.hg38)})
base <- "C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Minor internship"
rng_dir <- file.path(base, "PPAR_RXR_motifs_Hocomoco_Jaspar", "Human42_analysis", "hit_ranges_rds")
files <- list.files(rng_dir, pattern = "[.]rds$", full.names = TRUE)
gl <- lapply(files, function(f) { g <- readRDS(f); strand(g) <- "*"; g }); allr <- do.call(c, gl); rm(gl)
std <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)[paste0("chr", c(1:22, "X", "Y"))]; gsize <- sum(as.numeric(std))
red <- reduce(allr); cat(sprintf("hits %s | merged regions %s | covered bp %s | genome %s bp | coverage %.2f%%\n", format(length(allr), big.mark = ","), format(length(red), big.mark = ","), format(sum(as.numeric(width(red))), big.mark = ","), format(gsize, big.mark = ","), sum(as.numeric(width(red))) / gsize * 100))
cat(sprintf("merged region width: median %d, mean %.1f\n", median(width(red)), mean(width(red))))
