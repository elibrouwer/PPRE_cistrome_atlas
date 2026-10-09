# Rebuild the per-motif summary and the H3K27ac overlap table from the files written by part 1 (the run ended on an empty logs folder before saving them).
suppressPackageStartupMessages({library(GenomicRanges); library(rtracklayer); library(readxl); library(GenomeInfoDb)})
base <- "C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Minor internship"; B <- file.path(base, "PPAR_RXR_motifs_Hocomoco_Jaspar")
res_root <- file.path(B, "Human42_analysis"); rng_dir <- file.path(res_root, "hit_ranges_rds"); chip_dir <- file.path(res_root, "ChIPseeker_output"); res_dir <- file.path(res_root, "Objective1_results")
HCM <- read_xlsx(file.path(B, "HCM_differentially_regions.xlsx"), sheet = "Supplementary Table 2A"); PLN <- read_xlsx(file.path(B, "PLN_differentially_regions.xlsx"), skip = 1)
ci <- function(d, pat) grep(pat, names(d))[1]; PLN <- PLN[!is.na(PLN[[ci(PLN, "^Region start")]]), ]
g_hcm <- GRanges(HCM[["Chromosome number"]], IRanges(as.integer(HCM[["Region start"]]), as.integer(HCM[["Region end"]])))
g_pln <- GRanges(PLN[[ci(PLN, "^Chromosome")]], IRanges(as.integer(PLN[[ci(PLN, "^Region start")]]), as.integer(PLN[[ci(PLN, "^Region end")]])))
chain <- import.chain(file.path(B, "hg19ToHg38.over.chain")); g_hcm <- unlist(liftOver(g_hcm, chain)); g_pln <- unlist(liftOver(g_pln, chain))
S <- list(); O <- list()
for (f in list.files(rng_dir, pattern = "[.]rds$", full.names = TRUE)) {
  m <- sub("[.]rds$", "", basename(f)); gr <- readRDS(f); prom <- read.delim(file.path(chip_dir, paste0(m, ".txt")), sep = " ", header = TRUE)
  S[[m]] <- data.frame(N_amount = length(gr), Mean_width = mean(width(gr)), Motif_name = m, N_amount_genes = nrow(prom), Mean_dist_tss = mean(prom$distanceToTSS))
  oh <- findOverlaps(g_hcm, gr, ignore.strand = TRUE, minoverlap = 1); op <- findOverlaps(g_pln, gr, ignore.strand = TRUE, minoverlap = 1)
  O[[m]] <- data.frame(motif_id = m, total_hits = length(gr), overlap_hcm = length(queryHits(oh)), overlap_pln = length(queryHits(op)), regions_hcm_hit = length(unique(queryHits(oh))), regions_pln_hit = length(unique(queryHits(op))), regions_hcm_total = length(g_hcm), regions_pln_total = length(g_pln))
  cat(m, length(gr), nrow(prom), "\n")
}
write.csv(do.call(rbind, S), file.path(res_dir, "motif_summary_chipseeker_42.csv"), row.names = FALSE); write.csv(do.call(rbind, O), file.path(res_dir, "H3K27ac_overlap_42.csv"), row.names = FALSE); cat("saved\n")
