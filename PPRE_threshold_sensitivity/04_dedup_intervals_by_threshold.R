# For each threshold tier (p<=1e-4, 1e-5, 1e-6):
#   (a) raw motif-call counts, per motif and per family (JASPAR/HOCOMOCO/CIS-BP
#       overlapping predictions at the same locus are NOT yet collapsed here --
#       this is the "raw motif calls" table)
#   (b) deduplicated genomic intervals: overlapping calls from different
#       motifs/databases at the same locus are merged into one interval,
#       retaining which motif(s) and family(ies) support each interval --
#       this is the "unique genomic intervals" table
#
# Full heterodimer, PPAR half-site, and RXR half-site hits are kept separate
# throughout (both in the raw and deduplicated views) before any combined
# view is produced downstream.

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(dplyr)
})
source("./PPRE_threshold_sensitivity/00_config.R")

all_hits <- readRDS(file.path(out_data, "all_hits_p1e-4.rds"))
motif_family <- read.csv(motif_family_csv, stringsAsFactors = FALSE)

all_hits <- merge(all_hits, motif_family[, c("motif_id", "family")], by = "motif_id", all.x = TRUE)
stopifnot(!anyNA(all_hits$family))

raw_summary_list   <- list()
interval_list      <- list()
interval_summary_list <- list()

for (pt in p_thresholds) {
  tag <- paste0("p", format(pt, scientific = TRUE))
  cat(sprintf("\n=== threshold %s ===\n", tag))

  sub <- all_hits[pvalue <= pt]
  cat(sprintf("Raw hits kept: %s\n", format(nrow(sub), big.mark = ",")))

  ## ---- (a) raw hits per motif and per family ----
  raw_by_motif <- sub[, .(n_raw_hits = .N, family = family[1]), by = motif_id]
  raw_by_family <- sub[, .(n_raw_hits = .N, n_motifs_contributing = uniqueN(motif_id)), by = family]
  raw_by_motif[, threshold := pt]; raw_by_family[, threshold := pt]
  raw_summary_list[[tag]] <- list(by_motif = raw_by_motif, by_family = raw_by_family)

  ## ---- (b) deduplicated genomic intervals (merge overlapping calls, ignore strand) ----
  seqn <- ifelse(sub$sequence_name %in% c("MT", "M"), "chrM", paste0("chr", sub$sequence_name))
  gr <- GRanges(seqnames = seqn, ranges = IRanges(start = sub$start, end = sub$stop))
  std_chroms <- paste0("chr", c(1:22, "X", "Y", "M"))
  keep <- as.character(seqnames(gr)) %in% std_chroms
  gr <- gr[keep]; sub_std <- sub[keep]
  cat(sprintf("Hits on non-standard contigs dropped: %s\n", format(sum(!keep), big.mark = ",")))

  merged <- reduce(gr, ignore.strand = TRUE)
  hits_ov <- findOverlaps(gr, merged, ignore.strand = TRUE)
  sub_std[, interval_id := subjectHits(hits_ov)]

  interval_agg <- sub_std[, .(
    n_raw_hits           = .N,
    n_distinct_motifs    = uniqueN(motif_id),
    motif_ids            = paste(sort(unique(motif_id)), collapse = ";"),
    families_present     = paste(sort(unique(family)), collapse = ";"),
    n_heterodimer_hits   = sum(family == "heterodimer"),
    n_PPAR_half_hits     = sum(family == "PPAR_half"),
    n_RXR_half_hits      = sum(family == "RXR_half"),
    best_pvalue          = min(pvalue),
    best_score           = max(score)
  ), by = interval_id]

  interval_coords <- data.table(
    interval_id = seq_along(merged),
    chrom = as.character(seqnames(merged)),
    start = start(merged),
    end   = end(merged),
    width = width(merged)
  )
  interval_tbl <- merge(interval_coords, interval_agg, by = "interval_id")
  interval_tbl[, threshold := pt]

  cat(sprintf("Deduplicated intervals: %s (from %s raw hits, %.2fx collapse)\n",
              format(nrow(interval_tbl), big.mark = ","),
              format(nrow(sub_std), big.mark = ","),
              nrow(sub_std) / nrow(interval_tbl)))

  interval_list[[tag]] <- interval_tbl

  interval_summary_list[[tag]] <- interval_tbl[, .(
    n_intervals = .N,
    n_intervals_heterodimer_only = sum(families_present == "heterodimer"),
    n_intervals_PPAR_half_only   = sum(families_present == "PPAR_half"),
    n_intervals_RXR_half_only    = sum(families_present == "RXR_half"),
    n_intervals_multi_family     = sum(grepl(";", families_present))
  )][, threshold := pt]

  fwrite(interval_tbl, file.path(out_results, paste0("deduplicated_intervals_", tag, ".csv")))
}

raw_by_motif_all  <- rbindlist(lapply(raw_summary_list, `[[`, "by_motif"))
raw_by_family_all <- rbindlist(lapply(raw_summary_list, `[[`, "by_family"))
interval_summary_all <- rbindlist(interval_summary_list)

fwrite(raw_by_motif_all,  file.path(out_results, "raw_hits_per_motif_by_threshold.csv"))
fwrite(raw_by_family_all, file.path(out_results, "raw_hits_per_family_by_threshold.csv"))
fwrite(interval_summary_all, file.path(out_results, "deduplicated_interval_summary_by_threshold.csv"))

saveRDS(interval_list, file.path(out_data, "deduplicated_intervals_by_threshold.rds"))

cat("\n=== Summary across thresholds ===\n")
print(raw_by_family_all)
print(interval_summary_all)
message("Done. Per-threshold interval tables, and summary CSVs, written to ", out_results)
