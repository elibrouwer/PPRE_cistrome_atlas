# Load all FIMO hits (p<=1e-4, i.e. everything FIMO_results_no_max kept) for
# all 30 motifs, across all 8 shards, into one compact table.
#
# FIMO_results_no_max/{shard}/fimo.tsv files total ~2.5GB across the 8 shards.
# Read shard-by-shard with data.table::fread (fast, low peak memory), drop
# columns not needed downstream (motif_alt_id is mostly empty; q-value is not
# used -- see Objective1_fimo_threshold_enrichment.R's note that FIMO's
# q-values are compressed too tightly across this genome-wide multi-motif scan
# to be useful; matched_sequence is dropped here to save memory, since this
# pipeline's outputs are interval/gene-level, not sequence-logo-level).
#
# Because p<=1e-5 and p<=1e-6 are strict subsets of what's loaded here (FIMO's
# own output threshold was 1e-4), this single load supports all three
# threshold tiers -- no need to re-read the files per threshold.

suppressPackageStartupMessages({
  library(data.table)
})
source("./PPRE_threshold_sensitivity/00_config.R")

all_hits_rds <- file.path(out_data, "all_hits_p1e-4.rds")

if (file.exists(all_hits_rds)) {
  message("Already exists, skipping load: ", all_hits_rds)
} else {
  hits_list <- vector("list", length(fimo_shards))
  names(hits_list) <- fimo_shards

  for (shard in fimo_shards) {
    f <- file.path(fimo_dir, shard, "fimo.tsv")
    stopifnot(file.exists(f))
    t0 <- Sys.time()
    # fread has no comment.char argument; FIMO appends 3 trailing "# ..." lines
    # after the data with only 1 field instead of 10, so read with fill=TRUE
    # and drop the resulting ragged/NA rows afterward instead.
    dt <- fread(f, sep = "\t", header = TRUE, fill = TRUE)
    dt <- dt[!is.na(start) & !startsWith(motif_id, "#")]
    dt <- dt[, .(motif_id, sequence_name = as.character(sequence_name),
                 start, stop, strand, score, pvalue = `p-value`)]
    dt[, shard := shard]
    hits_list[[shard]] <- dt
    cat(sprintf("%s: %s rows in %.1f sec\n", shard, format(nrow(dt), big.mark = ","),
                as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }

  all_hits <- rbindlist(hits_list, use.names = TRUE)
  rm(hits_list); gc()
  cat(sprintf("\nTotal hits (p<=1e-4) across all 30 motifs: %s\n", format(nrow(all_hits), big.mark = ",")))

  # sanity check: every motif_id present should be one of the 30 known motifs
  motif_family <- read.csv(motif_family_csv, stringsAsFactors = FALSE)
  unknown_motifs <- setdiff(unique(all_hits$motif_id), motif_family$motif_id)
  if (length(unknown_motifs) > 0)
    stop("FIMO output contains motif_id(s) not in the family table: ", paste(unknown_motifs, collapse = ", "))

  saveRDS(all_hits, all_hits_rds)
  message("Wrote: ", all_hits_rds)
}
