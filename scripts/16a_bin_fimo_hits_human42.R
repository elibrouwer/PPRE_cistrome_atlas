# Collapse the 45.7 M unmasked FIMO hits (42 human-source motifs, p <= 1e-4) to 500 bp bins so that distal kernels (up to +-100 kb of a TSS) can be scored.
# Per hit: rel = FIMO score / best score of that motif (relative score, comparable across motifs). Per (chr, bin): max_rel over all motifs, het_rel (max over PPAR:RXR heterodimer motifs), n hits.
# A batch contains whole motifs, so the motif maximum is computed inside each batch.
suppressMessages(library(data.table))
fimo_root <- "C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Minor internship/PPAR_RXR_motifs_Hocomoco_Jaspar/FIMO_results_human42_unmasked"
inv <- fread("C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Minor internship/PPAR_RXR_motifs_Hocomoco_Jaspar/Motifs_update_2026-09-30/motif_inventory_gates.csv", encoding = "UTF-8")
het <- inv[category == "PPAR:RXR heterodimer", motif_id]
BIN <- 500L; chrs <- c(1:22, "X")
out <- file.path("C:/Users/brouw/PPRE_cistrome_atlas/results/human42_pipeline_full_local"); dir.create(out, showWarnings = FALSE)
acc <- list()
for (f in list.files(fimo_root, "fimo\\.tsv$", recursive = TRUE, full.names = TRUE)) {
  d <- fread(f, sep = "\t", select = c("motif_id", "sequence_name", "start", "stop", "score"), colClasses = list(character = c("motif_id", "sequence_name")))
  d <- d[sequence_name %in% chrs & !startsWith(motif_id, "#")]
  d[, `:=`(start = as.integer(start), stop = as.integer(stop), score = as.numeric(score))]
  d[, rel := score / max(score), by = motif_id]
  d[, `:=`(bin = as.integer(((start + stop) %/% 2) %/% BIN), is_het = motif_id %in% het)]
  a <- d[, .(max_rel = max(rel), het_rel = if (any(is_het)) max(rel[is_het]) else 0, n = .N), by = .(chr = sequence_name, bin)]
  acc[[f]] <- a; cat(basename(dirname(f)), nrow(d), "hits ->", nrow(a), "bins\n")
}
B <- rbindlist(acc)[, .(max_rel = max(max_rel), het_rel = max(het_rel), n = sum(n)), by = .(chr, bin)]
B[, chr := paste0("chr", chr)]
saveRDS(B, file.path(out, "fimo_human42_bins500.rds"))
cat("bins:", nrow(B), " hits:", sum(B$n), "\n")
