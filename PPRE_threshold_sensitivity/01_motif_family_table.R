# Export a standalone motif-family classification table.
#
# This classification previously existed only inline in one chunk of
# Objective1_motif_alignment_CISBP_additions.Rmd (`manual_alignment` data
# frame, "Family" column: BOTH / PPAR / RXR). It is transcribed here verbatim
# and exported as a reusable CSV, then cross-checked against the MOTIF lines
# actually present in the 8 .meme shard files used by FIMO_results_no_max, so
# a silent mismatch (motif renamed/added/removed upstream) fails loudly
# instead of producing a table that quietly no longer matches the FIMO output.

suppressPackageStartupMessages(library(dplyr))
source("./PPRE_threshold_sensitivity/00_config.R")

motif_family <- tribble(
  ~motif_id,                  ~source,   ~family,      ~gene_symbol,     ~spacing_note,
  "MA1148.1",                 "JASPAR",  "heterodimer","PPARA::RXRA",    NA,
  "MA1148.2",                 "JASPAR",  "heterodimer","PPARA::RXRA",    NA,
  "MA0065.1",                 "JASPAR",  "heterodimer","PPARG::RXRA",    NA,

  "MA0066.1",                 "JASPAR",  "PPAR_half",  "PPARG",          NA,
  "MA0066.2",                 "JASPAR",  "PPAR_half",  "PPARG",          NA,
  "MA1550.1",                 "JASPAR",  "PPAR_half",  "PPARD",          NA,
  "MA1550.2",                 "JASPAR",  "PPAR_half",  "PPARD",          NA,
  "PPARA.H13CORE.0.P.B",      "HOCOMOCO","PPAR_half",  "PPARA",          NA,
  "PPARA.H13CORE.1.P.B",      "HOCOMOCO","PPAR_half",  "PPARA",          NA,
  "PPARD.H13CORE.0.PSM.A",    "HOCOMOCO","PPAR_half",  "PPARD",          NA,
  "PPARG.H13CORE.0.P.B",      "HOCOMOCO","PPAR_half",  "PPARG",          NA,
  "PPARG.H13CORE.1.P.B",      "HOCOMOCO","PPAR_half",  "PPARG",          NA,
  "PPARA_M12541_3.10",        "CIS-BP",  "PPAR_half",  "PPARA",          NA,
  "PPARD_M05794_3.10",        "CIS-BP",  "PPAR_half",  "PPARD",          "inverted/mirror spacing vs. canonical DR1 half-site (flagged in Objective1_motif_alignment_CISBP_additions.Rmd)",
  "PPARD_M05796_3.10",        "CIS-BP",  "PPAR_half",  "PPARD",          "inverted/mirror spacing vs. canonical DR1 half-site (flagged in Objective1_motif_alignment_CISBP_additions.Rmd)",
  "PPARG_M12496_3.10",        "CIS-BP",  "PPAR_half",  "PPARG",          NA,
  "PPARG_M12497_3.10",        "CIS-BP",  "PPAR_half",  "PPARG",          NA,
  "PPARG_M12498_3.10",        "CIS-BP",  "PPAR_half",  "PPARG",          NA,

  "MA0855.1",                 "JASPAR",  "RXR_half",   "RXRB",           NA,
  "MA0856.1",                 "JASPAR",  "RXR_half",   "RXRG",           NA,
  "MA1555.1",                 "JASPAR",  "RXR_half",   "RXRB",           NA,
  "MA1556.1",                 "JASPAR",  "RXR_half",   "RXRG",           NA,
  "RXRA.H13CORE.0.PS.A",      "HOCOMOCO","RXR_half",   "RXRA",           NA,
  "RXRA.H13CORE.1.S.C",       "HOCOMOCO","RXR_half",   "RXRA",           NA,
  "RXRA.H13CORE.2.SM.B",      "HOCOMOCO","RXR_half",   "RXRA",           NA,
  "RXRA.H13CORE.3.P.B",       "HOCOMOCO","RXR_half",   "RXRA",           NA,
  "RXRB.H13CORE.0.P.C",       "HOCOMOCO","RXR_half",   "RXRB",           NA,
  "RXRB.H13CORE.1.SM.B",      "HOCOMOCO","RXR_half",   "RXRB",           NA,
  "RXRA_M03621_3.10",         "CIS-BP",  "RXR_half",   "RXRA",           "everted/DR0-like spacing vs. canonical DR1 half-site (flagged in Objective1_motif_alignment_CISBP_additions.Rmd)",
  "RXRA_M03623_3.10",         "CIS-BP",  "RXR_half",   "RXRA",           "everted/DR0-like spacing vs. canonical DR1 half-site (flagged in Objective1_motif_alignment_CISBP_additions.Rmd)"
)

stopifnot(nrow(motif_family) == 30, !anyDuplicated(motif_family$motif_id))

## ---- cross-check against the actual .meme shard files used by FIMO_results_no_max ----
meme_files <- file.path(base_dir, paste0("RXR_PPAR_combined_files_", fimo_shards, ".meme"))
stopifnot(all(file.exists(meme_files)))

motifs_in_meme <- unlist(lapply(meme_files, function(f) {
  lines <- readLines(f)
  m <- grep("^MOTIF", lines, value = TRUE)
  vapply(strsplit(trimws(m), "\\s+"), `[`, character(1), 2)
}))

missing_from_table <- setdiff(motifs_in_meme, motif_family$motif_id)
extra_in_table      <- setdiff(motif_family$motif_id, motifs_in_meme)

if (length(missing_from_table) > 0)
  stop("Motif(s) in the .meme shards have no family assignment: ",
       paste(missing_from_table, collapse = ", "))
if (length(extra_in_table) > 0)
  stop("Motif(s) in the family table are not present in the .meme shards used by FIMO_results_no_max: ",
       paste(extra_in_table, collapse = ", "))

cat("Cross-check OK: all", length(motifs_in_meme), "motifs in the 8 .meme shards match the family table exactly.\n")
cat("Family counts: ", paste(names(table(motif_family$family)), table(motif_family$family), sep = "=", collapse = ", "), "\n")

write.csv(motif_family, motif_family_csv, row.names = FALSE)
message("Wrote: ", motif_family_csv)
