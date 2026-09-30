"""Derive scripts/03b_fimo_chipseeker_PPRE_score_human42.Rmd from scripts/03_fimo_chipseeker_PPRE_score.Rmd.

Only the motif set is changed (42 human-source motifs instead of 30): FIMO input folder and batch names,
the two hard-coded motif-class tables and the data-type table (now built from
results/motif_selection/motif_inventory_human_only.csv), and the output folder. The analysis logic
(ChIPseeker window, sameStrand, scoring, GTEx overlap) is untouched.
Run from the repository root:  python scripts/00_motif_selection/make_03b_human42.py
"""
import re

src = open("scripts/03_fimo_chipseeker_PPRE_score.Rmd", encoding="utf-8").read()
out = src

def sub_once(pattern, repl, text, flags=0):
    new, n = re.subn(pattern, lambda m: repl, text, count=1, flags=flags)
    if n != 1:
        raise SystemExit(f"pattern not found: {pattern[:60]}")
    return new

out = sub_once(r'title: "FIMO Motif Analysis with ChIPseeker"', 'title: "FIMO Motif Analysis with ChIPseeker (42 human-source motifs)"', out)
out = sub_once(r"FIMO was run for all 30 motifs \(8 batch files.*?stored scores\.",
               "FIMO was run for the 42 human-source PPAR/RXR motifs (11 batch files, `RXR_PPAR_human42_files_1_4.meme` to `..._41_42.meme`) against the unmasked GRCh38 primary assembly, with an order-0 genome background, the default p-value threshold of 1e-4 and `--max-stored-scores 20000000` so FIMO does not discard hits at its default limit of 100,000 stored scores. The motif set is described in `results/motif_selection/`.",
               out, flags=re.S)
out = sub_once(r'no_max_fimo_dir <- "data/PPAR_RXR_motifs_Hocomoco_Jaspar/FIMO_results_no_max"',
               'no_max_fimo_dir <- "data/PPAR_RXR_motifs_Hocomoco_Jaspar/FIMO_results_human42_unmasked"', out)
out = sub_once(r'fimo_batches <- c\([^)]*\)',
               'fimo_batches <- c("1_4", "5_8", "9_12", "13_16", "17_20", "21_24", "25_28", "29_32", "33_36", "37_40", "41_42")', out)
out = sub_once(r'results_root <- "results"', 'results_root <- "results/human42_pipeline"', out)
# hits on either strand are assigned to the nearest TSS (sameStrand = FALSE): with sameStrand = TRUE the
# strongest promoter hits of e.g. CPT1A, CPT1B and HADHB (opposite strand) were dropped
out = sub_once(r'#ChIPseeker, specify the search window 3kb up and 1kb down, sameStrand = TRUE',
               '#ChIPseeker, specify the search window 3kb up and 1kb down, sameStrand = FALSE (either strand)', out)
out = sub_once(r'sameStrand = TRUE', 'sameStrand = FALSE', out)
out = out.replace("# Each of the 8 motif batches has its own fimo.tsv", "# Each of the 11 motif batches has its own fimo.tsv")
out = out.replace("# All 8 FIMO batch files (30 motifs) are annotated.", "# All 11 FIMO batch files (42 motifs) are annotated.")
out = out.replace("PPRE_FIMO_hits_30_motifs.bed", "PPRE_FIMO_hits_42_motifs.bed")

# motif metadata from the inventory, inserted right after the load-data chunk
meta = '''
```{r motif-metadata}
# Motif metadata (family, database, data type) comes from the motif selection, not from hand-typed tables.
motif_meta <- read.csv("results/motif_selection/motif_inventory_human_only.csv", check.names = FALSE,
                       fileEncoding = "UTF-8-BOM", stringsAsFactors = FALSE)
stopifnot(nrow(motif_meta) == 42)
motif_class <- setNames(
  ifelse(motif_meta$category == "PPAR:RXR heterodimer",
         paste0("Full-", toupper(motif_meta$name)), paste0("Half-", toupper(motif_meta$name))),
  motif_meta$motif_id)

simple_data_type <- function(x) dplyr::case_when(
  grepl("SMILE", x, ignore.case = TRUE) ~ "SMILE-seq",
  grepl("SELEX", x) ~ "SELEX",
  grepl("ChIP", x) ~ "ChIP-Seq",
  TRUE ~ x)
detailed_data_type <- function(db, data_type, source, native_id) {
  dplyr::case_when(
    db == "JASPAR" & data_type == "SMiLE-seq" ~ "SMILE-seq",
    db %in% c("JASPAR", "HOCOMOCO") ~ data_type,
    grepl("Kulakovskiy", source) ~ "ChIP-Seq",             # HOCOMOCO v11 models re-exported by CIS-BP
    grepl("Vorontsov", source) ~ "GHT-SELEX",
    grepl("Mathelier", source) ~ "SELEX",                    # JASPAR MA0066.2 re-exported by CIS-BP
    grepl("Isakova", source) ~ "SMILE-seq",
    grepl("Gerstein", source) ~ "ChIP-Seq",
    grepl("Methyl", native_id) ~ "HT-SELEX + Methyl-HT-SELEX",
    TRUE ~ "HT-SELEX")                                       # Jolma / Yin SELEX
}
```
'''
marker = 'dir.create(file.path(results_dir, "Scoring_FIMO"), showWarnings = FALSE, recursive = TRUE)\n```\n'
assert marker in out
out = out.replace(marker, marker + meta, 1)

# the two hard-coded motif_class definitions
pattern = r'motif_class <- c\(.*?"RXRA_M03623_3\.10" = "Half-RXRA"\s*\)'
out, n = re.subn(pattern, "# motif_class is defined in the motif-metadata chunk from the motif inventory", out, flags=re.S)
if n != 2:
    raise SystemExit(f"expected 2 motif_class blocks, replaced {n}")

# data-type table
pattern = r'measuring_method <- tibble::tribble\(.*?"RXRA_M03623_3\.10",\s*"HT-SELEX",\s*"SELEX"\s*\)'
new_mm = '''measuring_method <- motif_meta %>%
  dplyr::transmute(
    Motif_name = motif_id,
    Data_type_source = detailed_data_type(sub(" .*", "", database), data_type, source, native_id),
    Data_type_source_simple = simple_data_type(Data_type_source))'''
out, n = re.subn(pattern, lambda m: new_mm, out, flags=re.S)
if n != 1:
    raise SystemExit("measuring_method block not found")
out = out.replace('''      "SMILE-seq",
      "Transfac")), #Reorder the names''', '''      "SMILE-seq",
      "GHT-SELEX",
      "Transfac")), #Reorder the names''')
assert '"GHT-SELEX",\n      "Transfac"' in out

open("scripts/03b_fimo_chipseeker_PPRE_score_human42.Rmd", "w", encoding="utf-8", newline="\n").write(out)
print("written", len(out.splitlines()), "lines")
