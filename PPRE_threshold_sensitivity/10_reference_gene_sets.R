# Cache community-curated PPAR-related gene sets from MSigDB (msigdbr) for use
# as external validation targets -- NOT derived from this project's own FIMO
# output, so they don't circularly validate our own motif calls.
#
#   - PPAR_DR1_Q2 (MSigDB C3, TRANSFAC): genes with a computationally predicted
#     PPAR DR1 element near their TSS -- independently built years ago from a
#     different motif database/pipeline (TRANSFAC), the closest external
#     analogue to what this project's own DR1/FIMO analysis is trying to do.
#   - KEGG_PPAR_SIGNALING_PATHWAY (MSigDB C2): canonical KEGG pathway hsa03320
#     membership -- functional, not motif-based, so enrichment here is a
#     genuinely independent check.
#   - REACTOME_REGULATION_OF_LIPID_METABOLISM_BY_PPARALPHA (MSigDB C2):
#     curated Reactome pathway, PPARA-specific, directly relevant to the fatty-
#     acid-oxidation biology already used elsewhere in this project
#     (PPARA_chipseq_overlap's 5 genes: CPT1A, CPT1B, HADHA, HADHB, ACADVL).
#   - SANDERSON_PPARA_TARGETS (MSigDB C2): smaller experimentally-derived set.

suppressPackageStartupMessages({
  library(msigdbr)
  library(data.table)
})
source("./PPRE_threshold_sensitivity/00_config.R")

cache_rds <- file.path(out_data, "msigdbr_human_all.rds")
if (file.exists(cache_rds)) {
  m <- readRDS(cache_rds)
} else {
  m <- as.data.table(msigdbr(species = "Homo sapiens"))
  saveRDS(m, cache_rds)
}

target_sets <- c(
  "PPAR_DR1_Q2",
  "KEGG_PPAR_SIGNALING_PATHWAY",
  "REACTOME_REGULATION_OF_LIPID_METABOLISM_BY_PPARALPHA",
  "SANDERSON_PPARA_TARGETS"
)

gene_sets <- lapply(target_sets, function(gs) sort(unique(m[gs_name == gs, gene_symbol])))
names(gene_sets) <- target_sets

# This project's own project-confirmed PPARA targets (fatty-acid-oxidation
# genes with real ChIP-seq overlap checked in PPARA_chipseq_overlap/03_overlap_analysis.R)
gene_sets[["PROJECT_CONFIRMED_PPARA_FAO_GENES"]] <- c("CPT1A", "CPT1B", "HADHA", "HADHB", "ACADVL")

for (nm in names(gene_sets)) cat(sprintf("%-55s n=%d\n", nm, length(gene_sets[[nm]])))

saveRDS(gene_sets, file.path(out_data, "reference_gene_sets.rds"))
message("\nWrote: ", file.path(out_data, "reference_gene_sets.rds"))
