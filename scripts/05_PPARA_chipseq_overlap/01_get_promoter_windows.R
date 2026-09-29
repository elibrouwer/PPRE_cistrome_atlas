# Step 1: define promoter windows (-3000/+1000 of the TSS) for the known PPRE-responsive genes.
# The windows of all UCSC knownGene transcripts of a gene are merged, so alternative TSSs are
# included. This is the same window as in scripts 03 and 04.

suppressPackageStartupMessages({
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(GenomicRanges)
})

out_dir <- "results/PPARA_chipseq_overlap"

genes_symbols <- c("CPT1B", "HADHA", "HADHB", "CPT1A", "ACADVL")

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

entrez_map <- AnnotationDbi::select(org.Hs.eg.db, keys = genes_symbols,
                                     keytype = "SYMBOL", columns = "ENTREZID")
entrez_map <- entrez_map[!is.na(entrez_map$ENTREZID), ]

tx_by_gene <- transcriptsBy(txdb, by = "gene")

promo_list <- list()
for (i in seq_len(nrow(entrez_map))) {
  sym <- entrez_map$SYMBOL[i]
  eid <- entrez_map$ENTREZID[i]
  if (!(eid %in% names(tx_by_gene))) {
    cat("No transcripts found for", sym, "(entrez", eid, ") - skipping\n")
    next
  }
  tx_gr <- tx_by_gene[[eid]]
  promo <- promoters(tx_gr, upstream = 3000, downstream = 1000)
  promo <- reduce(promo)
  mcols(promo)$gene <- sym
  promo_list[[sym]] <- promo
}

promoters_all <- do.call(c, unname(promo_list))
promoters_all <- sort(promoters_all)

df <- data.frame(
  chrom = as.character(seqnames(promoters_all)),
  start = start(promoters_all) - 1L,   # BED is 0-based
  end   = end(promoters_all),
  gene  = mcols(promoters_all)$gene,
  score = 0,
  strand = as.character(strand(promoters_all))
)

write.table(df, file.path(out_dir, "gene_promoter_windows.bed"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

cat("\nPromoter windows (-3000/+1000 of TSS, union across isoforms):\n")
print(df)
cat("\nWrote", file.path(out_dir, "gene_promoter_windows.bed"), "\n")
