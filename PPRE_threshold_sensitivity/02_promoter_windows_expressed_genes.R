# Build -3000/+1000-of-TSS promoter windows (union across UCSC knownGene
# transcript isoforms per gene, same rule as PPARA_chipseq_overlap/01_get_promoter_windows.R
# and the ChIPseeker tssRegion=c(-3000,1000) calls elsewhere in this project),
# for the full cardiomyocyte-expressed gene sets rather than just the 5 genes
# PPARA_chipseq_overlap covered.
#
# Two gene sets, both already exported by Objective1_Expressed_genes_GTEx_TSS.Rmd
# (GTEx v9 snRNA-seq, Myocyte cells, Heart tissue) and re-used as-is here:
#   - "nonzero": any detectable expression (17,796 gene/Ensembl-ID rows)
#   - "top20":   top 20% by mean expression (3,536 gene/Ensembl-ID rows)

suppressPackageStartupMessages({
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(GenomicRanges)
  library(readr)
  library(dplyr)
})
source("./PPRE_threshold_sensitivity/00_config.R")

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
tx_by_gene <- transcriptsBy(txdb, by = "gene")   # names = Entrez IDs

build_promoter_bed <- function(gene_symbols, out_bed) {
  gene_symbols <- unique(gene_symbols)
  entrez_map <- AnnotationDbi::select(org.Hs.eg.db, keys = gene_symbols,
                                       keytype = "SYMBOL", columns = "ENTREZID")
  entrez_map <- entrez_map[!is.na(entrez_map$ENTREZID), ]
  entrez_map <- entrez_map[!duplicated(entrez_map$SYMBOL), ]  # keep first Entrez per symbol

  n_no_tx <- 0L
  promo_list <- vector("list", nrow(entrez_map))
  for (i in seq_len(nrow(entrez_map))) {
    eid <- entrez_map$ENTREZID[i]
    if (!(eid %in% names(tx_by_gene))) { n_no_tx <- n_no_tx + 1L; next }
    tx_gr <- tx_by_gene[[eid]]
    promo <- promoters(tx_gr, upstream = 3000, downstream = 1000)
    promo <- reduce(promo)
    mcols(promo)$gene <- entrez_map$SYMBOL[i]
    mcols(promo)$entrez_id <- eid
    promo_list[[i]] <- promo
  }
  promo_list <- promo_list[!vapply(promo_list, is.null, logical(1))]
  promoters_all <- sort(do.call(c, unname(promo_list)))

  # standard chromosomes only -- alt/patch scaffolds excluded (matches
  # PPARD_threshold_validation's documented convention: "alt/patch scaffolds
  # are excluded rather than reconciled")
  std_chroms <- paste0("chr", c(1:22, "X", "Y", "M"))
  promoters_all <- promoters_all[as.character(seqnames(promoters_all)) %in% std_chroms]

  df <- data.frame(
    chrom  = as.character(seqnames(promoters_all)),
    start  = start(promoters_all) - 1L,   # BED is 0-based
    end    = end(promoters_all),
    gene   = mcols(promoters_all)$gene,
    score  = 0,
    strand = as.character(strand(promoters_all))
  )
  write.table(df, out_bed, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

  cat(sprintf("%s: %d input gene symbols -> %d with an Entrez match -> %d with >=1 transcript -> %d promoter windows (%d genes with no txdb transcript, skipped)\n",
              basename(out_bed), length(gene_symbols), nrow(entrez_map),
              nrow(entrez_map) - n_no_tx, nrow(df), n_no_tx))
  invisible(df)
}

nonzero_genes <- read_csv(expr_nonzero_csv, show_col_types = FALSE)$HGNC_gene_name
top20_genes   <- read_csv(expr_top20_csv,   show_col_types = FALSE)$HGNC_gene_name

build_promoter_bed(nonzero_genes, promoter_windows_nonzero_bed)
build_promoter_bed(top20_genes,   promoter_windows_top20_bed)
