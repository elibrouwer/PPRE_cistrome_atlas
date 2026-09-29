# Reproduce the article's own gene-scoring method (TEPIC-style: Schmidt et al.
# 2017 -- FIMO motif strength weighted by exponential TSS-distance decay,
# summarized per gene as sum/mean/median/max) as faithfully as the Methods
# text supports, and test it -- not the alternative methods from 11-16 -- with
# the same rigor used throughout this analysis. Two changes from the article's
# current pipeline, both explicit and separately checkable:
#
#   1. Uses FIMO_results_no_max (untruncated, --max-stored-scores 20000000)
#      instead of the original default-settings FIMO run the article's
#      Methods describes -- filtered down to the SAME 22 JASPAR+HOCOMOCO
#      motifs the article used (excluding the 8 CIS-BP motifs added later in
#      this project), so the ONLY change from the article's stated pipeline
#      is fixing the truncation, not changing the motif set.
#   2. TSS-distance is computed per-transcript (not from this analysis's
#      earlier per-gene UNIONED promoter windows, which lose the exact TSS
#      position once multiple isoforms are merged and would make a distance-
#      decay term meaningless) -- each hit is assigned to its single nearest
#      transcript TSS genome-wide within -3000/+1000, mirroring ChIPseeker's
#      nearest-TSS annotation logic the article's Methods describes, rather
#      than this analysis's earlier isoform-unioned windows.
#
# The article does not state its exact decay constant beyond "exponential
# decay model" (TEPIC, Schmidt et al. 2017) -- decay_constant=1000bp is used
# here, chosen to be consistent with the article's own reported mean
# motif-to-TSS distance (739.7 +/- 61.1 bp). This is flagged explicitly as an
# assumption, not a literal reproduction of TEPIC's internal parameterization.

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(ranger)
  library(pROC)
  library(rtracklayer)
})
source("./PPRE_threshold_sensitivity/00_config.R")

DECAY_CONSTANT <- 1000  # bp; see note above

## ---- 1. restrict to the article's 22 motifs (exclude the 8 CIS-BP additions) ----
motif_family <- fread(motif_family_csv)
article_motifs <- motif_family[source != "CIS-BP", motif_id]
cat(sprintf("Using %d of 30 motifs (excluding %d CIS-BP additions) to match the article's motif set.\n",
            length(article_motifs), 30 - length(article_motifs)))

all_hits <- readRDS(file.path(out_data, "all_hits_p1e-4.rds"))
hits <- all_hits[motif_id %in% article_motifs]
cat(sprintf("Raw hits (p<=1e-4, article's 22 motifs, untruncated FIMO_results_no_max): %s\n", format(nrow(hits), big.mark=",")))

seqn <- ifelse(hits$sequence_name %in% c("MT","M"), "chrM", paste0("chr", hits$sequence_name))
std_chroms <- paste0("chr", c(1:22,"X","Y","M"))
keep <- seqn %in% std_chroms
hits <- hits[keep]; seqn <- seqn[keep]
hits_gr <- GRanges(seqn, IRanges(hits$start, hits$stop))
mcols(hits_gr)$pvalue <- hits$pvalue

## ---- 2. per-transcript TSS windows (NOT the isoform-unioned windows used elsewhere
## in this analysis -- unioning loses the exact TSS position a decay term needs) ----
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
tx_by_gene <- transcriptsBy(txdb, by = "gene")

expr_genes <- fread(expr_nonzero_csv)[, unique(HGNC_gene_name)]
entrez_map <- AnnotationDbi::select(org.Hs.eg.db, keys = expr_genes, keytype = "SYMBOL", columns = "ENTREZID")
entrez_map <- entrez_map[!is.na(entrez_map$ENTREZID) & entrez_map$ENTREZID %in% names(tx_by_gene), ]
entrez_map <- entrez_map[!duplicated(entrez_map$SYMBOL), ]

tx_all <- unlist(tx_by_gene[entrez_map$ENTREZID])
gene_lookup <- setNames(entrez_map$SYMBOL, entrez_map$ENTREZID)
mcols(tx_all)$gene <- gene_lookup[names(tx_all)]

tss <- ifelse(as.character(strand(tx_all)) == "-", end(tx_all), start(tx_all))
promo_tx <- GRanges(seqnames(tx_all),
                     IRanges(start = tss - 3000L, end = tss + 1000L),
                     strand = strand(tx_all), gene = mcols(tx_all)$gene, tss = tss)
cat(sprintf("Per-transcript promoter windows (NOT isoform-unioned): %d transcripts, %d genes\n",
            length(promo_tx), length(unique(promo_tx$gene))))

## ---- 3. assign each hit to its single NEAREST transcript TSS (mirrors ChIPseeker) ----
ov <- findOverlaps(hits_gr, promo_tx, ignore.strand = TRUE)
hit_mid <- (start(hits_gr) + end(hits_gr)) / 2
dist_to_tss <- abs(hit_mid[queryHits(ov)] - promo_tx$tss[subjectHits(ov)])

assign_tbl <- data.table(
  hit_idx = queryHits(ov), gene = promo_tx$gene[subjectHits(ov)],
  distance = dist_to_tss, pvalue = hits_gr$pvalue[queryHits(ov)]
)
setorder(assign_tbl, hit_idx, distance)
assign_tbl <- assign_tbl[!duplicated(hit_idx)]  # keep only the nearest gene per hit
cat(sprintf("Motif hits assigned to a nearest-TSS gene: %s (of %s promoter-overlapping hit-transcript pairs)\n",
            format(nrow(assign_tbl), big.mark=","), format(length(ov), big.mark=",")))

## ---- 4. TEPIC-style site score and gene-level summaries ----
assign_tbl[, site_score := -log10(pmax(pvalue, 1e-300)) * exp(-distance / DECAY_CONSTANT)]

gene_scores_tepic <- assign_tbl[, .(
  tepic_sum    = sum(site_score),
  tepic_mean   = mean(site_score),
  tepic_median = median(site_score),
  tepic_max    = max(site_score),
  n_sites      = .N
), by = gene]

# genes in the expressed universe with zero assigned sites (rare, but keep them at 0 for a fair ranking)
missing_genes <- setdiff(unique(promo_tx$gene), gene_scores_tepic$gene)
if (length(missing_genes) > 0) {
  gene_scores_tepic <- rbindlist(list(gene_scores_tepic, data.table(
    gene = missing_genes, tepic_sum = 0, tepic_mean = 0, tepic_median = 0, tepic_max = 0, n_sites = 0)))
}
fwrite(gene_scores_tepic, file.path(out_results, "article_tepic_style_gene_scores.csv"))
cat(sprintf("\nTEPIC-style gene scores computed for %d genes.\n", nrow(gene_scores_tepic)))

## ---- 5. the article's own face-validity check: where do its 4 focal genes rank? ----
focal_genes <- c("CPT1A", "HADHA", "HADHB", "ACADVL")
pct_rank <- function(x) { r <- rank(-x, ties.method = "average"); 100*(1-(r-1)/(length(x)-1)) }
focal_check <- data.table(gene = focal_genes)
for (col in c("tepic_sum","tepic_mean","tepic_median","tepic_max"))
  focal_check[[col]] <- pct_rank(gene_scores_tepic[[col]])[match(focal_genes, gene_scores_tepic$gene)]
cat("\n=== Article's own face-validity check: percentile rank (100=top) of its 4 focal genes ===\n")
print(focal_check)

message("\nStep 1 (scoring) done -- see 20_article_tepic_score_validation.R for the rigorous ChIP/cardiac validation.")
