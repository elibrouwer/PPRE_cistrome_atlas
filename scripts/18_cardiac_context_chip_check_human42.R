# Is there any cardiac-context ChIP-seq that can serve as a gold standard for the PPRE score (42 motifs)?
# Available in the benchmark cache (ChIP-Atlas): PPARA in iPSC-derived cardiac cells (372 peaks), PPARG in HUVEC (endothelial, 9,603 peaks), no RXRA cardiac group.
# (a) gene-level: genes with >= 1 peak within promoter (-3k/+1k) or +-10 kb of a TSS = positives; matched-background AUROC of the score variants (as in 16b/17).
# (b) peak-level (iPSC cardiac PPARA): fraction of peak centres whose 500 bp bin holds a strong heterodimer-motif hit vs GC-matched control windows.
BENCH <- "C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Minor internship/R_code/benchmark"
source(file.path(BENCH, "config.R")); source(file.path(BENCH, "helpers_hits.R")); source(file.path(BENCH, "helpers_gold.R"))
suppressPackageStartupMessages(library(GenomicRanges))
set.seed(1)
out <- "C:/Users/brouw/PPRE_cistrome_atlas/results/human42_pipeline/distal_kernels"
gi <- load_genes(); tss <- readRDS(file.path(cache_dir, "tss_manuscript.rds")); K <- readRDS(file.path(out, "distal_kernel_scores.rds"))
rd <- function(f) { d <- fread(file.path(cache_dir, f), header = FALSE); setnames(d, c("chr", "start", "end", "srx", "group", "cell", "score")); d[chr %in% STD_CHR] }
pa <- rd("chipatlas_PPARA_min.tsv"); ipsc <- pa[grepl("cardiac", cell)]
pg <- rd("chipatlas_PPARG_min.tsv"); huvec <- pg[cell == "HUVEC"]
cat("iPSC cardiac PPARA peaks:", nrow(ipsc), "  HUVEC PPARG peaks:", nrow(huvec), "\n")

# (a) gene-level
u <- tss[in_universe == TRUE]
mkwin <- function(up, dn) { w <- GRanges(u$chr, IRanges(ifelse(u$strand == "+", u$tss - up, u$tss - dn), ifelse(u$strand == "+", u$tss + dn, u$tss + up))); w <- trim(suppressWarnings(w)); w$ENSEMBL <- u$ENSEMBL; w }
W <- list(promoter = mkwin(3000, 1000), `+-10kb` = mkwin(10000, 10000))
gr_of <- function(d) GRanges(d$chr, IRanges(d$start + 1L, d$end))
peaks <- list(`iPSC cardiac PPARA` = gr_of(ipsc), `HUVEC PPARG (score>=250)` = gr_of(huvec[score >= 250]), `HUVEC PPARG (all)` = gr_of(huvec))
lab <- list(); for (pn in names(peaks)) for (wn in names(W)) lab[[paste0(pn, " | ", wn)]] <- unique(W[[wn]]$ENSEMBL[overlapsAny(W[[wn]], peaks[[pn]])])
comp <- data.table(label = names(lab), n_genes = lengths(lab)); print(comp); fwrite(comp, file.path(out, "cardiac_chip_gene_labels.csv"))
lab <- lab[lengths(lab) >= 30]; lab <- lapply(lab, function(g) if (length(g) > 800) sample(g, 800) else g)
if (length(lab)) {
  S <- K[c("prom_base", "prom_het", "ac_tau20k_max", "rankavg_prom_ac", "wide_tau20k_max", "flat10k_sum", "flat50k_sum")]
  ev <- make_evaluator(gi, lab, ndraw = 300L)
  R <- rbindlist(lapply(names(S), function(nm) ev(S[[nm]], nm)$table)); fwrite(R, file.path(out, "cardiac_chip_gene_level_auroc.csv"))
  print(dcast(R, gold_set + n_gold ~ score, value.var = "auroc_matched")[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])
}

# (b) peak-level, iPSC cardiac PPARA: strong heterodimer hit in the peak-centre bin vs GC-matched 1 kb control windows
B <- readRDS("C:/Users/brouw/PPRE_cistrome_atlas/results/human42_pipeline_full_local/fimo_human42_bins500.rds")
gc <- fread("C:/Users/brouw/PPRE_cistrome_atlas/results/human42_pipeline/clustering/gc_bins_1kb.tsv"); gc <- gc[n == 0]; gc[, cls := floor(gc / 25)]
key <- function(chr, bin) paste(chr, bin); Bk <- B[, .(k = key(chr, bin), het_rel, max_rel)]; setkey(Bk, k)
look <- function(chr, pos, col) { v <- Bk[J(key(chr, pos %/% 500L))][[col]]; ifelse(is.na(v), 0, v) }
mid <- (ipsc$start + ipsc$end) %/% 2L
gcw <- merge(data.table(chrom = ipsc$chr, bin = mid %/% 1000L), gc, by = c("chrom", "bin"), all.x = TRUE)
ctrl <- rbindlist(lapply(seq_len(nrow(ipsc)), function(i) { if (is.na(gcw$cls[i])) return(NULL); p <- gc[cls == gcw$cls[i]][sample(.N, 10)]; data.table(chr = p$chrom, pos = p$bin * 1000L + sample(0:999, nrow(p), TRUE)) }))
res <- rbindlist(lapply(c(0.6, 0.7, 0.8, 0.9), function(th) { a <- mean(look(ipsc$chr, mid, "het_rel") >= th); b <- mean(look(ctrl$chr, ctrl$pos, "het_rel") >= th)
  data.table(het_rel_threshold = th, peaks = a, gc_matched_controls = b, fold = a / b) }))
print(res); fwrite(res, file.path(out, "ipsc_cardiac_ppara_peak_level_heterodimer_hit.csv"))
