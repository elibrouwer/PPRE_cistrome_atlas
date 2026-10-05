# PPARbeta/delta ChIP-seq in human myofibroblasts (WPMY-1; Adhikary et al. 2011 PLoS ONE, ArrayExpress E-MTAB-371, file PPAR_DMSO_peaks.tsv, 4,542 MACS peaks, hg19 -> hg38 by liftOver).
# (a) peak level: do FIMO hits of the motif subsets separate peak windows (500 bp at the peak centre) from 3 GC-matched control windows each?
# (b) gene level: genes with a PPARD peak in the promoter (-3 kb/+1 kb) and genes that are PPARD-bound AND induced by GW501516 (+ siPPARD-sensitive) vs matched genes; matched-background AUROC of the promoter scores.
BENCH <- "C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Minor internship/R_code/benchmark"
source(file.path(BENCH, "config.R")); source(file.path(BENCH, "helpers_hits.R")); source(file.path(BENCH, "helpers_gold.R"))
suppressPackageStartupMessages({ library(GenomicRanges); library(rtracklayer) })
set.seed(26)
out <- "C:/Users/brouw/PPRE_cistrome_atlas/results/human42_pipeline/pparD_myofibroblast"; dir.create(out, showWarnings = FALSE)
fimo_root <- "C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Minor internship/PPAR_RXR_motifs_Hocomoco_Jaspar/FIMO_results_human42_unmasked"
inv <- fread("C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Minor internship/PPAR_RXR_motifs_Hocomoco_Jaspar/Motifs_update_2026-09-30/motif_inventory_gates.csv", encoding = "UTF-8")
pk <- fread(file.path(BENCH, "external_data/PPARD_myofibroblast/PPAR_DMSO_peaks.tsv"), check.names = FALSE)
names(pk)[grepl("GW501516_allstars vs allstars_DMSO.*Largest_Status", names(pk)) & grepl("^Primary", names(pk))] <- "gw_status"
names(pk)[grepl("DMSO_siPPARD vs allstars_DMSO.*Largest_Status", names(pk)) & grepl("^Primary", names(pk))] <- "sipparD_status"
setnames(pk, c("Primary Gene ID", "Primary Gene Distance", "Primary Gene Name"), c("ENSG19", "pdist", "gname"))
cat("peaks:", nrow(pk), " GW501516 status table:\n"); print(table(pk$gw_status)); print(table(pk$sipparD_status))
gr19 <- GRanges(pk$chr, IRanges(pk$start + 1L, pk$end)); mcols(gr19) <- pk[, .(ENSG19, pdist, gname, gw_status, sipparD_status)]
chain <- import.chain(chain_file); lf <- liftOver(gr19, chain); keep <- lengths(lf) == 1; gr <- unlist(lf[keep]); gr <- gr[seqnames(gr) %in% paste0("chr", c(1:22, "X"))]
cat("lifted to hg38:", length(gr), "of", nrow(pk), "\n")

# ---- (a) peak-level ----
gc <- fread("C:/Users/brouw/PPRE_cistrome_atlas/results/human42_pipeline/clustering/gc_bins_1kb.tsv"); gc[, `:=`(gcf = gc / (length - n), cls = floor(100 * gc / (length - n) / 2.5))]; gcok <- gc[n == 0 & length == 1000]
mid <- (start(gr) + end(gr)) %/% 2L
pw <- data.table(chr = as.character(seqnames(gr)), pos = mid, y = 1L); pw[, bin := pos %/% 1000L]; pw <- merge(pw, gc[, .(chr = chrom, bin, gcf, cls)], by = c("chr", "bin"), all.x = TRUE); pw <- pw[!is.na(cls)]
cw <- rbindlist(lapply(split(pw, pw$cls), function(d) { pool <- gcok[cls == d$cls[1]]; if (!nrow(pool)) return(NULL); s <- pool[sample(.N, 3 * nrow(d), replace = TRUE)]
  data.table(chr = s$chrom, pos = s$bin * 1000L + sample(250:749, nrow(s), TRUE), y = 0L, gcf = s$gcf, cls = s$cls) }))
cw <- cw[!overlapsAny(GRanges(cw$chr, IRanges(cw$pos - 250L, cw$pos + 250L)), gr)]
W <- rbind(pw[, .(chr, pos, y, gcf)], cw[, .(chr, pos, y, gcf)]); W[, id := .I]; gw <- GRanges(W$chr, IRanges(W$pos - 250L, W$pos + 250L))
cat("windows: peaks", sum(W$y), " controls", sum(W$y == 0), "\n")
feat <- list()
for (f in list.files(fimo_root, "fimo\\.tsv$", recursive = TRUE, full.names = TRUE)) {
  d <- fread(f, sep = "\t", select = c("motif_id", "sequence_name", "start", "stop", "score"), colClasses = list(character = c("motif_id", "sequence_name")))
  d <- d[sequence_name %in% as.character(c(1:22, "X")) & !startsWith(motif_id, "#")]; d[, `:=`(start = as.integer(start), stop = as.integer(stop), score = as.numeric(score))]; d[, rel := score / max(score), by = motif_id]
  ov <- findOverlaps(GRanges(paste0("chr", d$sequence_name), IRanges(d$start, d$stop)), gw)
  feat[[f]] <- data.table(id = subjectHits(ov), motif = d$motif_id[queryHits(ov)], rel = d$rel[queryHits(ov)])[, .(rel = max(rel)), by = .(id, motif)] }
F <- rbindlist(feat); motifs <- sort(unique(inv$motif_id[inv$motif_id %in% F$motif])); M <- matrix(0, nrow(W), length(motifs), dimnames = list(NULL, motifs)); M[cbind(F$id, match(F$motif, motifs))] <- F$rel
auc <- function(x, y) { r <- rank(x); n1 <- sum(y == 1); n0 <- sum(y == 0); (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
cg <- setNames(inv$category, inv$motif_id)[motifs]
sub <- list(`all 42` = rep(TRUE, length(motifs)), `PPAR-containing 14` = cg %in% c("PPAR", "PPAR:RXR heterodimer"), `PPAR-only 11` = cg == "PPAR", `heterodimer 3` = cg == "PPAR:RXR heterodimer", `RXR-only 28` = cg == "RXR")
pl <- data.table(model = c("GC covariate only", names(sub)), AUROC = c(auc(W$gcf, W$y), sapply(sub, function(s) auc(apply(M[, s, drop = FALSE], 1, max), W$y))))
het <- cg == "PPAR:RXR heterodimer"; pl[, `:=`(pct_peaks_het_hit_rel0.7 = 100 * mean(apply(M[W$y == 1, het, drop = FALSE], 1, max) >= 0.7), pct_controls_het_hit_rel0.7 = 100 * mean(apply(M[W$y == 0, het, drop = FALSE], 1, max) >= 0.7))]
print(pl[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)]); fwrite(pl, file.path(out, "pparD_myofibroblast_peak_level.csv"))
# per-motif AUROC
pm <- data.table(motif_id = motifs, category = cg, AUROC = apply(M, 2, function(x) auc(x, W$y)))[order(-AUROC)]; fwrite(pm, file.path(out, "pparD_myofibroblast_per_motif_auroc.csv")); print(head(pm, 8))

# ---- (b) gene level ----
gi <- load_genes(); tss <- readRDS(file.path(cache_dir, "tss_manuscript.rds")); u <- tss[in_universe == TRUE]
w <- GRanges(u$chr, IRanges(ifelse(u$strand == "+", u$tss - 3000, u$tss - 1000), ifelse(u$strand == "+", u$tss + 1000, u$tss + 3000))); w <- trim(suppressWarnings(w)); w$ENSEMBL <- u$ENSEMBL
prom_genes <- unique(w$ENSEMBL[overlapsAny(w, gr)])
sy <- toupper(mcols(gr)$gname); induced_up <- unique(gi$ENSEMBL[toupper(gi$SYMBOL) %in% sy[mcols(gr)$gw_status == "u" & abs(mcols(gr)$pdist) <= 10000]])
dep <- unique(gi$ENSEMBL[toupper(gi$SYMBOL) %in% sy[mcols(gr)$gw_status == "u" & mcols(gr)$sipparD_status == "d" & abs(mcols(gr)$pdist) <= 10000]])
gold <- list(`PPARD promoter peak` = prom_genes, `PPARD peak <=10kb + GW501516-induced` = induced_up, `... + siPPARD-sensitive` = dep)
gold <- lapply(gold, intersect, gi$ENSEMBL); print(sapply(gold, length)); gold <- gold[lengths(gold) >= 15]
Sub <- readRDS("C:/Users/brouw/PPRE_cistrome_atlas/results/human42_pipeline/distal_kernels/motif_subset_scores.rds"); C <- readRDS("C:/Users/brouw/PPRE_cistrome_atlas/results/human42_pipeline/chip_calibration/chip_calibrated_gene_scores.rds")
J <- readRDS(file.path(res_dir, "J_candidate_gene_table.rds"))[match(gi$ENSEMBL, ENSEMBL)]
S <- list(all_42 = Sub$all_42, PPAR_containing_14 = Sub$PPAR_containing_14, PPAR_only_11 = Sub$PPAR_only_11, heterodimer_3 = Sub$heterodimer_3, RXR_only_28 = Sub$RXR_only_28, chip_model_top3_mean = C$chip_model_top3_mean,
          `promoter GC` = gi$gc, `cardiomyocyte expression` = ifelse(is.na(J$log10_expression), log10(1e-4), J$log10_expression))
G <- lapply(gold, function(g) if (length(g) > 800) sample(g, 800) else g)
ev <- NULL; for (s in 1:30) { set.seed(900 + s); ev <- tryCatch(make_evaluator(gi, G, ndraw = 300L), error = function(e) NULL); if (!is.null(ev)) break }
R <- rbindlist(lapply(names(S), function(nm) ev(S[[nm]], nm)$table)); fwrite(R, file.path(out, "pparD_myofibroblast_gene_level_auroc.csv"))
print(dcast(R, gold_set + n_gold ~ score, value.var = "auroc_matched")[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])
print(R[score %in% c("heterodimer_3", "PPAR_containing_14", "all_42"), .(gold_set, score, auroc_matched = round(auroc_matched, 3), lo = round(lo, 3), hi = round(hi, 3))])
