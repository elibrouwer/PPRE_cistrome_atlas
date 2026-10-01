# Does a wider window / distal weighting / H3K27ac-linked weighting improve the gene-level PPRE score (42 human-source motifs)?
# Input: 500 bp bins of the unmasked FIMO hits (16a). Gold sets and matched-background AUROC evaluator come from ./benchmark (PPARgene targets, mouse-heart PPARA-response sets).
# Score variants (value per bin = relative FIMO score, max over motifs):
#   prom_base          -3 kb/+1 kb window, exp(-|d|/3000), max     (binned version of the manuscript score, reference)
#   prom_het           same, heterodimer motifs only
#   wide_tau{10k,20k,50k}_{max,sum}   +-100 kb window, exp(-|d|/tau), max or sum over bins (all TSS of the gene)
#   flat{10k,50k}_sum  sum within +-W without decay
#   excl_tau20k_{max,sum}  each bin credited only to its nearest TSS (MAESTRO-like exclusivity)
#   ac_tau20k_{max,sum}    only bins inside heart H3K27ac differential regions (HCM or PLN, lifted to hg38), +-100 kb, tau 20 kb
#   prom_plus_ac_max       max(prom_base, ac_tau20k_max)
#   rankavg_prom_ac        mean rank of prom_base and ac_tau20k_max
BENCH <- "C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Minor internship/R_code/benchmark"
source(file.path(BENCH, "config.R")); source(file.path(BENCH, "helpers_hits.R")); source(file.path(BENCH, "helpers_gold.R"))
suppressPackageStartupMessages({ library(GenomicRanges); library(GenomeInfoDb); library(rtracklayer) })
out <- "C:/Users/brouw/PPRE_cistrome_atlas/results/human42_pipeline/distal_kernels"; dir.create(out, showWarnings = FALSE)
gi <- load_genes(); tss <- readRDS(file.path(cache_dir, "tss_manuscript.rds"))
B <- readRDS("C:/Users/brouw/PPRE_cistrome_atlas/results/human42_pipeline_full_local/fimo_human42_bins500.rds")
B[, pos := bin * 500L + 250L]

# ---- H3K27ac regions (same procedure as E2) ----
hcm <- suppressWarnings(readxl::read_xlsx(hcm_xlsx, sheet = "Supplementary Table 2A")); g_hcm <- GRanges(hcm$`Chromosome number`, IRanges(hcm$`Region start`, hcm$`Region end`))
pln <- suppressWarnings(readxl::read_xlsx(pln_xlsx, skip = 1, col_types = "text")); pln <- pln[!is.na(suppressWarnings(as.numeric(pln$`Region start†`))) & !is.na(pln$Chromosome), ]
g_pln <- GRanges(pln$Chromosome, IRanges(as.numeric(pln$`Region start†`), as.numeric(pln$`Region end†`)))
fix_chr <- function(g) { seqlevelsStyle(g) <- "UCSC"; g }
chain <- import.chain(chain_file)
lift <- function(g) { l <- unlist(liftOver(fix_chr(g), chain)); l[seqnames(l) %in% STD_CHR] }
ac <- reduce(c(lift(g_hcm), lift(g_pln)))
cat("H3K27ac regions after liftOver:", length(ac), " covering", round(sum(width(ac)) / 1e6, 1), "Mb\n")
B[, ac := countOverlaps(GRanges(chr, IRanges(pos, pos)), ac) > 0]
cat("bins in H3K27ac regions:", sum(B$ac), "of", nrow(B), "\n")

# ---- bin -> TSS joins ----
U <- unique(tss[in_universe == TRUE & ENSEMBL %in% gi$ENSEMBL, .(ENSEMBL, chr, strand, tss)])
U[, `:=`(g = match(ENSEMBL, gi$ENSEMBL), lo = tss - 100000L, hi = tss + 100000L, sgn = ifelse(strand == "+", 1L, -1L))]
cat("TSS rows:", nrow(U), "\n")
J <- B[U, on = .(chr, pos >= lo, pos <= hi), nomatch = 0L, allow.cartesian = TRUE, .(g = i.g, d = (x.pos - i.tss) * i.sgn, max_rel, het_rel, ac)]
cat("joined rows:", nrow(J), "\n")
nG <- nrow(gi)
S <- list()
P <- J[d >= -3000 & d <= 1000]; kP <- exp(-abs(P$d) / 3000)
S$prom_base <- { r <- rep(0, nG); a <- P[, .(v = max(max_rel * kP)), by = g]; r[a$g] <- a$v; r }
S$prom_het  <- { r <- rep(0, nG); a <- P[, .(v = max(het_rel * kP)), by = g]; r[a$g] <- a$v; r }
kern <- function(dt, tau) exp(-abs(dt$d) / tau)
for (tau in c(10000, 20000, 50000)) {
  k <- kern(J, tau); tg <- paste0("tau", tau / 1000, "k")
  S[[paste0("wide_", tg, "_max")]] <- { r <- rep(0, nG); a <- J[, .(v = max(max_rel * k)), by = g]; r[a$g] <- a$v; r }
  S[[paste0("wide_", tg, "_sum")]] <- { r <- rep(0, nG); a <- J[, .(v = sum(max_rel * k)), by = g]; r[a$g] <- a$v; r }
}
for (W in c(10000, 50000)) { s <- J[abs(d) <= W]; S[[paste0("flat", W / 1000, "k_sum")]] <- { r <- rep(0, nG); a <- s[, .(v = sum(max_rel)), by = g]; r[a$g] <- a$v; r } }
Jac <- J[ac == TRUE]; k <- kern(Jac, 20000)
S$ac_tau20k_max <- { r <- rep(0, nG); a <- Jac[, .(v = max(max_rel * k)), by = g]; r[a$g] <- a$v; r }
S$ac_tau20k_sum <- { r <- rep(0, nG); a <- Jac[, .(v = sum(max_rel * k)), by = g]; r[a$g] <- a$v; r }
S$prom_plus_ac_max <- pmax(S$prom_base, S$ac_tau20k_max)
S$rankavg_prom_ac <- (frank(S$prom_base) + frank(S$ac_tau20k_max)) / 2
# nearest-TSS exclusive version: every bin credited to its nearest TSS (universe genes + competitor genes)
Tall <- unique(tss[, .(ENSEMBL, chr, strand, tss)]); Tall[, `:=`(g = match(ENSEMBL, gi$ENSEMBL), sgn = ifelse(strand == "+", 1L, -1L))]
setorder(Tall, chr, tss)
nearest <- function(chr_, pos_) {
  T <- Tall[chr == chr_]; i <- findInterval(pos_, T$tss); i0 <- pmax(i, 1L); i1 <- pmin(i + 1L, nrow(T))
  pick <- ifelse(abs(pos_ - T$tss[i0]) <= abs(T$tss[i1] - pos_), i0, i1); list(g = T$g[pick], d = (pos_ - T$tss[pick]) * T$sgn[pick])
}
E <- rbindlist(lapply(unique(B$chr), function(ch) { b <- B[chr == ch]; n <- nearest(ch, b$pos); data.table(g = n$g, d = n$d, max_rel = b$max_rel) }))
E <- E[!is.na(g) & abs(d) <= 100000]; ke <- exp(-abs(E$d) / 20000)
S$excl_tau20k_max <- { r <- rep(0, nG); a <- E[, .(v = max(max_rel * ke)), by = g]; r[a$g] <- a$v; r }
S$excl_tau20k_sum <- { r <- rep(0, nG); a <- E[, .(v = sum(max_rel * ke)), by = g]; r[a$g] <- a$v; r }
S <- S[!vapply(S, is.null, TRUE)]
cat("scores:", paste(names(S), collapse = ", "), "\n")
saveRDS(S, file.path(out, "distal_kernel_scores.rds"))

# ---- evaluation on the benchmark gold sets (matched-background AUROC) ----
gold <- load_gold_all(gi); gold <- gold[c(PG_SETS, HEART_SETS)]
cat("gold set sizes:", paste(names(gold), lengths(gold), collapse = "; "), "\n")
ev <- make_evaluator(gi, gold, ndraw = 300L)
tabs <- lapply(names(S), function(nm) ev(S[[nm]], nm)$table)
R <- rbindlist(tabs)
fwrite(R, file.path(out, "distal_kernel_auroc_by_goldset.csv"))
W <- dcast(R, score ~ gold_set, value.var = "auroc_matched")
W[, `:=`(mean_PPARgene = rowMeans(.SD[, ..PG_SETS]), mean_heart = rowMeans(.SD[, ..HEART_SETS]))]
W[, overall := (mean_PPARgene + mean_heart) / 2]; W <- W[order(-overall)]
fwrite(W, file.path(out, "distal_kernel_auroc_summary.csv"))
print(W[, .(score, mean_PPARgene = round(mean_PPARgene, 3), mean_heart = round(mean_heart, 3), overall = round(overall, 3))])

# paired bootstrap vs the near-TSS reference for the top scores
top <- setdiff(head(W$score, 5), "prom_base")
bc <- boot_compare(S[c("prom_base", top)], gold, gi, ND = 100L, B = 500L)
d <- sapply(top, function(nm) {
  dpg <- rowMeans(bc$ba[, nm, PG_SETS] - bc$ba[, "prom_base", PG_SETS]); dh <- rowMeans(bc$ba[, nm, HEART_SETS] - bc$ba[, "prom_base", HEART_SETS])
  c(PPARgene = mean(dpg), PPARgene_lo = quantile(dpg, .025), PPARgene_hi = quantile(dpg, .975), heart = mean(dh), heart_lo = quantile(dh, .025), heart_hi = quantile(dh, .975))
})
cat("\npaired bootstrap difference vs prom_base (500 replicates):\n"); print(round(t(d), 3))
fwrite(as.data.table(t(d), keep.rownames = "score"), file.path(out, "distal_kernel_paired_bootstrap.csv"))
