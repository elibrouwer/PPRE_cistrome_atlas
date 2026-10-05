# Registry of every gold-standard gene set used so far (gene ENSEMBL ids of the benchmark universe), saved once so later scoring experiments (27b, 28, 29) use identical labels.
#   A  literature-validated human PPRE genes (tier 1, tier 1+2; promoter window)
#   B  PPARgene targets + mouse-heart PPARA-response sets
#   C  genes with a promoter-proximal peak (-3 kb/+1 kb): ChIP-Atlas PPARG per cell group, PPARA liver, RXRA liver; PPARD myofibroblast (E-MTAB-371)
#   D  genes with a peak within 10 kb / 50 kb of a TSS (distal-inclusive; positives that are NOT promoter peaks only) for PPARG adipocyte / liver / digestive tract and PPARD myofibroblast
#   E  PPARD myofibroblast: peak <= 10 kb and induced by GW501516
BENCH <- "C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Minor internship/R_code/benchmark"
source(file.path(BENCH, "config.R")); source(file.path(BENCH, "helpers_hits.R")); source(file.path(BENCH, "helpers_gold.R"))
suppressPackageStartupMessages({ library(GenomicRanges); library(rtracklayer) })
set.seed(27)
out <- "C:/Users/brouw/PPRE_cistrome_atlas/results/human42_pipeline/score_exploration"; dir.create(out, showWarnings = FALSE)
gi <- load_genes(); tss <- readRDS(file.path(cache_dir, "tss_manuscript.rds")); u <- tss[in_universe == TRUE]
mkw <- function(up, dn) { w <- GRanges(u$chr, IRanges(ifelse(u$strand == "+", u$tss - up, u$tss - dn), ifelse(u$strand == "+", u$tss + dn, u$tss + up))); w <- trim(suppressWarnings(w)); w$ENSEMBL <- u$ENSEMBL; w }
Wp <- mkw(3000, 1000); W10 <- mkw(10000, 10000); W50 <- mkw(50000, 50000)
genes_in <- function(W, gr) unique(W$ENSEMBL[overlapsAny(W, gr)])
reg <- list()
cap <- function(g, n = 800) { g <- intersect(g, gi$ENSEMBL); if (length(g) > n) sample(g, n) else g }
# B
gold <- load_gold_all(gi); for (n in c(PG_SETS, HEART_SETS)) reg[[paste0("B | ", n)]] <- gold[[n]]
# A
lit <- fread("C:/Users/brouw/Downloads/literature_PPRE_human_prioritised.csv")[tier %in% 1:2 & is_human == "yes"]
cand <- function(s) unique(toupper(trimws(unlist(strsplit(gsub("[()/,]", " ", s), "\\s+")))))
lit[, ENSEMBL := vapply(Gene_symbol, function(s) { m <- gi$ENSEMBL[toupper(gi$SYMBOL) %in% cand(s)]; if (length(m)) m[1] else NA_character_ }, character(1))]
lit[, rel := ifelse(is.na(end_num), start_num, (start_num + end_num) / 2)]
use <- lit[!is.na(ENSEMBL) & window_final == "yes" & ref_class == "TSS" & !is.na(rel) & rel >= -3000 & rel <= 1000]
reg[["A | literature tier 1"]] <- unique(use[tier == 1, ENSEMBL]); reg[["A | literature tier 1+2"]] <- unique(use$ENSEMBL)
# C / D from ChIP-Atlas
rd <- function(f) { d <- fread(file.path(cache_dir, f), header = FALSE); setnames(d, c("chr", "start", "end", "srx", "group", "cell", "score")); d[chr %in% STD_CHR & score >= 250][, group := gsub("%20", " ", group)] }
ca <- list(PPARG = rd("chipatlas_PPARG_min.tsv"), PPARA = rd("chipatlas_PPARA_min.tsv"), RXRA = rd("chipatlas_RXRA_min.tsv"))
sets <- list(c("PPARG", "Adipocyte"), c("PPARG", "Liver"), c("PPARG", "Digestive tract"), c("PPARG", "Blood"), c("PPARG", "Others"), c("PPARA", "Liver"), c("RXRA", "Liver"))
for (s in sets) { d <- ca[[s[1]]][group == s[2]]; gr <- GRanges(d$chr, IRanges(d$start + 1L, d$end)); p <- genes_in(Wp, gr); reg[[paste0("C | ", s[1], " ", s[2], " promoter peak")]] <- cap(p)
  if (s[1] == "PPARG" && s[2] %in% c("Adipocyte", "Liver", "Digestive tract")) { reg[[paste0("D | ", s[1], " ", s[2], " peak 3kb-10kb only")]] <- cap(setdiff(genes_in(W10, gr), p))
    reg[[paste0("D | ", s[1], " ", s[2], " peak 10kb-50kb only")]] <- cap(setdiff(genes_in(W50, gr), genes_in(W10, gr))) } }
# PPARD myofibroblast
pk <- fread(file.path(BENCH, "external_data/PPARD_myofibroblast/PPAR_DMSO_peaks.tsv"), check.names = FALSE)
names(pk)[grepl("GW501516_allstars vs allstars_DMSO.*Largest_Status", names(pk)) & grepl("^Primary", names(pk))] <- "gw_status"
setnames(pk, c("Primary Gene Distance", "Primary Gene Name"), c("pdist", "gname"))
g19 <- GRanges(pk$chr, IRanges(pk$start + 1L, pk$end)); mcols(g19) <- pk[, .(pdist, gname, gw_status)]
lf <- liftOver(g19, import.chain(chain_file)); gr <- unlist(lf[lengths(lf) == 1]); gr <- gr[seqnames(gr) %in% paste0("chr", c(1:22, "X"))]
pp <- genes_in(Wp, gr); reg[["C | PPARD myofibroblast promoter peak"]] <- cap(pp)
reg[["D | PPARD myofibroblast peak 3kb-10kb only"]] <- cap(setdiff(genes_in(W10, gr), pp)); reg[["D | PPARD myofibroblast peak 10kb-50kb only"]] <- cap(setdiff(genes_in(W50, gr), genes_in(W10, gr)))
sy <- toupper(mcols(gr)$gname); reg[["E | PPARD peak <=10kb + GW501516-induced"]] <- intersect(unique(gi$ENSEMBL[toupper(gi$SYMBOL) %in% sy[mcols(gr)$gw_status == "u" & abs(mcols(gr)$pdist) <= 10000]]), gi$ENSEMBL)
reg <- lapply(reg, function(g) intersect(g, gi$ENSEMBL)); print(data.table(set = names(reg), n = lengths(reg)), nrows = 60)
saveRDS(reg, file.path(out, "gold_registry.rds"))
