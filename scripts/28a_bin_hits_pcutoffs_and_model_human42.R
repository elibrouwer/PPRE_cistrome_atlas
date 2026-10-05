# Bin the 42-motif FIMO hits (500 bp) for (i) stricter p-value cutoffs and (ii) the distal ChIP-calibrated model.
# Per (chr, bin): max relative score (score / best score of its motif) over all 42 / PPAR-containing 14 / heterodimer 3 motifs at p <= 1e-4, 3e-5, 1e-5, 1e-6,
# and lp = sum over motifs of coef_m x (max relative score of motif m in the bin, p <= 1e-4), with coef_m from the adipocyte-PPARG elastic net (19; seed 11, same fit as 20).
suppressPackageStartupMessages({ library(data.table); library(glmnet) })
set.seed(11)
full_local <- "C:/Users/brouw/PPRE_cistrome_atlas/results/human42_pipeline_full_local"
fimo_root <- "C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Minor internship/PPAR_RXR_motifs_Hocomoco_Jaspar/FIMO_results_human42_unmasked"
inv <- fread("C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Minor internship/PPAR_RXR_motifs_Hocomoco_Jaspar/Motifs_update_2026-09-30/motif_inventory_gates.csv", encoding = "UTF-8")
cal <- readRDS(file.path(full_local, "chip_calibration_features.rds")); W <- cal$W; M <- cal$M; motifs <- colnames(M); ad <- which(W$group == "Adipocyte")
cv <- cv.glmnet(cbind(M, gc = W$gcf)[ad, ], W$y[ad], family = "binomial", alpha = 0.5, nfolds = 5)
cf <- as.matrix(coef(cv, s = "lambda.1se")); coefs <- setNames(cf[, 1], rownames(cf)); b0 <- coefs[["(Intercept)"]]; bgc <- coefs[["gc"]]; coefs <- coefs[motifs]
cat("non-zero motif coefficients:", sum(coefs != 0), " intercept:", round(b0, 3), " gc coef:", round(bgc, 3), " median window GC:", round(median(W$gcf), 3), "\n")
saveRDS(list(coefs = coefs, b0 = b0, bgc = bgc, gc_med = median(W$gcf)), file.path(full_local, "adipocyte_model_coefs.rds"))
cat_of <- setNames(inv$category, inv$motif_id); pp <- c("PPAR", "PPAR:RXR heterodimer")
cuts <- c(`1e-4` = 1e-4, `3e-5` = 3e-5, `1e-5` = 1e-5, `1e-6` = 1e-6)
acc <- list()
for (f in list.files(fimo_root, "fimo\\.tsv$", recursive = TRUE, full.names = TRUE)) {
  d <- fread(f, sep = "\t", select = c("motif_id", "sequence_name", "start", "stop", "score", "p-value"), colClasses = list(character = c("motif_id", "sequence_name")))
  d <- d[sequence_name %in% as.character(c(1:22, "X")) & !startsWith(motif_id, "#")]; setnames(d, "p-value", "p")
  d[, `:=`(start = as.integer(start), stop = as.integer(stop), score = as.numeric(score), p = as.numeric(p))]; d[, rel := score / max(score), by = motif_id]
  cg <- cat_of[d$motif_id]; d[, `:=`(chr = paste0("chr", sequence_name), bin = ((start + stop) %/% 2L) %/% 500L, is_pp = cg %in% pp, is_het = cg == "PPAR:RXR heterodimer")]
  base <- NULL
  for (nm in names(cuts)) { s <- d[p <= cuts[[nm]]]
    a <- s[, .(all = max(rel), pp14 = max(rel * is_pp), het3 = max(rel * is_het)), by = .(chr, bin)]; setnames(a, c("all", "pp14", "het3"), paste0(c("all", "pp14", "het3"), "_", nm))
    base <- if (is.null(base)) a else merge(base, a, by = c("chr", "bin"), all.x = TRUE) }
  m <- d[, .(mx = max(rel)), by = .(chr, bin, motif_id)][, .(lp = sum(coefs[motif_id] * mx)), by = .(chr, bin)]
  base <- merge(base, m, by = c("chr", "bin"), all.x = TRUE); acc[[f]] <- base; cat(basename(dirname(f)), nrow(d), "hits ->", nrow(base), "bins\n")
}
B <- rbindlist(acc, fill = TRUE); cols <- setdiff(names(B), c("chr", "bin", "lp"))
B <- B[, c(lapply(.SD, max, na.rm = TRUE), list(lp = sum(lp, na.rm = TRUE))), by = .(chr, bin), .SDcols = cols]
for (j in cols) set(B, which(!is.finite(B[[j]])), j, 0)
saveRDS(B, file.path(full_local, "fimo_human42_bins500_pcuts_lp.rds")); cat("bins:", nrow(B), "\n")
