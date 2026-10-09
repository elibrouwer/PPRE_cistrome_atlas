# GSE262419 contrasts vs same-plate vehicle (DESIGN_V4.md Addendum 8). Sanity gate evaluated here, before any scoring.
suppressPackageStartupMessages(library(DESeq2))
proj <- "C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Minor internship/R_code"
g0 <- "C:/Users/brouw/AppData/Local/Temp/claude/C--Users-brouw-OneDrive---Universiteit-Utrecht-Master-Bioinformatics-Minor-internship-R-code/f7175609-6765-498e-8602-f582703f131b/scratchpad/gate0/"
hash <- read.csv(gzfile(paste0(g0, "GSE262419_hash.csv.gz")), stringsAsFactors = FALSE)
pad <- function(w) sub("^([A-Z]+)([0-9])$", "\\10\\2", w)
hash$col <- sprintf("Plate%d-%s", hash$Plate_ID, pad(hash$Well_ID))
contrasts <- list(list(14, "Pirinixic acid", 10), list(14, "Pirinixic acid", 1), list(10, "Troglitazone", 10), list(10, "Troglitazone", 1),
                  list(9, "Perfluorooctanoic acid", 10), list(8, "Perfluorononanoic acid", 10), list(12, "Perfluorodecanoic acid", 10), list(12, "Perfluorooctanesulfonic acid", 10))
known <- c("PDK4", "ANGPTL4", "CPT1A", "CPT1B")
cache <- list(); res_all <- list()
for (cc in contrasts) { tryCatch({
  p <- cc[[1]]; chem <- cc[[2]]; conc <- cc[[3]]
  key <- as.character(p)
  if (is.null(cache[[key]])) {
    m <- read.csv(gzfile(sprintf("%sGSE262419_Plate%02d.csv.gz", g0, p)), check.names = FALSE, stringsAsFactors = FALSE)
    sym <- sub("_[^_]*$", "", m$Genes); mat <- rowsum(as.matrix(m[, -1]), group = sym); cache[[key]] <- mat
  }
  mat <- cache[[key]]
  trt <- hash$col[hash$Plate_ID == p & hash$Chemical_name == chem & hash$Chemical_Concentration_uM == conc]
  veh <- hash$col[hash$Plate_ID == p & hash$Chemical_Index == "VEH"]
  trt <- intersect(trt, colnames(mat)); veh <- intersect(veh, colnames(mat))
  samples <- c(veh, trt)
  coldata <- data.frame(row.names = samples, cond = factor(c(rep("veh", length(veh)), rep("trt", length(trt))), levels = c("veh", "trt")))
  dds <- DESeqDataSetFromMatrix(countData = round(mat[, samples]), colData = coldata, design = ~ cond)
  dds <- estimateSizeFactors(dds); dds <- dds[rowMeans(counts(dds, normalized = TRUE)) >= 10, ]
  dds <- DESeq(dds, quiet = TRUE); r <- results(dds, name = "cond_trt_vs_veh")
  out <- data.frame(symbol = rownames(r), baseMean = r$baseMean, log2FC = r$log2FoldChange, stat = r$stat, padj = r$padj)
  tag <- paste0(gsub("[^A-Za-z0-9]", "", chem), "_", conc, "uM_plate", p)
  write.csv(out, file.path(proj, "results", "PPRE_score_v4", paste0("labels_gse262419_", tag, ".csv")), row.names = FALSE)
  ks <- out$stat[match(known, out$symbol)]
  gate <- sum(out$padj < 0.05, na.rm = TRUE) >= 50 && sum(ks > 2, na.rm = TRUE) >= 2
  cat(sprintf("%-45s trt=%d veh=%d genes=%d padj<0.05: %d | known stat: %s | GATE: %s\n", tag, length(trt), length(veh), nrow(out), sum(out$padj < 0.05, na.rm = TRUE),
              paste(known, round(ks, 1), collapse = " "), if (gate) "PASS" else "fail"))
}, error = function(e) cat("ERROR", cc[[2]], cc[[3]], "plate", cc[[1]], ":", conditionMessage(e), "
")) }
