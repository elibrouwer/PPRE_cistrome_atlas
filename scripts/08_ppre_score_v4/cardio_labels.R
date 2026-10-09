# Single-agonist vs Control DESeq2 labels for GSE160987 (DESIGN_V4.md Addendum 7). Written before any scoring.
suppressPackageStartupMessages(library(DESeq2))
proj <- "C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Minor internship/R_code"
raw <- read.csv(gzfile(file.path(proj, "PPARA_validation_GSE244905", "data", "GSE160987_20201015_Dubois_counts.csv.gz")), check.names = FALSE, stringsAsFactors = FALSE)
raw <- raw[!is.na(raw$gene_symbol) & complete.cases(raw[, -1]), ]
raw <- rowsum(as.matrix(raw[, -1]), group = raw$gene_symbol)
ctrl <- grep("^Control rep [0-9]+$", colnames(raw), value = TRUE)
for (ag in c("WY14643", "Rosiglitazone", "GW0742")) {
  trt <- grep(paste0("^", ag, " rep [0-9]+$"), colnames(raw), value = TRUE)
  cat(ag, "n =", length(trt), "vs control n =", length(ctrl), "\n")
  samples <- c(ctrl, trt)
  coldata <- data.frame(row.names = samples, cond = factor(c(rep("ctrl", length(ctrl)), rep("trt", length(trt))), levels = c("ctrl", "trt")))
  dds <- DESeqDataSetFromMatrix(countData = round(raw[, samples]), colData = coldata, design = ~ cond)
  dds <- estimateSizeFactors(dds); dds <- dds[rowMeans(counts(dds, normalized = TRUE)) >= 10, ]
  dds <- DESeq(dds, quiet = TRUE); res <- results(dds, name = "cond_trt_vs_ctrl")
  out <- data.frame(symbol = rownames(res), baseMean = res$baseMean, log2FC = res$log2FoldChange, stat = res$stat, padj = res$padj)
  write.csv(out, file.path(proj, "results", "PPRE_score_v4", paste0("labels_cardio_", ag, ".csv")), row.names = FALSE)
  known <- c("PDK4", "ANGPTL4", "CPT1A", "CPT1B")
  cat("tested", nrow(out), "| padj<0.05:", sum(out$padj < 0.05, na.rm = TRUE), "| known targets stat:", paste(known, round(out$stat[match(known, out$symbol)], 1), collapse = ", "), "\n")
}
