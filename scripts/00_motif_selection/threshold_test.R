suppressMessages({library(motifStack)})
motif_dir <- "C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Minor internship/PPAR_RXR_motifs_Hocomoco_Jaspar/Motifs_update_2026-09-30"
pfms <- importMatrix(file.path(motif_dir, "PPAR_RXR_human_only_combined.meme"))
names(pfms) <- sapply(pfms, function(m) m@name)
ma <- read.csv(file.path(motif_dir, "alignment_human42/manual_alignment_human42.csv"))
PPRE_mat <- matrix(c(1,1,1,0,0,1,0,0,0,0,1,1,1,0,0,0,0,1, 0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,1,0,
                     0,0,0,0,0,0,1,1,0,0,0,0,0,1,1,0,0,0, 0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,1,0,0), nrow = 4, byrow = TRUE)
rownames(PPRE_mat) <- c("A","C","G","T")
ppre <- new("pfm", mat = PPRE_mat, name = "PPRE motif")

revcomp_mat <- function(m) { m <- m[c("T","G","C","A"), ncol(m):1, drop = FALSE]; rownames(m) <- c("A","C","G","T"); m }
allowed <- rbind(c(1,0,1,0), c(0,0,1,0), c(0,0,1,0), c(0,0,0,1), c(0,1,0,0), c(1,0,0,0))
hs_scores <- function(m) { w <- ncol(m); if (w < 6) return(numeric(0))
  vapply(seq_len(w - 5), function(i) mean(vapply(1:6, function(k) sum(m[, i + k - 1] * allowed[k, ]), 0)), 0) }
hs_sites <- function(m, th = 0.6) { s <- hs_scores(m); cand <- which(s >= th); cand <- cand[order(-s[cand])]
  kept <- integer(0); for (p in cand) if (all(abs(p - kept) >= 6)) kept <- c(kept, p); sort(kept) }
col_ic <- function(m) { m <- sweep(m, 2, pmax(colSums(m), 1e-9), "/"); q <- pmax(m, 1e-12); 2 + colSums(m * log2(q)) }

evaluate <- function(al) {
  mats <- lapply(al, function(m) m@mat)
  nm <- sapply(al, function(m) m@name)
  keep <- !grepl("PPRE", nm)
  mats <- mats[keep]
  sites <- lapply(mats, hs_sites)
  has_plus <- mean(vapply(sites, length, 0) > 0)
  all_sites <- unlist(sites)
  tab <- sort(table(all_sites), decreasing = TRUE)
  top3 <- sum(head(tab, 3)) / sum(tab)
  mat3 <- simplify2array(mats); cons <- apply(mat3, c(1, 2), mean); cons <- sweep(cons, 2, colSums(cons), "/")
  c(n_rc = sum(grepl("\\(RC\\)$", nm[keep])), frac_with_plus_halfsite = round(has_plus, 3),
    frac_halfsites_in_top3_cols = round(top3, 3), consensus_total_IC = round(sum(col_ic(cons)), 2),
    max_col_IC = round(max(col_ic(cons)), 2))
}

res <- list()
pre <- pfms
for (i in seq_len(nrow(ma))) if (ma$Reverse_complement[i]) pre[[ma$Motif_name[i]]]@mat <- revcomp_mat(pre[[ma$Motif_name[i]]]@mat)
for (th in c(0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0)) {
  for (variant in c("default(revcomp=TRUE)", "pre-oriented(revcomp=FALSE)")) {
    lst <- if (grepl("pre", variant)) c(pre, list(ppre)) else c(pfms, list(ppre))
    names(lst) <- sapply(lst, function(m) m@name)
    rc <- if (grepl("pre", variant)) rep(FALSE, length(lst)) else rep(TRUE, length(lst))
    al <- tryCatch(DNAmotifAlignment(lst, threshold = th, minimalConsensus = 6, revcomp = rc), error = function(e) NULL)
    r <- if (is.null(al)) rep(NA, 5) else evaluate(al)
    res[[length(res) + 1]] <- data.frame(threshold = th, variant = variant, t(r))
    cat(sprintf("th=%.1f %-28s %s\n", th, variant, paste(names(r), r, collapse = " | ")), flush = TRUE)
  }
}
write.csv(do.call(rbind, res), "threshold_test_results.csv", row.names = FALSE)
