# Cross-field alternatives to the genomics-standard combination methods
# (Fisher/Stouffer/Poisson), tested honestly against the same real ChIP label
# used in 15 -- not just proposed, actually cross-validated.
#
#   A. Learning-to-rank (information retrieval / search ranking field):
#      instead of a hand-picked combination formula, train a random forest
#      (ranger) and a regularized logistic regression (glmnet, LASSO) on ALL
#      engineered features at once to predict real ChIP support, evaluated
#      with proper k-fold cross-validation (AUC-ROC via pROC). This can find
#      nonlinear feature interactions no linear/statistical formula would.
#      Includes an explicit "promoter size alone" baseline model -- if the ML
#      model does no better than size alone, that means whatever it "learned"
#      is still just window size, not motif content.
#   B. Nearest-neighbor "comparables" (case-based reasoning / sports
#      scouting): z-score every feature, take the centroid of the 5
#      confirmed PPARA target genes in that z-scored feature space, and rank
#      every other gene by similarity (negative Euclidean distance) to that
#      centroid. No p-values or combination formula at all -- purely "which
#      genes look most like the known positives, overall."
#
# Both are validated the same two ways used throughout this analysis: (1)
# fgsea against PPAR_DR1_Q2 (community motif-based prediction, weak/circular
# form of validation) and (2) the real liver PPARA ChIP-seq label from 15
# (the meaningful bar, since it's actual experimental data, not another
# sequence-based prediction).

suppressPackageStartupMessages({
  library(data.table)
  library(ranger)
  library(glmnet)
  library(pROC)
  library(fgsea)
  library(GenomicRanges)
})
source("./PPRE_threshold_sensitivity/00_config.R")

scores <- fread(file.path(out_results, "gene_scoring_size_normalized.csv"))
dr1    <- fread(file.path(out_results, "dr1_hexamer_geometry_per_gene.csv"))[, .(gene, n_DR1_pairs)]
expr   <- fread(expr_nonzero_csv)[, .(gene = HGNC_gene_name, expression_value)]
expr   <- unique(expr, by = "gene")  # collapse any 1-to-many symbol->Ensembl duplicates

dat <- Reduce(function(a,b) merge(a,b,by="gene",all.x=TRUE), list(scores, dr1, expr))
dat[is.na(n_DR1_pairs), n_DR1_pairs := 0]
dat[, DR1_pairs_per_kb := n_DR1_pairs / (promoter_bp/1000)]
dat[is.na(expression_value), expression_value := 0]

## ---- real ChIP label (same as 15) ----
bed <- fread(promoter_windows_nonzero_bed, header = FALSE, col.names = c("chrom","start0","end","gene","score","strand"))
promo_gr <- GRanges(bed$chrom, IRanges(bed$start0+1L, bed$end), gene = bed$gene)
chip_bed <- file.path(repo_dir, "PPARA_chipseq_overlap", "chipatlas_liver_peaks_raw.bed")
raw <- fread(chip_bed, header = FALSE, col.names = c("chrom","start","end"))
std_chroms <- paste0("chr", c(1:22,"X","Y","M"))
raw <- raw[chrom %in% std_chroms]
chip_gr <- reduce(GRanges(raw$chrom, IRanges(raw$start+1L, raw$end)))
chip_supported_genes <- unique(promo_gr$gene[queryHits(findOverlaps(promo_gr, chip_gr, ignore.strand = TRUE))])
dat[, chip_supported := as.integer(gene %in% chip_supported_genes)]

project_genes <- c("CPT1A","CPT1B","HADHA","HADHB","ACADVL")
gene_sets <- readRDS(file.path(out_data, "reference_gene_sets.rds"))

feature_cols <- c("fisher_n_loci_combined","fisher_n_loci_heterodimer","fisher_n_loci_PPAR_half","fisher_n_loci_RXR_half",
                   "promoter_bp","n_windows",
                   "poisson_fold_combined","poisson_fold_heterodimer","poisson_fold_PPAR_half","poisson_fold_RXR_half",
                   "n_DR1_pairs","DR1_pairs_per_kb","expression_value")

## =========================== A. Learning-to-rank ===========================
set.seed(42)
K <- 5
dat[, fold := sample(rep(1:K, length.out = .N))]
dat[, `:=`(pred_rf = NA_real_, pred_size_only = NA_real_)]

for (k in 1:K) {
  train <- dat[fold != k]; test_idx <- dat[, which(fold == k)]
  rf <- ranger(x = train[, ..feature_cols], y = factor(train$chip_supported),
               probability = TRUE, num.trees = 500, seed = 42)
  pred <- predict(rf, data = dat[test_idx, ..feature_cols])$predictions[, "1"]
  dat[test_idx, pred_rf := pred]

  size_model <- glm(chip_supported ~ promoter_bp, data = train, family = binomial())
  dat[test_idx, pred_size_only := predict(size_model, newdata = dat[test_idx], type = "response")]
}

auc_rf   <- auc(dat$chip_supported, dat$pred_rf, quiet = TRUE)
auc_size <- auc(dat$chip_supported, dat$pred_size_only, quiet = TRUE)
cat(sprintf("Random forest (all features), 5-fold CV AUC vs. real ChIP label: %.3f\n", as.numeric(auc_rf)))
cat(sprintf("Promoter-size-ONLY baseline,  5-fold CV AUC vs. real ChIP label: %.3f  (0.5 = no better than chance)\n", as.numeric(auc_size)))

## variable importance from a single fit on all data (for interpretation only, not for the CV numbers above)
rf_full <- ranger(x = dat[, ..feature_cols], y = factor(dat$chip_supported),
                   probability = TRUE, num.trees = 1000, importance = "impurity", seed = 42)
imp <- sort(importance(rf_full), decreasing = TRUE)
cat("\nRandom forest variable importance (impurity), full-data fit:\n"); print(imp)

## LASSO logistic regression, for an interpretable linear-combination alternative
x_mat <- as.matrix(dat[, ..feature_cols]); x_mat <- scale(x_mat)
cv_lasso <- cv.glmnet(x_mat, dat$chip_supported, family = "binomial", alpha = 1)
lasso_coef <- coef(cv_lasso, s = "lambda.1se")
cat("\nLASSO logistic regression coefficients (lambda.1se, standardized features):\n")
print(lasso_coef)
dat[, pred_lasso := as.numeric(predict(cv_lasso, newx = x_mat, s = "lambda.1se", type = "response"))]
auc_lasso <- auc(dat$chip_supported, dat$pred_lasso, quiet = TRUE)
cat(sprintf("\nLASSO logistic regression AUC (in-sample, for reference): %.3f\n", as.numeric(auc_lasso)))

## =========================== B. Nearest-neighbor to known positives ===========================
z_mat <- scale(as.matrix(dat[, ..feature_cols]))
rownames(z_mat) <- dat$gene
centroid <- colMeans(z_mat[rownames(z_mat) %in% project_genes, , drop = FALSE])
dist_to_centroid <- sqrt(rowSums(sweep(z_mat, 2, centroid, "-")^2))
dat[, similarity_score := -dist_to_centroid[match(gene, rownames(z_mat))]]  # higher = more similar

## =========================== validate everything the same way ===========================
set.seed(1)
jitter <- function(x) x + rnorm(length(x), 0, 1e-9 * (max(x,na.rm=TRUE)-min(x,na.rm=TRUE)+1e-9))
candidate_scores <- list(
  ML_random_forest_CV     = dat$pred_rf,
  ML_lasso_logistic       = dat$pred_lasso,
  similarity_to_known_pos = dat$similarity_score
)
names(candidate_scores) <- names(candidate_scores)
val_rows <- list()
for (nm in names(candidate_scores)) {
  stat <- setNames(jitter(candidate_scores[[nm]]), dat$gene)
  r_dr1  <- as.data.table(fgsea(list(PPAR_DR1_Q2 = gene_sets[["PPAR_DR1_Q2"]]), stat, minSize=5, maxSize=500, eps=0, scoreType="pos"))
  r_chip <- as.data.table(fgsea(list(chip_supported = chip_supported_genes), stat, minSize=5, maxSize=length(chip_supported_genes)+100, eps=0, scoreType="pos"))
  val_rows[[nm]] <- data.table(method = nm,
                                NES_PPAR_DR1_Q2 = r_dr1$NES, padj_PPAR_DR1_Q2 = r_dr1$padj,
                                NES_real_ChIP = r_chip$NES, padj_real_ChIP = r_chip$padj)
}
val_tbl <- rbindlist(val_rows)
fwrite(val_tbl, file.path(out_results, "ml_and_similarity_validation.csv"))
cat("\n=== Validation: ML and similarity-based methods, both against DR1_Q2 (weak/circular) and real ChIP (meaningful) ===\n")
print(val_tbl)

## where do the known genes land under each?
pct_rank <- function(x) { r <- rank(-x, ties.method="average"); 100*(1-(r-1)/(length(x)-1)) }
known_check <- data.table(gene = project_genes,
                           RF_percentile = pct_rank(dat$pred_rf)[match(project_genes, dat$gene)],
                           similarity_percentile = pct_rank(dat$similarity_score)[match(project_genes, dat$gene)])
cat("\n=== Project-confirmed genes' percentile under the new methods ===\n"); print(known_check)

fwrite(dat[, .(gene, pred_rf, pred_lasso, similarity_score, chip_supported)],
       file.path(out_results, "ml_and_similarity_gene_scores.csv"))
message("\nDone.")
