# Pipeline for the 42-motif list, part 2: summary table, gene lists, Venn, H3K27ac overlap plot, new score list and ORA.
# Corresponds to the chunks summary-table, write-table-xlsx, gene-overlap, add-GTEX-nozero, venn-diagram, overlap-regions and run-ORA of the notebook.
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(stringr); library(readxl); library(writexl); library(ggplot2); library(ggvenn); library(org.Hs.eg.db); library(clusterProfiler)})
base <- "C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Minor internship"
res_root <- file.path(base, "PPAR_RXR_motifs_Hocomoco_Jaspar", "Human42_analysis"); res_dir <- file.path(res_root, "Objective1_results"); chip_dir <- file.path(res_root, "ChIPseeker_output")
v4 <- file.path(base, "R_code", "results", "PPRE_score_v4"); gtex <- file.path(base, "GTEX_single_nuclear_data")
log_ <- function(...) cat(sprintf(...), "\n")
# ---- 1. motif summary table
S <- read.csv(file.path(res_dir, "motif_summary_chipseeker_42.csv")); M <- read.csv(file.path(v4, "motif_summary_42.csv"))
tab <- S %>% left_join(M %>% select(motif, receptor, type, database), by = c("Motif_name" = "motif")) %>%
  transmute(`Motif name` = Motif_name, `Motif gene` = receptor, Family = type, Database = database, Peaks = N_amount, `Amount of genes` = N_amount_genes, `Average to TSS` = Mean_dist_tss, `Mean Width` = Mean_width)
ms <- function(x) sprintf("%.1f (%.1f)", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
sumrow <- tibble(`Motif name` = "Mean (SD)", `Motif gene` = "-", Family = "-", Database = "-", Peaks = ms(tab$Peaks), `Amount of genes` = ms(tab$`Amount of genes`), `Average to TSS` = ms(tab$`Average to TSS`), `Mean Width` = ms(tab$`Mean Width`))
tab_chr <- tab %>% mutate(Peaks = formatC(Peaks, format = "d", big.mark = ","), `Amount of genes` = formatC(`Amount of genes`, format = "d", big.mark = ","), `Average to TSS` = formatC(`Average to TSS`, format = "f", digits = 1), `Mean Width` = formatC(`Mean Width`, format = "f", digits = 1))
write_xlsx(list("PPRE table" = bind_rows(tab_chr, sumrow)), file.path(res_dir, "motif_results_42.xlsx"))
log_("motifs %d | total hits %s | mean hits per motif %s | mean width %s | mean distance to TSS %s", nrow(S), format(sum(S$N_amount), big.mark = ","), ms(S$N_amount), ms(S$Mean_width), ms(S$Mean_dist_tss))
# ---- 2. gene lists from the ChIPseeker promoter annotations
files <- list.files(chip_dir, pattern = "\\.txt$", full.names = TRUE)
chip_all <- do.call(rbind, lapply(files, function(f) read.delim(f, sep = " ", header = TRUE)))
gt <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(na.omit(chip_all$ENSEMBL)), keytype = "ENSEMBL", columns = c("ENSEMBL", "GENETYPE"))
chip_pc <- chip_all %>% left_join(gt, by = "ENSEMBL", relationship = "many-to-many") %>% filter(GENETYPE == "protein-coding")
genes_pc <- chip_pc %>% select(seqnames, motif_id, annotation, geneId, transcriptId, ENSEMBL, SYMBOL, pvalue, distanceToTSS) %>% distinct(transcriptId, motif_id, .keep_all = TRUE) %>%
  add_count(transcriptId, name = "n_motifs") %>% distinct(ENSEMBL, .keep_all = TRUE) %>% arrange(desc(n_motifs))
write.csv(genes_pc, file.path(res_dir, "PPRE_genes_protein_coding_42.csv"), row.names = FALSE)
cm <- read.csv(file.path(gtex, "GTEx_non_zeros_ensembl.csv")); cm_genes <- genes_pc %>% filter(ENSEMBL %in% cm$ensembl_gene_id)
write.csv(cm_genes, file.path(res_dir, "Motifs_GTEX_cm_expressed_42.csv"), row.names = FALSE)
log_("promoter rows %s | unique protein-coding genes with a promoter PPRE %d | of these expressed in cardiomyocytes (GTEx non-zero) %d | CM-expressed genes in GTEx list %d",
     format(nrow(chip_all), big.mark = ","), nrow(genes_pc), nrow(cm_genes), length(unique(cm$ensembl_gene_id)))
# ---- 3. Venn diagram
g8 <- read.csv(file.path(gtex, "GTEx_8_genes.csv"))
vl <- list("GTEX genes" = unique(g8$gene_name), "Expressed cardiomyocyte genes" = sort(unique(na.omit(cm$HGNC_gene_name))), "All Motif PPRE genes" = sort(unique(na.omit(genes_pc$SYMBOL))))
pdf(file.path(res_dir, "venn_diagram_42.pdf"), width = 8.27 / 2, height = 11.69 / 4); print(ggvenn(vl, show_percentage = FALSE, fill_color = c("white", "white", "white"), stroke_size = 1, set_name_size = 2.8, text_size = 2.8)); dev.off()
# ---- 4. H3K27ac overlap
O <- read.csv(file.path(res_dir, "H3K27ac_overlap_42.csv"))
nh <- O$regions_hcm_total[1]; np <- O$regions_pln_total[1]
O <- O %>% mutate(pct_hcm = regions_hcm_hit / nh * 100, pct_pln = regions_pln_hit / np * 100)
log_("H3K27ac regions after liftOver: HCM %d, PLN %d | percent of HCM regions overlapping at least one motif hit: %.2f +- %.2f SD (range %.1f-%.1f) | PLN: %.2f +- %.2f SD (range %.1f-%.1f)", nh, np, mean(O$pct_hcm), sd(O$pct_hcm), min(O$pct_hcm), max(O$pct_hcm), mean(O$pct_pln), sd(O$pct_pln), min(O$pct_pln), max(O$pct_pln))
long <- O %>% pivot_longer(c(total_hits, overlap_hcm, overlap_pln))
st <- ggplot(long, aes(motif_id, value, fill = name)) + geom_bar(stat = "identity", position = "stack") + theme_bw() + scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = c(total_hits = "#4A90E2", overlap_pln = "#50E3C2", overlap_hcm = "#D0021B"), labels = c(total_hits = "Total Hits", overlap_pln = "Overlap PLN", overlap_hcm = "Overlap HCM")) +
  theme(panel.border = element_blank(), panel.grid = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_blank()) + labs(x = "Motif ID", y = "Amount of FIMO peaks")
ggsave(file.path(res_dir, "H3K27ac_overlap_stackplot_42.pdf"), st, width = 12, height = 4); write.csv(O, file.path(res_dir, "H3K27ac_overlap_percent_42.csv"), row.names = FALSE)
# ---- 5. final score list and ORA
sc <- read.csv(file.path(v4, "PPRE_gene_scores_sitecount.csv"))
pc_expr <- read.csv(file.path(gtex, "GTEx_non_zeros_prt_coding.csv")); id_col <- grep("ensembl", names(pc_expr), ignore.case = TRUE, value = TRUE)[1]
sc_cm <- sc %>% filter(ENSEMBL %in% pc_expr[[id_col]]) %>% arrange(desc(PPRE_score_percentile), desc(n_sites))
write.csv(sc_cm, file.path(res_dir, "PPRE_score_sitecount_CM_expressed_protein_coding.csv"), row.names = FALSE)
write.table(data.frame(ensembl_genes = sc_cm$ENSEMBL, symbol = sc_cm$SYMBOL, score = sc_cm$PPRE_score_percentile), file.path(res_dir, "PPRE_sitecount_ranked.bed"), sep = "\t", quote = FALSE, row.names = FALSE)
log_("score list: %d genes expressed in cardiomyocytes (protein-coding)", nrow(sc_cm))
ora <- function(df, tag, use_universe = TRUE) {
  top <- head(df, 5000); ent <- bitr(top$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db); uni <- bitr(df$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  eg <- enrichGO(gene = ent$ENTREZID, universe = if (use_universe) uni$ENTREZID else NULL, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2, minGSSize = 10, readable = TRUE)
  r <- as.data.frame(eg); write.csv(r, file.path(res_dir, paste0("ORA_GOBP_", tag, ".csv")), row.names = FALSE); log_("ORA %s: %d significant GO BP terms; top: %s", tag, nrow(r), paste(head(r$Description, 5), collapse = "; "))
  if (nrow(r) > 0) { pdf(file.path(res_dir, paste0("ORA_dotplot_", tag, ".pdf")), width = 8, height = 7); print(dotplot(eg, showCategory = 20)); dev.off() }
}
ora(sc %>% arrange(desc(PPRE_score_percentile), desc(n_sites)), "all_genes_scored_universe_scored"); ora(sc_cm, "CM_expressed_universe_CM")
ora(sc %>% arrange(desc(PPRE_score_percentile), desc(n_sites)), "all_genes_scored_default_universe", FALSE); ora(sc_cm, "CM_expressed_default_universe", FALSE)
