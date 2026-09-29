library(anndata)
library(sceasy)
library(Seurat)
library(BPCells)
library(dplyr)
library(patchwork)
library(cowplot)
library(RColorBrewer)
library(ggplot2)
library(gridExtra)
library(grid)
library(cowplot)
library(ggplot2)
library(patchwork)
library(ggtext)
library(magrittr)
library(stringr)
library(ggVennDiagram)
library(readxl)
library(org.Hs.eg.db)
library(biomaRt)
library(DESeq2)
library(tidyr)

fibroblast_sample1 <- read.csv("data/Cultered_fibroblast/GSM8157769_HFF-50-CTR1.csv")
fibroblast_sample2 <- read.csv("data/Cultered_fibroblast/GSM8157770_HFF-50-CTR2.csv")


combined_fibroblast <- fibroblast_sample1 %>%
  left_join(fibroblast_sample2, keep = FALSE) 


rownames(combined_fibroblast) = combined_fibroblast$X
combined_fibroblast %<>%
  dplyr::select(-c(X, X.1))

library_sizes <- colSums(combined_fibroblast)
tp10k <- t(t(combined_fibroblast) / library_sizes) * 10000
tp10k_1 <- tp10k + 1

log10k_1 <- log10(tp10k_1)

fibroblast_seurat <- CreateSeuratObject(log10k_1)
fibroblast_seurat <- NormalizeData(fibroblast_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
fibroblast_seurat <- ScaleData(fibroblast_seurat)
calcNormFactors(dge)

fibro_gtex_seurat <- merge(GTEX_Seurat_subset, fibroblast_seurat, merge.data = TRUE)

fibro_gtex_seurat@meta.data %<>%
  mutate(Broad.cell.type = replace_na(Broad.cell.type, "Fibroblasts_cultered"))

fibro_gtex_seurat <- JoinLayers(fibro_gtex_seurat, layers= "counts")


bold_labels <- function(labels) {
  sapply(labels, function(label) {
    if (label %in% Combined_pam_pept_kin$hgnc_symbol) {
      paste0("**", label, "**") # Add bold markdown syntax
    } else {
      label
    }
  })
}

heatmap_GTEX_Fib_pam_peptides <- SCpubr::do_ExpressionHeatmap(
  fibro_gtex_seurat,
  slot = 'counts',
  features = Pam_peptide_data$hgnc_symbol,
  axis.text.x.angle = 90,
  group.by = c("Broad.cell.type"),
  groups.order = list( Broad.cell.type = sort(levels(factor(fibro_gtex_seurat@meta.data$Broad.cell.type)), decreasing = TRUE)),
  min.cutoff = 0,
  diverging.palette = "RdYlBu",
  use_viridis = FALSE)


heatmap_GTEX_Fib_pam_peptides <- heatmap_GTEX_Fib_pam_peptides + scale_fill_gradientn(
  colors = colorRampPalette(c("black", "lightblue", "darkblue"))(100),
  values = c(0, 0.001, 1),       # Adjust distribution
  limits = c(0, 3)
) + labs(title = "Expression of Pam-peptides in GTEx snRNAseq", y = "Cell type", x = "Pam gene petides") + scale_x_discrete(labels = bold_labels)  + theme(axis.title.x = element_text(size = 14, face = "bold", hjust = 0.5), axis.text.x = ggtext::element_markdown(size = 12), plot.title = element_text(size = 16, face = "bold", hjust = 0.5))


heatmap_GTEX_Fib_pam_kinase <- SCpubr::do_ExpressionHeatmap(
  fibro_gtex_seurat,
  slot = 'counts',
  features = Pam_kinase_data$`Kinase Name`,
  axis.text.x.angle = 90,
  group.by = c("Broad.cell.type"),
  legend.position = "none",
  groups.order = list( Broad.cell.type = sort(levels(factor(fibro_gtex_seurat@meta.data$Broad.cell.type)), decreasing = TRUE)),
  min.cutoff = 0,
  diverging.palette = "RdYlBu",
  use_viridis = FALSE)


heatmap_GTEX_Fib_pam_kinase <- heatmap_GTEX_Fib_pam_kinase + scale_fill_gradientn(
  colors = colorRampPalette(c("black", "lightblue", "darkblue"))(100),
  values = c(0, 0.001, 1),       # Adjust distribution
  limits = c(0, 3)
) + labs(title = "Expression of Pam-kinase in GTEx snRNAseq 112 list", y = "Cell type", x = "Pam gene kinase")  + scale_x_discrete(labels = bold_labels) +theme(axis.title.x = element_text(size = 14, face = "bold", hjust = 0.5), axis.text.x = ggtext::element_markdown(size = 12), plot.title = element_text(size = 16, face = "bold", hjust = 0.5))

pdf("data/Pam_chip/Pam_heatmap_GTEX_Fib_combined_newscale.pdf", height = 15, width = 30)
heatmap_arranged <- arrangeGrob(patchworkGrob(heatmap_GTEX_Fib_pam_peptides) , patchworkGrob(heatmap_GTEX_Fib_pam_kinase), nrow=2)
plot(heatmap_arranged)
dev.off()

pdf("data/Pam_chip/Summary_Pam_gene_table.pdf")

write.csv(Summary_Pam_gene_table, "data/Pam_chip/Summary_Pam_gene_table.csv", row.names = FALSE)

pdf("data/Pam_chip/Pam_heatmap_GTEX_Fib_Peptides.pdf",height =  20, width = 28)
heatmap_GTEX_Fib_pam_final_peptides
dev.off()

pdf("data/Pam_chip/Pam_heatmap_GTEX_Fib_kinases.pdf",height =  20, width = 28)
heatmap_GTEX_Fib_pam_kinase
dev.off()

pdf("data/Pam_chip/Pam_upsetplt.pdf", width = 8, height = 5)
Upsetplt
dev.off()
