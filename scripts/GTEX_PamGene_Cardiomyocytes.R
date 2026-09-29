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
library(gplots)
library(EnhancedVolcano)


GTEX_file = "data/GTEX_single_nuclear_data/GTEx_8_tissues_snRNAseq_atlas_071421.public_obs.h5ad"
GTEX_metadata_path = "data/GTEX_single_nuclear_data/GTEx_8_annotation.csv"
GTEX_genes = "data/GTEX_single_nuclear_data/GTEx_8_genes.csv"

Pam_peptide_path = "data/Pam_chip/PamGene_kinases.csv"
Pam_peptide_data <- read.csv(Pam_peptide_path) %>%
  distinct(hgnc_symbol, .keep_all = TRUE) #removing duplicates 

Pam_kinase_path = "data/Pam_chip/112_Kinases.xlsx"
Pam_kinase_data <- read_excel(Pam_kinase_path, sheet = "Sheet2")

Pam_kinase_ensembl_path = "data/Pam_chip/Pam_112_kinases_with_ensembl.csv"
Pam_kinase_data_ensembl <- read.csv(Pam_kinase_ensembl_path)

GTEX_metadata <- read.csv(GTEX_metadata_path)
GTEX_genes <- read.csv(GTEX_genes)
 

parse.data <- open_matrix_anndata_hdf5(GTEX_file)

write_matrix_dir(mat = parse.data, dir = "data/GTEX_single_nuclear_data/GTEx_8_tissues_parsed")

GTEX_mat <- open_matrix_dir(dir = "data/GTEX_single_nuclear_data/GTEx_8_tissues_parsed")


# section to use if ensembl names -----------------------------------------
hgnc_rows <- rownames(GTEX_mat)

gene_data <- tibble(
  gene_name = hgnc_rows,
  row_index = seq_along(hgnc_rows) # Keep track of original order
)
aligned_genes <- GTEX_genes %>%
  left_join(gene_data, by = c( "X" = "gene_name"))
aligned_indices <- aligned_genes$row_index[!is.na(aligned_genes$row_index)]
rownames(GTEX_mat) <- aligned_genes$gene_ids


# SeuratObject ------------------------------------------------------------
GTEX_object <- CreateSeuratObject(GTEX_mat, meta.data = GTEX_metadata)
raw_counts <- GetAssayData(object = GTEX_object, layer = "counts")
matrix <- as.matrix(raw_counts)

# Number of genes and samples to sample
num_genes <- 1000
num_samples <- 6000

# Randomly sample rows and columns
subset_matrix <- matrix[
  sample(1:nrow(matrix), num_genes),
  sample(1:ncol(matrix), num_samples)
]
a <-summary(subset_matrix)
median(subset_matrix)
hist(as.vector(subset_matrix), breaks = 50, main = "Distribution of counts of GTEx", xlab = "Expression Values")

cardiac_cells <- GTEX_object@meta.data[grepl("Heart", GTEX_object@meta.data$Tissue.Site.Detail, ignore.case = TRUE), ]
skin_uv_cells <- GTEX_object@meta.data[grepl("skin", GTEX_object@meta.data$Tissue.Site.Detail, ignore.case = TRUE), ] %>%
  filter(Broad.cell.type == "Fibroblast")

subset_df <- rbind(cardiac_cells, skin_uv_cells) %>%
  distinct(X) #remove duplicates

GTEX_Seurat_subset <- subset(GTEX_object, cells = rownames(subset_df))
GTEX_Seurat_subset@assays$RNA@layers[["counts"]] <- as(GTEX_Seurat_subset@assays$RNA@layers[["counts"]], "dgCMatrix")

all.genes <- rownames(GTEX_Seurat_subset)
GTEX_Seurat_subset %<>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData( features = all.genes)


GTEX_Seurat_subset <- RunPCA(GTEX_Seurat_subset, features = VariableFeatures(object = GTEX_Seurat_subset), dims = 1:10)
GTEX_Seurat_subset %<>%
  FindNeighbors( dims = 1:10) %>%
  FindClusters( resolution = 0.5) %>%
  RunUMAP( dims = 1:10)

DimPlot(GTEX_Seurat_subset, reduction = "umap")

GTEX_Seurat_subset@meta.data %<>%
  mutate(Broad.cell.type = if_else(
    str_starts(Tissue.Site.Detail, "Skin"),
    str_c(Broad.cell.type, "_Skin"),
    Broad.cell.type))

GTEX_Seurat_subset@meta.data$Broad.cell.type

Idents(GTEX_Seurat_subset) <- "Broad.cell.type"
Cardiomyocyte_markers <- FindMarkers(GTEX_Seurat_subset, ident.1 = "Myocyte (cardiac, cytoplasmic)", ident.2 = "Myocyte (cardiac)")

Cardiomyocyte_markers$Significance <- ifelse(Cardiomyocyte_markers$p_val_adj < 0.05 & abs(Cardiomyocyte_markers$avg_log2FC) > 0.5, "Significant", "Not Significant")
Cardiomyocyte_markers$cluster <- Idents(GTEX_Seurat_subset)[match(rownames(Cardiomyocyte_markers), rownames(GTEX_Seurat_subset@assays$RNA@features))]

top_genes <- Cardiomyocyte_markers %>% 
  filter(p_val_adj < 0.05) %>%
  head(50)

Differentially_expressed_genes_Heatmap <- DoHeatmap(GTEX_Seurat_subset, features = rownames(top_genes), group.by = "Broad.cell.type")

Differentially_expressed_genes_Heatmap <- Differentially_expressed_genes_Heatmap  + scale_fill_gradientn( colors = colorpanel(100,"grey100","grey90","darkgreen"))

EnhancedVolcano(Cardiomyocyte_markers , 
                rownames(Cardiomyocyte_markers ),
                selectLab = rownames(top_genes),
                title = 'Classic versus cytoplasmic cardiomyocytes',
                x ="avg_log2FC", 
                y ="p_val_adj")
pdf("data/GTEX_single_nuclear_data/Differentially_expressed_genes_Cytoplasmic_Classic.pdf", width = 20, height = 10)
Differentially_expressed_genes_Heatmap
dev.off()


# venn diagram ------------------------------------------------------------
Combined_pam_pept_kin <- Pam_peptide_data %>%
  left_join(Pam_kinase_data, by = c("hgnc_symbol" = "Kinase Name")) %>%
  filter(!is.na(`Kinase Uniprot ID`)) #overlapping genes between peptides and kinases

venn <- ggVennDiagram(list(Pam_peptide_data$hgnc_symbol, Pam_kinase_data$`Kinase Name`), label = "both", label_geom = "label", category.names = c("Pam peptides", "Pam kinase"))
final_venn <- venn + scale_x_continuous(expand = expansion(mult = .2)) +   scale_fill_gradient(low = "white", high = "white") + labs(title = "Overlap between the peptides and kinases") + theme(legend.position = "none" , plot.title = element_text(size = 16, face = "bold", hjust = 0.5))

pdf("data/Pam_chip/Pam_venn_diagram.pdf")
final_venn
dev.off()
# violin plots ------------------------------------------------------------
pdf("data/Pam_chip/PamGene_kinases_gtex_violin.pdf")
for (i in 1:nrow(Pam_peptide_data)){
  print(i)
  try({
    plot <- VlnPlot(
      GTEX_Seurat_subset, 
      features = Pam_peptide_data$hgnc_symbol[i], 
      group.by = "Broad.cell.type",
      split.by = "Tissue.Site.Detail")
    print(plot)
  })
}
dev.off()




# Expression heatmaps -----------------------------------------------------
#creating a function to conditionally make the xlabs bold!
bold_labels <- function(labels) {
  sapply(labels, function(label) {
    if (label %in% Combined_pam_pept_kin$hgnc_symbol) {
      paste0("**", label, "**") # Add bold markdown syntax
    } else {
      label
    }
  })
}

heatmap_pam_peptides <- SCpubr::do_ExpressionHeatmap(
  GTEX_Seurat_subset,
  slot = 'counts',
  features = Pam_peptide_data$hgnc_symbol,
  axis.text.x.angle = 45,
  group.by = c("Broad.cell.type"),
  groups.order = list( Broad.cell.type = sort(levels(factor(GTEX_Seurat_subset@meta.data$Broad.cell.type)), decreasing = TRUE)),
  min.cutoff = 0,
  diverging.palette = "RdYlBu",
  use_viridis = FALSE)


heatmap_pam_peptides <- heatmap_pam_peptides + scale_fill_gradientn(
  colors = colorRampPalette(c("grey", "lightblue", "darkblue"))(100),
  values = c(0, 0.1, 1),       # Adjust distribution
  #limits = c(0, 3)
) + labs(title = "Expression of Pam-peptides in GTEx snRNAseq", y = "Cell type", x = "Pam gene petides") + scale_x_discrete(labels = bold_labels)  + theme(axis.title.x = element_text(size = 14, face = "bold", hjust = 0.5), axis.text.x = ggtext::element_markdown(size = 12), plot.title = element_text(size = 16, face = "bold", hjust = 0.5))


heatmap_pam_kinase <- SCpubr::do_ExpressionHeatmap(
  GTEX_Seurat_subset,
  slot = 'counts',
  features = Pam_kinase_data$`Kinase Name`,
  axis.text.x.angle = 45,
  group.by = c("Broad.cell.type"),
  legend.position = "none",
  groups.order = list( Broad.cell.type = sort(levels(factor(GTEX_Seurat_subset@meta.data$Broad.cell.type)), decreasing = TRUE)),
  min.cutoff = 0,
  diverging.palette = "RdYlBu",
  use_viridis = FALSE)


heatmap_pam_kinase <- heatmap_pam_kinase + scale_fill_gradientn(
  colors = colorRampPalette(c("grey", "lightblue", "darkblue"))(100),
  values = c(0, 0.1, 1),       # Adjust distribution
  limits = c(0, 3)
) + labs(title = "Expression of Pam-kinase in GTEx snRNAseq 112 list", y = "Cell type", x = "Pam gene kinase")  + scale_x_discrete(labels = bold_labels) +theme(axis.title.x = element_text(size = 14, face = "bold", hjust = 0.5), axis.text.x = ggtext::element_markdown(size = 12), plot.title = element_text(size = 16, face = "bold", hjust = 0.5))

pdf("data/Pam_chip/Pam_heatmap.pdf", height = 15, width = 30)
heatmap_arranged <- arrangeGrob(patchworkGrob(heatmap_pam_kinase) , patchworkGrob(heatmap_pam_peptides), nrow=2)
plot(heatmap_arranged)
dev.off()



# Overlap Pamgene (peptides + kinases) with PPRE regulated genes i --------
PPRE_responsive_genes_path = "data/PPRE/PPRE_responsive_genes.txt"
PPRE_bed_file_path = "data/PPRE/Bedfile_PPREs.xlsx"
PPRE_responsive_genes <- read.delim(PPRE_responsive_genes_path, sep = ",")
PPRE_responsive_genes <- colnames(PPRE_responsive_genes) %>%
  data.frame()

PPRE_bed_file_data <- read_excel(PPRE_bed_file_path, sheet = "Sheet1")
#cardiomyocyte specific
combined_peptide_PPRE_response <- Pam_peptide_data %>%
  left_join(PPRE_responsive_genes, by = c("hgnc_symbol" = "."), keep = TRUE) 

combined_kinase_PPRE_response <- Pam_kinase_data %>%
  right_join(PPRE_responsive_genes, by = c("Kinase Name" = "."), keep = TRUE)

venn_PPRE_peptides <- ggVennDiagram(list(Pam_peptide_data$hgnc_symbol, PPRE_responsive_genes$.), label = "both", label_geom = "label", category.names = c("Pam peptides", "PPRE genes"))
final_venn_PPRE_peptides <- venn_PPRE_peptides + scale_x_continuous(expand = expansion(mult = .2)) +   scale_fill_gradient(low = "white", high = "white") + labs(title = "Overlap between the peptides and PPRE influenced genes in cardiomyocytes") + theme(legend.position = "none" , plot.title = element_text(size = 16, face = "bold", hjust = 0.5))

venn_PPRE_kinases <- ggVennDiagram(list(Pam_kinase_data$`Kinase Name`, PPRE_responsive_genes$.), label = "both", label_geom = "label", category.names = c("Pam kinases", "PPRE genes"))
final_venn_PPRE_kinases <- venn_PPRE_kinases + scale_x_continuous(expand = expansion(mult = .2)) +   scale_fill_gradient(low = "white", high = "white") + labs(title = "Overlap between the kinases and PPRE influenced genes in cardiomyocytes") + theme(legend.position = "none" , plot.title = element_text(size = 16, face = "bold", hjust = 0.5))

pdf("data/Pam_chip/Pam_venn_diagram_PPRE_cardiomyocytes.pdf", width = 10)
venn_PPRE_arranged <- arrangeGrob(final_venn_PPRE_peptides , final_venn_PPRE_kinases, nrow=2)
plot(venn_PPRE_arranged)
dev.off()
#large bedfile
combined_peptide_PPRE_response <- Pam_peptide_data %>%
  left_join(PPRE_bed_file_data, by = c("hgnc_symbol" = "gene_name"), keep = TRUE) 

combined_kinase_PPRE_response <- Pam_kinase_data %>%
  left_join(PPRE_bed_file_data, by = c("Kinase Name" = "gene_name"))

venn_PPRE_bed_peptides <- ggVennDiagram(list(Pam_peptide_data$hgnc_symbol, PPRE_bed_file_data$gene_name), label = "both", label_geom = "label", category.names = c("Pam peptides", "PPRE genes"))
final_venn_PPRE_bed_peptides <- venn_PPRE_bed_peptides + scale_x_continuous(expand = expansion(mult = .2)) +   scale_fill_gradient(low = "white", high = "white") + labs(title = "Overlap between the peptides and PPRE influenced genes") + theme(legend.position = "none" , plot.title = element_text(size = 16, face = "bold", hjust = 0.5))

venn_PPRE_bed_kinases <- ggVennDiagram(list(Pam_kinase_data$`Kinase Name`, PPRE_bed_file_data$gene_name), label = "both", label_geom = "label", category.names = c("Pam kinases", "PPRE genes"))
final_venn_PPRE_bed_kinases <- venn_PPRE_bed_kinases + scale_x_continuous(expand = expansion(mult = .2)) +   scale_fill_gradient(low = "white", high = "white") + labs(title = "Overlap between the kinases and PPRE influenced genes") + theme(legend.position = "none" , plot.title = element_text(size = 16, face = "bold", hjust = 0.5))

pdf("data/Pam_chip/Pam_venn_diagram_PPRE_bed_cardiomyocytes.pdf", width = 10)
venn_PPRE_arranged <- arrangeGrob(final_venn_PPRE_bed_peptides , final_venn_PPRE_bed_kinases, nrow=2)
plot(venn_PPRE_arranged)
dev.off()

#MA11481
PPRE_MA11481_PPREs_path = "data/PPRE/MA11481_PPREs.xlsx"
PPRE_MA11481_PPREs_data <- read_excel(PPRE_MA11481_PPREs_path, sheet = "Sheet1")

combined_peptide_PPRE_MA11481_PPREs <- Pam_peptide_data %>%
  left_join(PPRE_MA11481_PPREs_data, by = c("ensembl_gene_id" = "ensembl_gene_id"), keep= TRUE) %>%
  filter(!is.na(ensembl_gene_id.y))

write.csv(combined_peptide_PPRE_MA11481_PPREs, "data/Pam_chip/Motif_venns/Pam_peptide_MA11481_overlap.csv")

combined_kinase_PPRE_MA11481_PPREs <- Pam_kinase_data_ensembl %>%
  left_join(PPRE_MA11481_PPREs_data, by = c("ensembl_gene_id" = "ensembl_gene_id"), keep= TRUE) %>%
  filter(!is.na(ensembl_gene_id.y))

write.csv(combined_kinase_PPRE_MA11481_PPREs, "data/Pam_chip/Motif_venns/Pam_kinase_MA11481_overlap.csv")

venn_PPRE_MA11481_PPREs_peptides <- ggVennDiagram(list(Pam_peptide_data$ensembl_gene_id, PPRE_MA11481_PPREs_data$ensembl_gene_id), label = "both", label_geom = "label", category.names = c("Pam peptides", "PPRE genes"))
final_PPRE_MA11481_PPREs_peptides <- venn_PPRE_MA11481_PPREs_peptides + scale_x_continuous(expand = expansion(mult = .2)) +   scale_fill_gradient(low = "white", high = "white") + labs(title = "Overlap between the peptides and PPRE MA11481 motif") + theme(legend.position = "none" , plot.title = element_text(size = 16, face = "bold", hjust = 0.5))

vennPPRE_MA11481_PPREs_kinases <- ggVennDiagram(list(Pam_kinase_data_ensembl$ensembl_gene_id, PPRE_MA11481_PPREs_data$ensembl_gene_id), label = "both", label_geom = "label", category.names = c("Pam kinases", "PPRE genes"))
final_PPRE_MA11481_PPREs_bed_kinases <- vennPPRE_MA11481_PPREs_kinases + scale_x_continuous(expand = expansion(mult = .2)) +   scale_fill_gradient(low = "white", high = "white") + labs(title = "Overlap between the kinases and PPRE MA11481 motif") + theme(legend.position = "none" , plot.title = element_text(size = 16, face = "bold", hjust = 0.5))

pdf("data/Pam_chip/Motif_venns/Pam_venn_diagram_PPRE_MA11481_cardiomyocytes.pdf", width = 10)
venn_PPRE_arranged <- arrangeGrob(final_PPRE_MA11481_PPREs_peptides , final_PPRE_MA11481_PPREs_bed_kinases, nrow=2)
plot(venn_PPRE_arranged)
dev.off()

#MA00651
PPRE_MA00651_PPREs_path = "data/PPRE/MA00651_PPREs.xlsx"
PPRE_MA00651_PPREs_data <- read_excel(PPRE_MA00651_PPREs_path, sheet = "Sheet1")

combined_peptide_PPRE_MA00651_PPREs <- Pam_peptide_data %>%
  left_join(PPRE_MA00651_PPREs_data, by = c("ensembl_gene_id" = "ensembl_gene_id"), keep= TRUE) %>%
  filter(!is.na(ensembl_gene_id.y))

write.csv(combined_peptide_PPRE_MA00651_PPREs, "data/Pam_chip/Motif_venns/Pam_peptide_MA00651_overlap.csv")

combined_kinase_PPRE_MA00651_PPREs <- Pam_kinase_data_ensembl %>%
  left_join(PPRE_MA00651_PPREs_data, by = c("ensembl_gene_id" = "ensembl_gene_id"), keep= TRUE) %>%
  filter(!is.na(ensembl_gene_id.y))

write.csv(combined_kinase_PPRE_MA00651_PPREs, "data/Pam_chip/Motif_venns/Pam_kinase_MA00651_overlap.csv")


venn_PPRE_MA00651_PPREs_peptides <- ggVennDiagram(list(Pam_peptide_data$ensembl_gene_id, PPRE_MA00651_PPREs_data$ensembl_gene_id), label = "both", label_geom = "label", category.names = c("Pam peptides", "PPRE genes"))
final_venn_PPRE_MA00651_PPREs_peptides <- venn_PPRE_MA00651_PPREs_peptides + scale_x_continuous(expand = expansion(mult = .2)) +   scale_fill_gradient(low = "white", high = "white") + labs(title = "Overlap between the peptides and PPRE MA00651 motif") + theme(legend.position = "none" , plot.title = element_text(size = 16, face = "bold", hjust = 0.5))

venn_PPRE_MA00651_PPREs_kinases <- ggVennDiagram(list(Pam_kinase_data_ensembl$ensembl_gene_id, PPRE_MA00651_PPREs_data$ensembl_gene_id), label = "both", label_geom = "label", category.names = c("Pam kinases", "PPRE genes"))
final_PPRE_MA00651_PPREs_kinases <- venn_PPRE_MA00651_PPREs_kinases + scale_x_continuous(expand = expansion(mult = .2)) +   scale_fill_gradient(low = "white", high = "white") + labs(title = "Overlap between the kinases and PPRE MA00651 motif") + theme(legend.position = "none" , plot.title = element_text(size = 16, face = "bold", hjust = 0.5))

pdf("data/Pam_chip/Motif_venns/Pam_venn_diagram_PPRE_MA00651_cardiomyocytes.pdf", width = 10)
venn_PPRE_arranged <- arrangeGrob(final_venn_PPRE_MA00651_PPREs_peptides , final_PPRE_MA00651_PPREs_kinases, nrow=2)
plot(venn_PPRE_arranged)
dev.off()

install.packages("UpSetR")
library(UpSetR)
upset_list <- list( Pam_kinase = Pam_kinase_data_ensembl$ensembl_gene_id, Pam_peptide = Pam_peptide_data$ensembl_gene_id , MA00651 = PPRE_MA00651_PPREs_data$ensembl_gene_id, MA11481 = PPRE_MA11481_PPREs_data$ensembl_gene_id, PPRE_genes = PPRE_responsive_genes_ensembl$ensembl_gene_id)
heado <- fromList(upset_list)
devtools::install_github("krassowski/complex-upset") 
library(ComplexUpset)

Upsetplt <- upset(heado, intersect = c("MA11481","MA00651", "PPRE_genes", "Pam_kinase", "Pam_peptide"), intersections = list(c("MA00651", "Pam_peptide"), 
                                                                                                                   c("MA11481", "Pam_peptide"),
                                                                                                                   c("MA00651", "Pam_kinase"),
                                                                                                                   c("MA11481", "Pam_kinase"),
                                                                                                                   c("MA00651", "MA11481", "Pam_kinase"),
                                                                                                                   c("MA11481", "MA00651", "Pam_peptide"),
                                                                                                                 c("PPRE_genes", "Pam_peptide"), 
                                                                                                                   c("PPRE_genes", "Pam_kinase"),
                                                                                                                   c("PPRE_genes", "MA00651"), 
                                                                                                                   c("PPRE_genes", "MA11481"),                                                                                              
                                                                                                                   c("Pam_kinase", "Pam_peptide")), 
      sort_intersections = FALSE, 
      sort_sets = FALSE, set_sizes=(upset_set_size() + ylab('Amount of genes')),
      stripes='white',
      height_ratio=0.4,
      width_ratio=0.2)

pdf("data/PPAR_RXR_motifs_Hocomoco_Jaspar/Objective1_results/Scoring_FIMO/PPRE_score_histograms.pdf", width = 11, height = 6)  # adjust size if you want
for (nm in names(plot_list)) {
  print(plot_list[[nm]])}
dev.off()