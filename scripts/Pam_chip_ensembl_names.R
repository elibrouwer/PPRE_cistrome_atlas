#import pdf in R
library(pdftools)
library(tidyverse)
library(stringr)
library(biomaRt)
Pam_file_path = "data/Objective3_MultiOmics_overlaps_PPRE/Pam_chip/Pam_converted.csv"

Pam_table <- read.csv(Pam_file_path) %>%
  filter(No != "No") #filter the new columnames from the new page out

mart <- useEnsembl("ensembl", dataset="hsapiens_gene_ensembl", host = 'https://may2025.archive.ensembl.org')

a <- listFilters(mart)
Ensembl_list <- getBM(mart = mart,
           attribute = c('hgnc_symbol', 'ensembl_gene_id' ,'uniprot_gn_id'),
           filters = 'uniprot_gn_id',
           values = Pam_table$Uniprot.Accession)

Pam_table %>%
  left_join(Ensembl_list, by = c("Uniprot.Accession"= "uniprot_gn_id"), relationship = "many-to-many") -> finished_pam_table

write.csv(finished_pam_table, "data/Pam_chip/Pam_with_ensembl.csv", row.names = FALSE )            


# PAM_kinases -------------------------------------------------------------
Pam_kinase_path = "data/Pam_chip/112_Kinases.xlsx"
Pam_kinase_data <- read_excel(Pam_kinase_path, sheet = "Sheet2") 


mart <- useEnsembl("ensembl", dataset="hsapiens_gene_ensembl", host = 'https://jan2024.archive.ensembl.org')


Ensembl_list <- getBM(mart = mart,
                      attribute = c('hgnc_symbol', 'ensembl_gene_id' ,'uniprot_gn_id'),
                      filters = 'uniprot_gn_id',
                      values = Pam_kinase_data$`Kinase Uniprot ID`)

Pam_kinase_data %>%
  left_join(Ensembl_list, by = c("Kinase Uniprot ID"= "uniprot_gn_id"), relationship = "many-to-many") -> Pam_kinase_data_ensembl

write.csv(Pam_kinase_data_ensembl, "data/Pam_chip/Pam_112_kinases_with_ensembl.csv", row.names = FALSE )            


# PPRE ensembl ------------------------------------------------------------

PPRE_responsive_genes <- read.delim(PPRE_responsive_genes_path, sep = ",")

PPRE_responsive_genes <- colnames(PPRE_responsive_genes) %>%
  data.frame()

PPRE_responsive_genes_ensembl <- getBM(mart = mart,
                      attribute = c('hgnc_symbol', 'ensembl_gene_id' ),
                      filters = 'hgnc_symbol',
                      values = PPRE_responsive_genes$.)


write.csv(PPRE_responsive_genes_ensembl, "data/PPRE/PPRE_responsive_genes_ensembl.csv", row.names = FALSE )            

pdf("data/Objective3_MultiOmics_overlaps_PPRE/Pam_chip/Pam_heatmap_GTEX_Fib_combined_newscale1.pdf", height = 15, width = 30)
heatmap_arranged <- wrap_plots(patchworkGrob(heatmap_GTEX_Fib_pam_final_peptides) , patchworkGrob(heatmap_GTEX_Fib_pam_kinase), nrow=2)
plot(heatmap_arranged)
dev.off()

pdf("data/Objective3_MultiOmics_overlaps_PPRE/Pam_chip/Pam_heatmap_2026_GTEX_Fib_combined.pdf", height = 15, width = 30)
heatmap_arranged <- wrap_plots(patchworkGrob(heatmap_GTEX_Fib_pam_final_peptides_2026) , patchworkGrob(heatmap_GTEX_Fib_pam_kinase), nrow=2)
plot(heatmap_arranged)
dev.off()

pdf("data/Objective3_MultiOmics_overlaps_PPRE/Pam_chip/Pam_heatmap_2026_GTEX_Fib_kinases.pdf", height = 15, width = 30)
heatmap_arranged <- wrap_plots(patchworkGrob(heatmap_GTEX_Fib_pam_final_peptides_2026) , patchworkGrob(heatmap_GTEX_Fib_pam_final_peptides), nrow=2)
plot(heatmap_arranged)
dev.off()
