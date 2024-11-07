#load Libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(magrittr)
library(XLConnect)
library(data.table)
library(harmony)
library(sctransform)

#Load Seurat Object
load("~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Objects/mef2c_v13_E85_scTrans_v3_harmony_subset.Robj")

#Set Working Directory
setwd("~/Desktop/Working_Directory")


#Analysis of NR2F2 Targets Upregulated in KO IFT-CMs

#Get list of all genes in the Seurat object -> 21185
all_genes <- rownames(mef2c_v13_E85_v3_harmony_subset)

#Get the list of Nr2f2 targets from the .tsv file downloaded from Harmonizome database
#Ensure only unique gene names 
data_nr <- read.delim("Nr2f2_Targets_Harmonizome_ENCODE_extracted.txt", header = F, sep = "\t")
Nr2f2_targets_raw <- as.character(data_nr$V1)
Nr2f2_targets <- unique(Nr2f2_targets_raw)

#Now we need to remove from the Nr2f2 targets list any genes that don't appear 
#in our Seurat object due to differences in genome annotations
match_Nr2f2_targets <- Nr2f2_targets[Nr2f2_targets %in% all_genes]
print(length(match_Nr2f2_targets))

#Total Nr2f2 targets -> 8163

#Get the list of up-regulated DEGs E8.5 IFT-CMs from the .txt file
#Or can repeat DEG analysis with FindMarkers()
data_DEGs <- read.delim("Mef2c_v13_E85_harmony_subset_markers_KOvWT_IFTCM_up_only.txt", header = TRUE, sep = "\t")
DEGs_raw <- data_DEGs$X
DEGs <- DEGs_raw[DEGs_raw != '']

#Find the number of Nr2f2 targets in the DEGs list -> 163
matches <- intersect(match_Nr2f2_targets, DEGs)
print(length(matches))

#Write out list of Nr2f2 targets upregulated
Nr2f2_targets_up_E85 <- data.frame(
  gene = matches
)
write.table(Nr2f2_targets_up_E85, "TableS7A_Nr2f2_Harmonizome_Targets_Upregulated_E85_IFT-CMs.tsv", sep="\t", row.names=FALSE, col.names=T, quote = FALSE)

#Get non-Nr2f2-targets, not up-regulated -> 12850
non_target_DEGs <- DEGs[!DEGs %in% matches]
non_Nr2f2_targets <- all_genes[!all_genes %in% match_Nr2f2_targets]
non_Nr2f2_targets_not_up <- non_Nr2f2_targets[!non_Nr2f2_targets %in% non_target_DEGs]
print(length(non_Nr2f2_targets_not_up))

#Fisher's Exact Test - for Fig6C
(tab <- matrix(c(163, 172, 8000, 12850), nrow = 2))
fisher.test(tab)


#NR2F2 Target Gene Expression

#Calculate module score using the Nr2f2 target genes that are up-regulated
mef2c_v13_E85_v3_harmony_subset <- AddModuleScore(mef2c_v13_E85_v3_harmony_subset, 
                                                  features = list(matches),
                                                  nbin = 10,
                                                  seed = 1,
                                                  name = "Nr2f2_Harmonizome_Targets_DEGs_Score"
)
Idents(mef2c_v13_E85_v3_harmony_subset) <- "cell_type_pool_x_genotype"
rp_Nr2f2 <- RidgePlot(object = mef2c_v13_E85_v3_harmony_subset,
                      features = "Nr2f2_Harmonizome_Targets_DEGs_Score1",
                      idents = c("CMs-IFT_WT", "CMs-IFT_KO", "CMs-V_WT", "CMs-V_KO", "CMs-OFT_WT", "CMs-OFT_KO")
)
ggsave('Fig6D_RidgePlot_E85_Nr2f2_Harmonizome_Targets_Up_Module_Score.pdf', plot = rp_Nr2f2, device = 'pdf', width = 7, height = 4, dpi = 300)


#For annotating CellOracle plots in Fig6I

#Get the list of CellOralce Nr2f2 Connections from the .txt files
CO_Nr2f2_WT <- read.delim("CellOracle_Nr2f2Connections_IFT-WT.txt", header = TRUE, sep = "\t")
CO_Nr2f2_WT_names <- as.character(CO_Nr2f2_WT$target)
CO_Nr2f2_KO <- read.delim("CellOracle_Nr2f2Connections_IFT-KO.txt", header = TRUE, sep = "\t")
CO_Nr2f2_KO_names <- as.character(CO_Nr2f2_KO$target)

#Look for Nr2f2 targets in connections
share_Nr2f2targets_CO_WT <- intersect(match_Nr2f2_targets, CO_Nr2f2_WT_names)
share_Nr2f2targets_CO_KO <- intersect(match_Nr2f2_targets, CO_Nr2f2_KO_names)

#Look for DEGs in connections
share_DEGs_CO_WT <- intersect(DEGs, CO_Nr2f2_WT_names)
share_DEGs_CO_KO <- intersect(DEGs, CO_Nr2f2_KO_names)


#Analysis of GATA4 Targets Upregulated in KO IFT-CMs

#Get list of all genes in the Seurat object -> 21185
all_genes <- rownames(mef2c_v13_E85_v3_harmony_subset)

#Get the list of Gata4 targets from the .tsv file downloaded from the TFLink database
#Ensure only unique gene names and
#Note: some target rows were annotated as "-" with no gene name, need to remove those
data <- read.delim("TFLink_targets_of_Gata4_Q08369.tsv", header = TRUE, sep = "\t")
Gata4_targets_raw <- as.character(data$Name.Target)
cleaned_vector <- Gata4_targets_raw[Gata4_targets_raw != "-"]
Gata4_targets <- unique(cleaned_vector)

#Now we need to remove from the Gata4 targets list any genes that don't appear 
#in our Seurat object due to differences in genome annotations
match_Gata4_targets <- Gata4_targets[Gata4_targets %in% all_genes]
print(length(match_Gata4_targets))

#Total Gata4 targets -> 12625

#Get the list of up-regulated DEGs in E8.5 IFT-CMs from the .txt file
#Or can repeat DEG analysis with FindMarkers()
data_DEGs <- read.delim("Mef2c_v13_E85_harmony_subset_markers_KOvWT_IFTCM_up_only.txt", header = TRUE, sep = "\t")
DEGs_raw <- data_DEGs$X
DEGs <- DEGs_raw[DEGs_raw != '']

#Find the number of Gata4 targets in the DEGs list -> 288
matches <- intersect(match_Gata4_targets, DEGs)
print(length(matches))

#Write out list of Gata4 targets upregulated
Gata4_targets_up_E85 <- data.frame(
  gene = matches
)
write.table(Gata4_targets_up_E85, "TableS7B_Gata4_TFLink_Targets_Upregulated_E85_IFT-CMs.tsv", sep="\t", row.names=FALSE, col.names=T, quote = FALSE)

#Get non-Gata4-targets, not up-regulated -> 8513
non_target_DEGs <- DEGs[!DEGs %in% matches]
non_Gata4_targets <- all_genes[!all_genes %in% match_Gata4_targets]
non_Gata4_targets_not_up <- non_Gata4_targets[!non_Gata4_targets %in% non_target_DEGs]
print(length(non_Gata4_targets_not_up))


#Fisher's Exact Test - for Fig6C
(tab <- matrix(c(288, 47, 12337, 8513), nrow = 2))
fisher.test(tab)


#Gata4 Target Expression

#Calculate module score using the Gata4 target genes that are up-regulated
mef2c_v13_E85_v3_harmony_subset <- AddModuleScore(mef2c_v13_E85_v3_harmony_subset, 
                                      features = list(matches),
                                      nbin = 10,
                                      seed = 1,
                                      name = "Gata4_Targets_DEGs_Score"
)
Idents(mef2c_v13_E85_v3_harmony_subset) <- "cell_type_pool_x_genotype"
rp_Gata4 <- RidgePlot(object = mef2c_v13_E85_v3_harmony_subset,
                      features = "Gata4_Targets_DEGs_Score1",
                      idents = c("CMs-IFT_WT", "CMs-IFT_KO", "CMs-V_WT", "CMs-V_KO", "CMs-OFT_WT", "CMs-OFT_KO")
)
ggsave('Fig6D_RidgePlot_E85_Gata4_TFLink_Targets_Up_Module_Score.pdf', plot = rp_Gata4, device = 'pdf', width = 7, height = 4, dpi = 300)


#For annotating CellOracle plots in Fig6J

#Get the list of CellOralce Gata4 Connections from the .txt files
CO_Gata4_WT <- read.delim("CellOracle_Gata4Connections_IFT-WT.txt", header = TRUE, sep = "\t")
CO_Gata4_WT_names <- as.character(CO_Gata4_WT$target)
CO_Gata4_KO <- read.delim("CellOracle_Gata4Connections_IFT-KO.txt", header = TRUE, sep = "\t")
CO_Gata4_KO_names <- as.character(CO_Gata4_KO$target)

#Look for Gata4 targets in connections
share_Gata4targets_CO_WT <- intersect(match_Gata4_targets, CO_Gata4_WT_names)
share_Gata4targets_CO_KO <- intersect(match_Gata4_targets, CO_Gata4_KO_names)

#Look for DEGs in connections
share_DEGs_CO_WT <- intersect(DEGs, CO_Gata4_WT_names)
share_DEGs_CO_KO <- intersect(DEGs, CO_Gata4_KO_names)


###----------Save Seurat Object---------------------------------------
save(mef2c_v13_E85_v3_harmony_subset, file = "~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Objects/mef2c_v13_E85_scTrans_v3_harmony_subset.Robj")

