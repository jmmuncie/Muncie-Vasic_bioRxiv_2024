#Load libraries
library(ArchR); library(Seurat); library(dplyr)
library(GenomeInfoDb); library(ensembldb)
library(ggplot2); library(patchwork)
set.seed(1234)
library(BSgenome) 
library(hexbin)

#Set ArchR Threads to 4:
addArchRThreads(threads = 4, force = FALSE)

#Load ArchR project
proj_Mef2c_v13_E85_subset <- loadArchRProject("~/Mef2c_ArchR_working/Mef2c_v13_E85_subset")

#Set working directory
setwd("~/Mef2c_ArchR_working")


#Plot UMAPs
pal = c("pSHF" = "#FEE500", 
        "PhM" = "#C06CAB", 
        "CMs_IFT" = "#20A39E", 
        "LPM" = "#D8A767", 
        "NA" = "#999999", 
        "PostM" = "#208A42", 
        "C16" = "#89288F", 
        "CMs_OFT" = "#D1495B", 
        "aSHF" = "#90D5E4", 
        "CMs_V" = "#EDAE49")
UMAP_cell_types <- plotEmbedding(ArchRProj = proj_Mef2c_v13_E85_subset, 
                                 name = "Cluster_Labels_Subset_celltypegroup", 
                                 embedding = "UMAP_Harmony_Combined_sub", 
                                 size = 1,
                                 labelAsFactors=F, 
                                 labelMeans=T,
                                 labelSize = 3,
                                 pal = pal,
                                 baseSize = 10,
                                 colorTitle = ""
                                 ) + 
                                 theme(legend.text = element_text(size=10)) + 
                                 guides(colour = guide_legend(override.aes = list(size=3)))
cols = c("E8.5_KO1" = "#9C27B0",
         "E8.5_WT1" = "#F3E5F5", 
         "E8.5_WT2" = "#D1C4E9", 
         "E8.5_KO2" = "#4527A0"
         )
UMAP_sample_labels <- plotEmbedding(ArchRProj = proj_Mef2c_v13_E85_subset, 
                                    name = "SampleNames", 
                                    embedding = "UMAP_Harmony_Combined_sub", 
                                    size = 1, 
                                    labelAsFactors=F, 
                                    labelMeans=F, 
                                    pal = cols,
                                    colorTitle = ""
                                    ) + 
                                    theme(legend.text = element_text(size=10)) + 
                                    guides(colour = guide_legend(override.aes = list(size=3)))
plotPDF(UMAP_cell_types, name = "Fig3A_UMAP_E85_subset_cell_types", ArchRProj = proj_Mef2c_v13_E85_subset, addDOC = F)
plotPDF(UMAP_sample_labels, name = "Fig3A_UMAP_E85_subset_sample_labels", ArchRProj = proj_Mef2c_v13_E85_subset, addDOC = F)
