#Load libraries
library(ArchR); library(Seurat); library(dplyr)
library(GenomeInfoDb); library(ensembldb)
library(ggplot2); library(patchwork)
set.seed(1234)
library(BSgenome) 
library(hexbin)

#Set ArchR Threads to 4:
addArchRThreads(threads = 4, force = FALSE)

#Load ArchR projects
proj_Mef2c_v13_E775 <- loadArchRProject("~/Mef2c_ArchR_working/Mef2c_v13_E775_ATAC_and_GEX")
proj_Mef2c_v13_E85 <- loadArchRProject("~/Mef2c_ArchR_working/Mef2c_v13_E85_ATAC_and_GEX")
proj_Mef2c_v13_E9 <- loadArchRProject("~/Mef2c_ArchR_working/Mef2c_v13_E9_ATAC_and_GEX")

#Set working directory
setwd("~/Mef2c_ArchR_working")

#Plot UMAPs

#E7.75
cols_E775 = c("E7.75_WT1" = "lightskyblue",
              "E7.75_KO1" = "blue",
              "E7.75_KO2" = "darkblue",
              "E7.75_WT2" = "turquoise1"
              )
UMAP_E775_cell_types <- plotEmbedding(
  ArchRProj = proj_Mef2c_v13_E775, 
  name = "Cluster_Labels", 
  embedding = "UMAP_Harmony_Combined", 
  size = 0.5, 
  labelAsFactors=F, 
  labelMeans=T,
  labelSize = 3,
  baseSize = 10,
  colorTitle = ""
  ) + 
  theme(legend.text = element_text(size=10)) + guides(colour = guide_legend(override.aes = list(size=3)))
UMAP_E775_sample_labels <- plotEmbedding(
  ArchRProj = proj_Mef2c_v13_E775, 
  name = "SampleNames", 
  embedding = "UMAP_Harmony_Combined", 
  size = 0.5, 
  labelAsFactors=F, 
  labelMeans=F, 
  pal = cols_E775,
  colorTitle = ""
  ) + 
  theme(legend.text = element_text(size=10)) + guides(colour = guide_legend(override.aes = list(size=3)))
plotPDF(UMAP_E775_cell_types, name = "FigS3A_UMAP_E775_cell_types", ArchRProj = proj_Mef2c_v13_E775, addDOC = F)
plotPDF(UMAP_E775_sample_labels, name = "FigS3A_UMAP_E775_sample_labels", ArchRProj = proj_Mef2c_v13_E775, addDOC = F)

#E8.5
cols_E85 = c("E8.5_KO1" = "#9C27B0",
             "E8.5_WT1" = "#F3E5F5", 
             "E8.5_WT2" = "#D1C4E9", 
             "E8.5_KO2" = "#4527A0"
             )
UMAP_E85_cell_types <- plotEmbedding(
  ArchRProj = proj_Mef2c_v13_E85, 
  name = "Cluster_Labels", 
  embedding = "UMAP_Harmony_Combined", 
  size = 0.5, 
  labelAsFactors=F, 
  labelMeans=T,
  labelSize = 3,
  baseSize = 10,
  colorTitle = ""
  ) + 
  theme(legend.text = element_text(size=10)) + guides(colour = guide_legend(override.aes = list(size=3)))
UMAP_E85_sample_labels <- plotEmbedding(
  ArchRProj = proj_Mef2c_v13_E85, 
  name = "SampleNames", 
  embedding = "UMAP_Harmony_Combined", 
  size = 0.5, 
  labelAsFactors=F, 
  labelMeans=F, 
  pal = cols_E85,
  colorTitle = ""
  ) + 
  theme(legend.text = element_text(size=10)) + guides(colour = guide_legend(override.aes = list(size=3)))
plotPDF(UMAP_E85_cell_types, name = "FigS3B_UMAP_E85_cell_types", ArchRProj = proj_Mef2c_v13_E85, addDOC = F)
plotPDF(UMAP_E85_sample_labels, name = "FigS3B_UMAP_E85_sample_labels", ArchRProj = proj_Mef2c_v13_E85, addDOC = F)

#E9
cols_E9 = c("E9.0_WT1" = "pink",
            "E9.0_KO1" = "red", 
            "E9.0_KO2" = "red4", 
            "E9.0_WT2" = "plum1"
            )
UMAP_E9_cell_types <- plotEmbedding(
  ArchRProj = proj_Mef2c_v13_E9, 
  name = "Cluster_Labels", 
  embedding = "UMAP_Harmony_Combined", 
  size = 0.5, 
  labelAsFactors=F, 
  labelMeans=T,
  labelSize = 3,
  baseSize = 10,
  colorTitle = ""
  ) + 
  theme(legend.text = element_text(size=10)) + guides(colour = guide_legend(override.aes = list(size=3)))
UMAP_E9_sample_labels <- plotEmbedding(
  ArchRProj = proj_Mef2c_v13_E9, 
  name = "SampleNames", 
  embedding = "UMAP_Harmony_Combined", 
  size = 0.5, 
  labelAsFactors=F, 
  labelMeans=F, 
  pal = cols_E9,
  colorTitle = ""
  ) + 
  theme(legend.text = element_text(size=10)) + guides(colour = guide_legend(override.aes = list(size=3)))
plotPDF(UMAP_E9_cell_types, name = "FigS3C_UMAP_E9_cell_types", ArchRProj = proj_Mef2c_v13_E9, addDOC = F)
plotPDF(UMAP_E9_sample_labels, name = "FigS3C_UMAP_E9_sample_labels", ArchRProj = proj_Mef2c_v13_E9, addDOC = F)
