#Load libraries
library(ArchR); library(Seurat); library(dplyr)
library(GenomeInfoDb); library(ensembldb)
library(ggplot2); library(patchwork)
set.seed(1234)
library(BSgenome) 
library(hexbin)

#Load ArchR project
proj_Mef2c_v13_E775_subset <- loadArchRProject("~/Mef2c_ArchR_working/Mef2c_v13_E775_subset")

#Set working directory
setwd("~/Mef2c_ArchR_working")


#Get DARs and Generate Volcano Plot
markerTest_CPs <- getMarkerFeatures(
  ArchRProj = proj_Mef2c_v13_E775_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "CellTypeByGenotype",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "CPs_x_KO",
  bgdGroups = "CPs_x_WT"
  )
VP <- plotMarkers(
  seMarker = markerTest_CPs, 
  name = "CPs_x_KO", 
  cutOff = "FDR <= 0.15 & abs(Log2FC) >= 0.5", 
  plotAs = "Volcano"
  ) + 
  theme(legend.text = element_text(size=12)) + 
  theme(axis.title = element_text(size = 12)) + 
  theme(axis.text = element_text(size = 12))    
plotPDF(VP, name = "FigS3D_VolcanoDARs_E775_CPs_KOvWT", width = 8, height = 5, ArchRProj = proj_Mef2c_v13_E775_subset, addDOC = FALSE)
