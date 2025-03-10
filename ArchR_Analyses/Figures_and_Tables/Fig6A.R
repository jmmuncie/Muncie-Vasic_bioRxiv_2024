#Load libraries
library(ArchR); library(Seurat); library(dplyr)
library(GenomeInfoDb); library(ensembldb)
library(ggplot2); library(patchwork)
set.seed(1234)
library(BSgenome) 
library(hexbin)

#Load ArchR project
proj_Mef2c_v13_E85_subset <- loadArchRProject("~/Mef2c_ArchR_working/Mef2c_v13_E85_subset")

#Set working directory
setwd("~/Mef2c_ArchR_working")


#Get DARs and Plot Motif Enrichments

#IFT-CMs
#Look at motifs enriched in peaks with +Log2FC in markerTest_CMs 
#i.e. the peaks that gain accessibility in MEF2C KO
markerTest_CMs_IFT <- getMarkerFeatures(
  ArchRProj = proj_Mef2c_v13_E85_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "CellTypeByGenotype",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "CMs_IFT_x_KO",
  bgdGroups = "CMs_IFT_x_WT"
)
motifsKO_IFT <- peakAnnoEnrichment(
  seMarker = markerTest_CMs_IFT,
  ArchRProj = proj_Mef2c_v13_E85_subset,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.15 & Log2FC >= 0.5"
)
#To prepare for plotting, create simple df containg
#the motif names, the corrected p-values, and the significance rank
dfmKO_IFT <- data.frame(TF = rownames(motifsKO_IFT), mlog10Padj = assay(motifsKO_IFT)[,1])
dfmKO_IFT <- dfmKO_IFT[order(dfmKO_IFT$mlog10Padj, decreasing = TRUE),]
dfmKO_IFT$rank <- seq_len(nrow(dfmKO_IFT))
#Plot the enriched motifs
ggKO_IFT <- ggplot(dfmKO_IFT, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 3) +
  ggrepel::geom_label_repel(
    data = dfmKO_IFT[rev(seq_len(6)), ], aes(x = rank, y = mlog10Padj, label = TF), 
    size = 4,
    nudge_x = 1,
    nudge_y = 0,
    color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(P-adj) Motif Enrichment") + 
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet")) + theme(axis.title = element_text(size = 12)) + theme(axis.text = element_text(size = 12)) 
plotPDF(ggKO_IFT, name = "Fig6A_MotifEnrich_E85_IFTCMs_GainedDARs_vierstra", width = 7, height = 7, ArchRProj = proj_Mef2c_v13_E85_subset, addDOC = F)
