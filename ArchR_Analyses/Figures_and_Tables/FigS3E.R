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


#Get DARs and Plot Motif Enrichments

#OFT-CMs
#Look at motifs enriched in peaks with -Log2FC in markerTest_CMs 
#i.e. the peaks that lose accessibility in MEF2C KO
markerTest_CMs_OFT <- getMarkerFeatures(
  ArchRProj = proj_Mef2c_v13_E85_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "CellTypeByGenotype",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "CMs_OFT_x_KO",
  bgdGroups = "CMs_OFT_x_WT"
  )
motifsWT_OFT <- peakAnnoEnrichment(
  seMarker = markerTest_CMs_OFT,
  ArchRProj = proj_Mef2c_v13_E85_subset,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.15 & Log2FC <= -0.5"
  )
#To prepare for plotting, create simple df containing
#the motif names, the corrected p-values, and the significance rank
dfmWT_OFT <- data.frame(TF = rownames(motifsWT_OFT), mlog10Padj = assay(motifsWT_OFT)[,1])
dfmWT_OFT <- dfmWT_OFT[order(dfmWT_OFT$mlog10Padj, decreasing = TRUE),]
dfmWT_OFT$rank <- seq_len(nrow(dfmWT_OFT))
#Plot the enriched motifs
ggWT_OFT <- ggplot(dfmWT_OFT, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 3) +
  ggrepel::geom_label_repel(
    data = dfmWT_OFT[rev(seq_len(6)), ], aes(x = rank, y = mlog10Padj, label = TF), 
    size = 4,
    nudge_x = 1,
    nudge_y = 0,
    color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(P-adj) Motif Enrichment") + 
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet")) + theme(axis.title = element_text(size = 12)) + theme(axis.text = element_text(size = 12)) 
plotPDF(ggWT_OFT, name = "FigS3E_MotifEnrich_E85_OFTCMs_LostDARs_vierstra", width = 7, height = 7, ArchRProj = proj_Mef2c_v13_E85_subset, addDOC = F)

#V-CMs
#Look at motifs enriched in peaks with -Log2FC in markerTest_CMs 
#i.e. the peaks that lose accessibility in MEF2C KO
markerTest_CMs_V <- getMarkerFeatures(
  ArchRProj = proj_Mef2c_v13_E85_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "CellTypeByGenotype",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "CMs_V_x_KO",
  bgdGroups = "CMs_V_x_WT"
  )
motifsWT_V <- peakAnnoEnrichment(
  seMarker = markerTest_CMs_V,
  ArchRProj = proj_Mef2c_v13_E85_subset,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.15 & Log2FC <= -0.5"
  )
#To prepare for plotting, create simple df containg
#the motif names, the corrected p-values, and the significance rank
dfmWT_V <- data.frame(TF = rownames(motifsWT_V), mlog10Padj = assay(motifsWT_V)[,1])
dfmWT_V <- dfmWT_V[order(dfmWT_V$mlog10Padj, decreasing = TRUE),]
dfmWT_V$rank <- seq_len(nrow(dfmWT_V))
#Plot the enriched motifs
ggWT_V <- ggplot(dfmWT_V, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 3) +
  ggrepel::geom_label_repel(
    data = dfmWT_V[rev(seq_len(6)), ], aes(x = rank, y = mlog10Padj, label = TF), 
    size = 4,
    nudge_x = 1,
    nudge_y = 0,
    color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(P-adj) Motif Enrichment") + 
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet")) + theme(axis.title = element_text(size = 12)) + theme(axis.text = element_text(size = 12)) 
plotPDF(ggWT_V, name = "FigS3E_MotifEnrich_E85_VCMs_LostDARs_vierstra", width = 7, height = 7, ArchRProj = proj_Mef2c_v13_E85_subset, addDOC = F)

#IFT-CMs
#Look at motifs enriched in peaks with -Log2FC in markerTest_CMs 
#i.e. the peaks that lose accessibility in MEF2C KO
markerTest_CMs_IFT <- getMarkerFeatures(
  ArchRProj = proj_Mef2c_v13_E85_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "CellTypeByGenotype",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "CMs_IFT_x_KO",
  bgdGroups = "CMs_IFT_x_WT"
  )
motifsWT_IFT <- peakAnnoEnrichment(
  seMarker = markerTest_CMs_IFT,
  ArchRProj = proj_Mef2c_v13_E85_subset,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.15 & Log2FC <= -0.5"
)
#To prepare for plotting, create simple df containg
#the motif names, the corrected p-values, and the significance rank
dfmWT_IFT <- data.frame(TF = rownames(motifsWT_IFT), mlog10Padj = assay(motifsWT_IFT)[,1])
dfmWT_IFT <- dfmWT_IFT[order(dfmWT_IFT$mlog10Padj, decreasing = TRUE),]
dfmWT_IFT$rank <- seq_len(nrow(dfmWT_IFT))
#Plot the enriched motifs
ggWT_IFT <- ggplot(dfmWT_IFT, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 3) +
  ggrepel::geom_label_repel(
    data = dfmWT_IFT[rev(seq_len(6)), ], aes(x = rank, y = mlog10Padj, label = TF), 
    size = 4,
    nudge_x = 1,
    nudge_y = 0,
    color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(P-adj) Motif Enrichment") + 
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet")) + theme(axis.title = element_text(size = 12)) + theme(axis.text = element_text(size = 12)) 
plotPDF(ggWT_IFT, name = "FigS3E_MotifEnrich_E85_IFTCMs_LostDARs_vierstra", width = 7, height = 7, ArchRProj = proj_Mef2c_v13_E85_subset, addDOC = F)
