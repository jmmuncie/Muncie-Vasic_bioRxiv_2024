#Load Libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(magrittr)
library(XLConnect)
library(data.table)
library(harmony)
library(sctransform)

#Set working directory
setwd("~/Desktop/Working_Directory")

#Load Seurat Objects
#E7.75 Subset object
load("~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Objects/mef2c_v13_E775_subset_scTrans.Robj")
#E8.5 Subset object
load("~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Objects/mef2c_v13_E85_scTrans_v3_harmony_subset.Robj")
#E9 Subset object
load("~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Objects/mef2c_v13_E9_subset_scTrans.Robj")


#Run FindAllMarkers and Output Lists of Cluster Marker Genes

#E7.75
#Should return markers for 21 clusters, corresponding to res 2.0
Idents(mef2c_v13_E775_subset) <- "SCT_snn_res.2"
mef2c_E775_subset.markers <- FindAllMarkers(mef2c_v13_E775_subset, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
mef2c_E775_subset.markers <- subset(mef2c_E775_subset.markers, p_val_adj < 0.05)
mef2c_E775_subset.markers <- mef2c_E775_subset.markers[order(mef2c_E775_subset.markers$avg_log2FC, decreasing = TRUE), ]
#Save to spreadsheet
fwrite(mef2c_E775_subset.markers, row.names = TRUE, file = "TableS2A_E775_subset_markers.csv")

#E8.5
#Should return markers for 17 clusters, corresponding to res 1
Idents(mef2c_v13_E85_v3_subset) <- "SCT_snn_res.1"
mef2c_E85_sub.markers <- FindAllMarkers(mef2c_v13_E85_v3_subset, assay = "SCT", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.6)
mef2c_E85_sub.markers <- subset(mef2c_E85_sub.markers, p_val_adj < 0.05)
mef2c_E85_sub.markers <- mef2c_E85_sub.markers[order(mef2c_E85_sub.markers$avg_log2FC, decreasing = TRUE), ]
#Save to spreadsheet
fwrite(mef2c_E85_sub.markers, row.names = TRUE, file = "TableS2B_E85_subset_markers.csv")

#E9
#Should return markers for 12 clusters, corresponding to res 0.6
Idents(mef2c_v13_E9_subset) <- "SCT_snn_res.0.6"
mef2c_E9_subset.markers <- FindAllMarkers(mef2c_v13_E9_subset, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.7)
mef2c_E9_subset.markers <- subset(mef2c_E9_subset.markers, p_val_adj < 0.05)
mef2c_E9_subset.markers <- mef2c_E9_subset.markers[order(mef2c_E9_subset.markers$avg_log2FC, decreasing = TRUE), ]
#Save to spreadsheet
fwrite(mef2c_E9_subset.markers, row.names = TRUE, file = "TableS2C_E9_subset_markers.csv")