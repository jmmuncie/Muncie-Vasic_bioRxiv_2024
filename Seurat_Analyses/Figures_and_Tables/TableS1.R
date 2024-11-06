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
#E7.75 Full Dataset object
load("~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Objects/mef2c_v13_E775_scTrans.Robj")
#E8.5 Full Dataset object
load("~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Objects/mef2c_v13_E85_scTrans_v3_harmony.Robj")
#E9 Full Dataset object
load("~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Objects/mef2c_v13_E9_scTrans.Robj")


#Run FindAllMarkers and Output Lists of Cluster Marker Genes

#E7.75
#Should return markers for 17 clusters, corresponding to res 0.5
Idents(mef2c_v13_E775) <- "SCT_snn_res.0.5"
mef2c_E775.markers <- FindAllMarkers(mef2c_v13_E775, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
mef2c_E775.markers <- subset(mef2c_E775.markers, p_val_adj < 0.05)
mef2c_E775.markers <- mef2c_E775.markers[order(mef2c_E775.markers$avg_log2FC, decreasing = TRUE), ]
#Save to spreadsheet
fwrite(mef2c_E775.markers, row.names = TRUE, file = "TableS1A_E775_markers.csv")

#E8.5
#Should return markers for 18 clusters, corresponding to res 0.4
Idents(mef2c_v13_E85_v3) <- "SCT_snn_res.0.4"
mef2c_E85.markers <- FindAllMarkers(mef2c_v13_E85_v3, assay = "SCT", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.6)
mef2c_E85.markers <- subset(mef2c_E85.markers, p_val_adj < 0.05)
mef2c_E85.markers <- mef2c_E85.markers[order(mef2c_E85.markers$avg_log2FC, decreasing = TRUE), ]
#Save to spreadsheet
fwrite(mef2c_E85.markers, row.names = TRUE, file = "TableS1B_E85_markers.csv")

#E9
#Should return markers for 22 clusters, corresponding to res 0.6
Idents(mef2c_v13_E9) <- "SCT_snn_res.0.6"
mef2c_E9.markers <- FindAllMarkers(mef2c_v13_E9, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.6)
mef2c_E9.markers <- subset(mef2c_E9.markers, p_val_adj < 0.05)
mef2c_E9.markers <- mef2c_E9.markers[order(mef2c_E9.markers$avg_log2FC, decreasing = TRUE), ]
#Save to spreadsheet
fwrite(mef2c_E9.markers, row.names = TRUE, file = "TableS1C_E9_markers.csv")