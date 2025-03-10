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


#Plot UMAPs

#E7.75
Idents(mef2c_v13_E775) <- "cell_type"
FigS2A_UMAP_E775_full_cell_types <- DimPlot(object = mef2c_v13_E775, 
                                            reduction = "umap", 
                                            label = T, repel = F, 
                                            label.size = 3.5, pt.size = 0.1
                                            ) + NoLegend()
FigS2A_UMAP_E775_full_sample_labels <- DimPlot(object = mef2c_v13_E775, 
                                               reduction = "umap", 
                                               group.by = "Sample_Name", 
                                               cols = c("lightskyblue","turquoise1","blue","darkblue"), 
                                               shuffle = T, seed = 1, pt.size = 0.1
                                               ) + NoLegend()
ggsave('FigS2A_UMAP_E775_full_cell_types.pdf', plot = FigS2A_UMAP_E775_full_cell_types + labs(title = 'Mef2c E7.75: 21237 cells'), device = 'pdf', width = 5, height = 5, dpi = 300)
ggsave('FigS2A_UMAP_E775_full_sample_labels.pdf', plot = FigS2A_UMAP_E775_full_sample_labels + labs(title = 'Mef2c E7.75: 21237 cells'), device = 'pdf', width = 5, height = 5, dpi = 300)

#E8.5
Idents(mef2c_v13_E85_v3) <- "harmony_cell_type"
FigS2B_UMAP_E85_full_cell_types <- DimPlot(object = mef2c_v13_E85_v3, 
                                           reduction = "umap", 
                                           label = T, repel = T, 
                                           label.size = 4, pt.size = 0.1
                                           ) + NoLegend()
FigS2B_UMAP_E85_full_sample_labels <- DimPlot(object = mef2c_v13_E85_v3, 
                                      reduction = "umap", 
                                      group.by = "Sample_Name", 
                                      cols = c("#F3E5F5", "#D1C4E9", "#9C27B0", "#4527A0"), 
                                      shuffle = T, seed = 1, pt.size = 0.1
                                      ) + NoLegend()
ggsave('FigS2B_UMAP_E85_full_cell_types.pdf', plot = FigS2B_UMAP_E85_full_cell_types + labs(title = 'Mef2c E8.5: 28373 cells'), device = 'pdf', width = 5, height = 5, dpi = 300)
ggsave('FigS2B_UMAP_E85_full_sample_labels.pdf', plot = FigS2B_UMAP_E85_full_sample_labels + labs(title = 'Mef2c E8.5: 28373 cells'), device = 'pdf', width = 5, height = 5, dpi = 300)

#E9
Idents(mef2c_v13_E9) <- "cell_type"
FigS2C_UMAP_E9_full_cell_types <- DimPlot(object = mef2c_v13_E9, 
                                          reduction = "umap", 
                                          label = T, repel = F, 
                                          label.size = 4, pt.size = 0.1
                                          ) + NoLegend()
FigS2C_UMAP_E9_full_sample_labels <- DimPlot(object = mef2c_v13_E9, 
                                             reduction = "umap", 
                                             group.by = "Sample_Name", 
                                             cols = c("pink", "plum1", "red", "red4"), 
                                             shuffle = T, seed = 1, pt.size = 0.1
                                             ) + NoLegend()
ggsave('FigS2C_UMAP_E9_full_cell_types.pdf', plot = FigS2C_UMAP_E9_full_cell_types + labs(title = 'Mef2c E9: 47482 cells'), device = 'pdf', width = 5, height = 5, dpi = 300)
ggsave('FigS2C_UMAP_E9_full_sample_labels.pdf', plot = FigS2C_UMAP_E9_full_sample_labels + labs(title = 'Mef2c E9: 47482 cells'), device = 'pdf', width = 5, height = 5, dpi = 300)
