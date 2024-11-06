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


#Plot UMAPs

#E7.75
Idents(mef2c_v13_E775_subset) <- "cell_type_subset"
Fig2A_UMAP_E775_subset_cell_types <- DimPlot(object = mef2c_v13_E775_subset, 
                                             reduction = "umap", 
                                             label = T, repel = T, 
                                             label.size = 4, pt.size = 0.5
                                             ) + NoLegend()
Fig2A_UMAP_E775_subset_sample_labels <- DimPlot(object = mef2c_v13_E775_subset, 
                                                reduction = "umap", 
                                                group.by = "Sample_Name", 
                                                cols = c("lightskyblue","turquoise1","blue","darkblue"), 
                                                shuffle = T, seed = 1, pt.size = 0.5
                                                ) + NoLegend()
ggsave('Fig2A_UMAP_E775_subset_cell_types.pdf', plot = Fig2A_UMAP_E775_subset_cell_types + labs(title = 'Mef2c E7.75 Subset: 3356 cells'), device = 'pdf', width = 5, height = 5, dpi = 300)
ggsave('Fig2A_UMAP_E775_subset_sample_labels.pdf', plot = Fig2A_UMAP_E775_subset_sample_labels + labs(title = 'Mef2c E7.75 Subset: 3356 cells'), device = 'pdf', width = 5, height = 5, dpi = 300)

#E8.5
Idents(mef2c_v13_E85_v3_subset) <- "harmony_cell_type_subset"
Fig2B_UMAP_E85_subset_cell_types <- DimPlot(object = mef2c_v13_E85_v3_subset, 
                                            reduction = "umap", 
                                            label = T, repel = T, 
                                            label.size = 4, pt.size = 0.1
                                            ) + NoLegend()
Fig2B_UMAP_E85_subset_sample_labels <- DimPlot(object = mef2c_v13_E85_v3_subset, 
                                               reduction = "umap", 
                                               group.by = "Sample_Name", 
                                               cols = c("#F3E5F5", "#D1C4E9", "#9C27B0", "#4527A0"), 
                                               shuffle = T, seed = 1, pt.size = 0.1
                                               ) + NoLegend()
ggsave('Fig2B_UMAP_E85_subset_cell_types.pdf', plot = Fig2B_UMAP_E85_subset_cell_types + labs(title = 'Mef2c E8.5 Subset: 9240 cells'), device = 'pdf', width = 5, height = 5, dpi = 300)
ggsave('Fig2B_UMAP_E85_subset_sample_labels.pdf', plot = Fig2B_UMAP_E85_subset_sample_labels + labs(title = 'Mef2c E8.5 Subset: 9240 cells'), device = 'pdf', width = 5, height = 5, dpi = 300)

#E9
Idents(mef2c_v13_E9_subset) <- "cell_type_subset"
Fig2C_UMAP_E9_subset_cell_types <- DimPlot(object = mef2c_v13_E9_subset, 
                                           reduction = "umap", 
                                           label = T, repel = T, 
                                           label.size = 4, pt.size = 0.5
                                           ) + NoLegend()
Fig2C_UMAP_E9_subset_sample_labels <- DimPlot(object = mef2c_v13_E9_subset, 
                                              reduction = "umap", 
                                              group.by = "Sample_Name", 
                                              cols = c("pink", "plum1", "red", "red4"), 
                                              shuffle = T, seed = 1, pt.size = 0.5
                                              ) + NoLegend()
ggsave('Fig2C_UMAP_E9_subset_cell_types.pdf', plot = Fig2C_UMAP_E9_subset_cell_types + labs(title = 'Mef2c E9 Subset: 9259 cells'), device = 'pdf', width = 5, height = 5, dpi = 300)
ggsave('Fig2C_UMAP_E9_subset_sample_labels.pdf', plot = Fig2C_UMAP_E9_subset_sample_labels + labs(title = 'Mef2c E9 Subset: 9259 cells'), device = 'pdf', width = 5, height = 5, dpi = 300)
