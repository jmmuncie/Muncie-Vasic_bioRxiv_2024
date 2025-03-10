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

#E7.75 Ridge Plot
Idents(mef2c_v13_E775_subset) <- "cell_type_pool_x_genotype"
rp1 <- RidgePlot(object = mef2c_v13_E775_subset,
                 features = "Mef2c",
                 idents = "CMs/FHF_WT"
                 )
ggsave('FigS1A_RidgePlot_E775_CMsFHF_Mef2c.pdf', plot = rp1, device = 'pdf', width = 7, height = 4, dpi = 300)

#E8.5 Ridge Plot
Idents(mef2c_v13_E85_v3_subset) <- "cell_type_pool_x_genotype"
rp2 <- RidgePlot(object = mef2c_v13_E85_v3_subset,
                 features = "Mef2c",
                 idents = c("CMs-IFT_WT", "CMs-V_WT", "CMs-OFT_WT")
                 )
ggsave('FigS1A_RidgePlot_E85_IFT_V_OFT_Mef2c.pdf', plot = rp2, device = 'pdf', width = 7, height = 4, dpi = 300)

#E9 Ridge Plot
Idents(mef2c_v13_E9_subset) <- "cell_type_pool_x_genotype"
rp3 <- RidgePlot(object = mef2c_v13_E9_subset,
                 features = "Mef2c",
                 idents = c("CMs-A_WT", "CMs-AVC_WT", "CMs-V_WT", "CMs-OFT_WT")
                 )
ggsave('FigS1A_RidgePlot_E9_A_AVC_V_OFT_Mef2c.pdf', plot = rp3, device = 'pdf', width = 7, height = 4, dpi = 300)
