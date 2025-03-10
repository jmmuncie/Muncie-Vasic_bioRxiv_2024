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


#Make DotPlots for Genes of Interest

#E7.75
mef2c_v13_E775_subset$cell_type_pool_x_genotype <- factor(mef2c_v13_E775_subset$cell_type_pool_x_genotype, 
                                                          levels = c("CMs/FHF_WT", "CMs/FHF_KO", "SHF_WT","SHF_KO", 
                                                                     "JCF_WT", "JCF_KO", "CrM_WT", "CrM_KO",
                                                                     "ExM_WT","ExM_KO", "SoM_WT", "SoM_KO", 
                                                                     "LPM_WT", "LPM_KO", "NMPs_WT", "NMPs_KO",
                                                                     "PrxM_WT", "PrxM_KO", "KPs_WT", "KPs_KO",
                                                                     "HSCs_WT", "HSCs_KO")
                                                          )
Idents(mef2c_v13_E775_subset) <- "cell_type_pool_x_genotype"
Fig2G_DotPlot_E775 <- DotPlot(mef2c_v13_E775_subset, 
                              features = c("Tnni1", "Myl7", "Myh7", "Cacna1c", "Myom1", "Actc1", "Myh6", "Ttn"), 
                              idents = c("CMs/FHF_WT", "CMs/FHF_KO", "SHF_WT", "SHF_KO", "JCF_WT", "JCF_KO")
                              ) + coord_flip()
ggsave('Fig2G_DotPlot_E775.pdf', plot = Fig2G_DotPlot_E775, device = 'pdf', width = 6.5, height = 4.5, dpi = 300)

#E8.5
mef2c_v13_E85_v3_subset$cell_type_pool_x_genotype <- factor(mef2c_v13_E85_v3_subset$cell_type_pool_x_genotype, 
                                                            levels = c("pSHF_WT", "pSHF_KO", "aSHF_WT", "aSHF_KO",
                                                                       "CMs-IFT_WT", "CMs-IFT_KO", "CMs-V_WT", "CMs-V_KO",   
                                                                       "CMs-OFT_WT", "CMs-OFT_KO", "PhM_WT", "PhM_KO", 
                                                                       "LPM_WT", "LPM_KO", "PostM_WT", "PostM_KO",
                                                                       "MixM_WT", "MixM_KO", "C16_WT", "C16_KO") 
                                                            )
Idents(mef2c_v13_E85_v3_subset) <- "cell_type_pool_x_genotype"
Fig2H_DotPlot_E85 <- DotPlot(mef2c_v13_E85_v3_harmony_subset, 
                             features = c("Tdgf1", "Wnt2", "Gata4", "Tbx5", "Nkx2-5", "Nppa", "Ttn", "Tnnt2"), 
                             idents = c("CMs-IFT_WT", "CMs-IFT_KO", "CMs-V_WT", "CMs-V_KO", "CMs-OFT_WT", "CMs-OFT_KO")
                             ) + coord_flip()
ggsave('Fig2H_DotPlot_E85.pdf', plot = Fig2H_DotPlot_E85, device = 'pdf', width = 6.5, height = 4.5, dpi = 300)

#E9
mef2c_v13_E9_subset$cell_type_pool_x_genotype <- factor(mef2c_v13_E9_subset$cell_type_pool_x_genotype, 
                                                        levels = c("SHF_WT","SHF_KO",
                                                                   "Pe_WT", "Pe_KO", 
                                                                   "VP_WT", "VP_KO",
                                                                   "CMs-A_WT", "CMs-A_KO",
                                                                   "CMs-AVC_WT", "CMs-AVC_KO",
                                                                   "CMs-V_WT", "CMs-V_KO",
                                                                   "CMs-OFT_WT","CMs-OFT_KO", 
                                                                   "PhM_WT", "PhM_KO",
                                                                   "C11_WT", "C11_KO")
                                                        )
Idents(mef2c_v13_E9_subset) <- "cell_type_pool_x_genotype"
Fig2I_DotPlot_E9 <- DotPlot(mef2c_v13_E9_subset, 
                            features = c("Tdgf1", "Wnt2", "Gata4", "Tbx5", "Nkx2-5", "Nppa", "Ttn", "Tnnt2"), 
                            idents = c("CMs-A_WT", "CMs-A_KO","CMs-AVC_WT","CMs-AVC_KO", "CMs-V_WT", "CMs-V_KO", "CMs-OFT_WT", "CMs-OFT_KO")
                            ) + coord_flip()
ggsave('Fig2I_DotPlot_E9.pdf', plot = Fig2I_DotPlot_E9, device = 'pdf', width = 7, height = 4.5, dpi = 300)
