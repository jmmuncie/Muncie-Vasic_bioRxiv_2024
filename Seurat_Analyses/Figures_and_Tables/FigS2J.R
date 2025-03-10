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

#Make DotPlots for Top Differentially Expressed Genes

#E7.75
#Top 10 up/down for CMs/FHF plotted for CMs/FHF, SHF, and JCF
mef2c_v13_E775_subset$cell_type_pool_x_genotype <- factor(mef2c_v13_E775_subset$cell_type_pool_x_genotype, 
                                                          levels = c("CMs/FHF_WT", "CMs/FHF_KO", "SHF_WT","SHF_KO", 
                                                                     "JCF_WT", "JCF_KO", "CrM_WT", "CrM_KO",
                                                                     "ExM_WT","ExM_KO", "SoM_WT", "SoM_KO", 
                                                                     "LPM_WT", "LPM_KO", "NMPs_WT", "NMPs_KO",
                                                                     "PrxM_WT", "PrxM_KO", "KPs_WT", "KPs_KO",
                                                                     "HSCs_WT", "HSCs_KO")
                                                          )
Idents(mef2c_v13_E775_subset) <- "cell_type_pool_x_genotype"
FigS2J_DP_E775 <- DotPlot(mef2c_v13_E775_subset, 
                          features = c("Apoc2", "Ttr", "Sfrp5", "Spink1", "Rbp4", "Peg10", "Adam12", "Igf2", "Apoa1", "Sox4",
                                       "Epha7", "Mid1", "Nebl", "Myom1", "Actc1", "Unc13c", "Tnnc1", "Erbb4", "Myh6", "Ttn"), 
                          idents = c("CMs/FHF_WT", "CMs/FHF_KO", "SHF_WT", "SHF_KO", "JCF_WT", "JCF_KO")
                          ) + coord_flip()
ggsave('FigS2J_DotPlot_E775.pdf', plot = FigS2J_DP_E775, device = 'pdf', width = 6, height = 8, dpi = 300)

#E8.5
#Top 4 up/down for each group: IFT, V, OFT CMs
mef2c_v13_E85_v3_subset$cell_type_pool_x_genotype <- factor(mef2c_v13_E85_v3_subset$cell_type_pool_x_genotype, 
                                                            levels = c("pSHF_WT", "pSHF_KO", "aSHF_WT", "aSHF_KO",
                                                                       "CMs-IFT_WT", "CMs-IFT_KO", "CMs-V_WT", "CMs-V_KO",   
                                                                       "CMs-OFT_WT", "CMs-OFT_KO", "PhM_WT", "PhM_KO", 
                                                                       "LPM_WT", "LPM_KO", "PostM_WT", "PostM_KO",
                                                                       "MixM_WT", "MixM_KO", "C16_WT", "C16_KO") 
                                                            )
Idents(mef2c_v13_E85_v3_subset) <- "cell_type_pool_x_genotype"
FigS2J_DP_E85 <- DotPlot(mef2c_v13_E85_v3_subset, 
                         features = c("Igfbp5", "Nfia", "Tenm2", "Cacna2d2", "Tenm3", "Robo2", "Igf2", "Sox4", "Dlk1", "Nnat", "Ccbe1", "Scd2",
                                      "Myl7", "Myom1", "Palld", "Mid1", "Hspb1", "Actn2", "Myh7", "Myl3", "Pde4d", "Ttn", "Tnnc1", "Myh6"), 
                         idents = c("CMs-IFT_WT", "CMs-IFT_KO", "CMs-V_WT", "CMs-V_KO", "CMs-OFT_WT", "CMs-OFT_KO")
                         ) + coord_flip()
ggsave('FigS2J_DotPlot_E85.pdf', plot = FigS2J_DP_E85, device = 'pdf', width = 6, height = 9, dpi = 300)

#E9
#Top 3 up/down for each group: A, AVC, V, OFT CMs
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
FigS2J_DP_E9 <- DotPlot(mef2c_v13_E9_subset, 
                        features = c("Jarid2", "Dpp6", "Nfia", "Tenm3", "Edil3", "Robo2", "Kcnh7", "Igf2r", "Cacna1d", "Slc2a1", "Magi1", "Foxp1",
                                     "Mid1", "Tnnc1", "Myl3", "Ccnd2", "Myh7", "Nppa", "Cacna1c", "Actn2", "Myl2", "Ctnna3", "Myl7", "Rpl29"), 
                        idents = c("CMs-A_WT", "CMs-A_KO","CMs-AVC_WT","CMs-AVC_KO", "CMs-V_WT", "CMs-V_KO", "CMs-OFT_WT", "CMs-OFT_KO")
                        ) + coord_flip()
ggsave('FigS2J_DotPlot_E9.pdf', plot = FigS2J_DP_E9, device = 'pdf', width = 7, height = 9, dpi = 300)

