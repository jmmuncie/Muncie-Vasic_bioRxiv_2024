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


#Make DotPlots for Marker Genes

#E7.75
#Create cluster labels for dot plot
Idents(mef2c_v13_E775_subset) <- "SCT_snn_res.2"
new.cluster.ids <- c("C0-CrM1", "C1-SoM1", "C2-ExM1", "C3-NMPs",
                     "C4-ExM2", "C5-PrxM", "C6-CrM2", "C7-SoM2",
                     "C8-ExM3", "C9-ExM4", "C10-LPM", "C11-SHF1",
                     "C12-ExM5", "C13-ExM6", "C14-JCF", "C15-SoM3",
                     "C16-CMs/FHF", "C17-HSCs", "C18-SHF2", "C19-CrM3", "C20-KPs"
)
names(new.cluster.ids) <- levels(mef2c_v13_E775_subset$SCT_snn_res.2)
mef2c_v13_E775_subset <- RenameIdents(mef2c_v13_E775_subset, new.cluster.ids)
#Plot and save dot plot of key markers 
features1 = c(
  "Cped1", "Sobp",                      #CrM1
  "Fst","Meox1",                        #SoM1
  "Hand1", "Ahnak", "Bmp4",             #ExM1
  "T", "Nkx1-2", "Sox2",                #NMPs 
  "Tdo2", "Bnc2",                       #ExM2
  "Tbx1", "Foxp2",                      #PrxM
  "Nrxn3", "Pitx2",                     #CrM2
  "Notch1", "Aldh1a2",                  #SoM2
  "Adamts17", "Msx1",                   #ExM3
  "Smoc2", "Tek",                       #ExM4
  "Epha7", "Gata4", "Gata6",            #LPM
  "Isl1", "Fgf8", "Fgf10",              #SHF1
  "Morc4", "Ctnnd2",                    #ExM5
  "Pmp22", "Tagln",                     #ExM6
  "Mab21l2",                            #JCF
  "Tcf15", "Cer1",                      #SoM3
  "Tnnt2", "Nkx2-5", "Tbx5",            #CMs/FHF
  "Runx1", "Itga4", "Cd44",             #HSCs
  "Mef2c", "Hand2",                     #SHF2
  "Otx2", "Eya1",                       #CrM3
  "Slc2a2", "Cubn", "Amn"               #KPs
)
DP_E775_subset <- DotPlot(mef2c_v13_E775_subset, features = features1) + RotatedAxis()
ggsave('FigS2G_DotPlot_E775_Subset_cluster_markers.pdf', plot = DP_E775_subset, device = 'pdf', width = 17, height = 7, dpi = 300)

#E8.5
#Create cluster labels for dot plot
Idents(mef2c_v13_E85_v3_subset) <- "SCT_snn_res.1"
new.cluster.ids <- c("C0-CMs-V", "C1-pSHF", "C2-MixM1", "C3-PhM1",
                     "C4-MixM2", "C5-CMs-IFT1", "C6-PhM2", "C7-CMs-OFT",
                     "C8-LPM1", "C9-PhM3", "C10-CMs-IFT2", "C11-aSHF",
                     "C12-LPM2", "C13-CMs-IFT3", "C14-PostM", "C15-MixM3",
                     "C16")
names(new.cluster.ids) <- levels(mef2c_v13_E85_v3_subset$SCT_snn_res.1)
mef2c_v13_E85_v3_subset <- RenameIdents(mef2c_v13_E85_v3_subset, new.cluster.ids)
#Plot and save dot plot of key markers 
features2 = c(
  "Ttn", "Nkx2-5", "Myl2", "Irx4",         #C0-CMs-V
  "Aldh1a2", "Hoxb1",                      #C1-pSHF
  "Vegfc", "Pdgfra",                       #C2-MixM1
  "Fst", "Ebf1", "Tbx1",                   #C3-PhM1
  #C4-MixM2 - no DEGs show up
  "Tbx5", "Wnt2",                          #C5-CMs-IFT1
  "Cdh11", "Foxd1",                        #C6-PhM2
  "Tdgf1", "Rgs5", "Nrp2",                 #C7-CMs-OFT
  "Dlk1", "Bmp4", "Hand1",                 #C8-LPM1
  "Nrg1", "Eda",                           #C9-PhM3
  "Mef2c", "Gata4",                        #C10-CMs-IFT2
  "Isl1", "Fgf8", "Fgf10",                 #C11-aSHF
  "Pmp22", "Nid1",                         #C12-LPM2
  "Tbx20", "Vsnl1",                        #C13-CMs-IFT3
  "Lef1", "Bmp5",                          #C14-PostM
  "Tnc", "Twist1",                         #C15-MixM3
  "Dcc", "Fgf14", "Nrcam"                  #C16
)
DP_E85_subset <- DotPlot(mef2c_v13_E85_v3_subset, features = features2) + RotatedAxis()
ggsave('FigS2H_DotPlot_E85_Subset_cluster_markers.pdf', plot = DP_E85_subset, device = 'pdf', width = 18, height = 6, dpi = 300)

#E9
#Create cluster labels for dot plot
Idents(mef2c_v13_E9_subset) <- "SCT_snn_res.0.6"
new.cluster.ids <- c("C0-SHF1", "C1-CMs-A1", "C2-CMs-V", "C3-VP1", "C4-CMs-AVC",
                     "C5-SHF2", "C6-CMs-OFT", "C7-PhM", "C8-CMs-A2", "C9-Pe",
                     "C10-VP2", "C11"
)
names(new.cluster.ids) <- levels(mef2c_v13_E9_subset$SCT_snn_res.0.6)
mef2c_v13_E9_subset <- RenameIdents(mef2c_v13_E9_subset, new.cluster.ids)
#Plot and save dot plot of key markers 
features3 = c(
  "Fgf10", "Nrg1", "Isl1",                         #SHF1
  "Angpt1", "Nr2f2", "Vsnl1",                      #Atrial CMs 1
  "Ttn", "Tnnt2", "Myl7",                          #CM markers
  "Myl2", "Irx4",                                  #Ventricular CMs
  "Wnt2", "Tbx5", "Nav3",                          #VP1
  "Rgs6", "Tbx2", "Rspo3",                         #AVC CMs
  "Cdh6", "Bmp5", "Mecom",                         #SHF2
  "Nrp2", "Tdgf1", "Sema3c", "Wnt11",              #OFT CMs
  "Ebf1", "Cdh11", "Tbx1",                         #PhM
  "Kcnq3", "Ankrd1",                               #Atrial CMs 2
  "Wt1", "Tbx18",                                  #Pe
  "Ror1", "Col1a1", "Podxl",                       #VP2     
  "Fgf14", "Nrcam", "Dcc"                          #C11
)
DP_E9_subset <- DotPlot(mef2c_v13_E9_subset, features = features3) + RotatedAxis()
ggsave('FigS2I_DotPlot_E9_Subset_cluster_markers.pdf', plot = DP_E9_subset, device = 'pdf', width = 15, height = 6, dpi = 300)


