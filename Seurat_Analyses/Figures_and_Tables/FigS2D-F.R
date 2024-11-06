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

#Make DotPlots for Marker Genes

#E7.75
#Create cluster labels for dot plot
Idents(mef2c_v13_E775) <- "SCT_snn_res.0.5"
new.cluster.ids <- c("C0-PostM", "C1-Ect1", "C2-VE", "C3-Ect2/PS",
                     "C4-Ect3", "C5-Blood", "C6-Allantois", "C7-AntM", "C8-ExM",
                     "C9-Endo", "C10-ExEct", "C11-Ect4", "C12-DE",
                     "C13-CPs/CMs", "C14-AxM", "C15-T/P", "C16"
)
names(new.cluster.ids) <- levels(mef2c_v13_E775$SCT_snn_res.0.5)
mef2c_v13_E775 <- RenameIdents(mef2c_v13_E775, new.cluster.ids)
#Plot and save dot plot of key markers 
features1 = c(
  "Pdgfra", "Lef1",                     #Posterior Mesoderm
  "Pax2", "Dcc",                        #Ectoderm 1
  "Ttr","Cubn",                         #Visceral Endoderm
  "T", "Wnt3a",                         #Ectoderm 2 & PS
  "Npas3", "Nrcam",                     #Ectoderm 3
  "Hba-a1", "Hba-a2", "Gata1",          #Blood
  "Tbx4",                               #Allantois
  "Foxp2", "Foxc2", "Tbx1",             #Anterior Mesoderm
  "Hand1", "Bmp4", "Ahnak",             #Extraembryonic Mesoderm
  "Cdh5", "Kdr", "Nfatc1",              #Endothelium (and endocardium)
  "Tfap2a", "Wnt6",                     #Extraembryonic Ectoderm
  "Sox2", "Otx2",                       #Ectoderm 4
  "Sox17", "Trh", "Epcam",              #Definitive Endoderm
  "Nkx2-5", "Smarcd3", "Tbx5",          #Cardiac Progenitors/Myocytes
  "Mef2c","Tnnt2",
  "Nog", "Shh",                         #Axial Mesoderm
  "Gjb3", "Wnt7b",                      #Trophoblast/Placental
  "Bmper", "Grm8"                       #C16
)
DP_E775 <- DotPlot(mef2c_v13_E775, features = features1) + RotatedAxis()
ggsave('FigS2D_DotPlot_E775_cluster_markers.pdf', plot = DP_E775, device = 'pdf', width = 15, height = 5, dpi = 300)

#E8.5
#Create cluster labels for dot plot
Idents(mef2c_v13_E85_v3) <- "SCT_snn_res.0.4"
new.cluster.ids <- c("C0-NT1", "C1-SHF1/MixM", "C2-En1", "C3-En2",
                     "C4-SoM", "C5-CMs1", "C6-SHF2/PhM", "C7-NT2", 
                     "C8-NT3", "C9-SHF3/LPM", "C10-CMs2/CPs", "C11-NC",
                     "C12-NT4", "C13-NT5", "C14-Endo", "C15-En3",
                     "C16-Noto", "C17-YS")
names(new.cluster.ids) <- levels(mef2c_v13_E85_v3$SCT_snn_res.0.4)
mef2c_v13_E85_v3 <- RenameIdents(mef2c_v13_E85_v3, new.cluster.ids)
#Plot and save dot plot of key markers 
features2 = c(
  "Nrcam", "Fgf14", "Pax6",                #C0-NT1
  "Pdgfra", "Aldh1a2", "Isl1", "Fgf8",     #C1-SHF1/MixM
  "Epcam", "Cldn6", "Bcam",                #C2-En1
  "Spink1", "Spint2",                      #C3-En2
  "Meox1", "Fst", "Foxd1",                 #C4-SoM
  "Ttn", "Tnnt2", "Myocd",                 #C5-CMs1
  "Ebf1", "Tbx1", "Foxp2",                 #C6-SHF2/PhM
  "Rmst", "Pax5", "Dcc",                   #C7-NT2
  "Npas3", "Car10",                        #C8-NT3
  "Bmp4", "Hand1", "Dlk1",                 #C9-SHF3/LPM
  "Tbx5", "Mef2c",                         #C10-CMs2/CPs
  "Sox10", "Tfap2b", "Ets1",               #C11-NC
  "Nrg1", "Nrp2", "Opcml",                 #C12-NT4
  "Rfx4", "Ntn1", "Foxa2",                 #C13-NT5
  "Flt1", "Kdr", "Cdh5", "Nfatc1",         #C14-Endothelium
  "Ahnak", "Pmp22", "Wnt6",                #C15-En3
  "T", "Shh",                              #C16-Notochord
  "Ttr", "Afp", "Cubn"                     #C17-YS
)
DP_E85 <- DotPlot(mef2c_v13_E85_v3, features = features2) + RotatedAxis()
ggsave('FigS2E_DotPlot_E85_cluster_markers.pdf', plot = DP_E85, device = 'pdf', width = 18, height = 6, dpi = 300)

#E9
#Create cluster labels for dot plot
Idents(mef2c_v13_E9) <- "SCT_snn_res.0.6"
new.cluster.ids <- c("C0-En1", "C1-MixM", "C2-NT1", "C3-SoM1", 
                     "C4-CMs", "C5-NT2", "C6-SHF", "C7-NT3",
                     "C8-Pe", "C9-NC", "C10-SpM", "C11-SoM2",
                     "C12-Endo", "C13-En2", "C14-En3", "C15-NT4",
                     "C16-Me", "C17-En4", "C18-Blood", "C19-En5",
                     "C20-YS", "C21-Noto"
)
names(new.cluster.ids) <- levels(mef2c_v13_E9$SCT_snn_res.0.6)
mef2c_v13_E9 <- RenameIdents(mef2c_v13_E9, new.cluster.ids)
#Plot and save dot plot of key markers 
features3 = c(
  "Epcam", "Spint2", "Cldn6",        #Endoderm 1
  "Pdgfra", "Dlk1", "Prrx1",         #Mixed Meso
  "Sox2", "Npas3", "Fgf14",          #Neural Tube 1
  "Meox1", "Fst",                    #Somitic Mesoderm 1
  "Ttn", "Tnnt2", "Myl7",            #CMs
  "Dcc", "Nlgn1",                    #Neural Tube 2
  "Isl1", "Tbx1", "Fgf8", "Fgf10",   #SHF
  "Nrcam", "Nrg3",                   #Neural Tube 3
  "Tbx18", "Wt1",                    #Proepicardium
  "Tfap2b", "Sox10", "Ets1",         #Neural Crest
  "Dach2", "Foxp2", "Foxf1",         #Splanchnic Meso
  "Sim1", "Vegfc",                   #Somitic Mesoderm 2
  "Cdh5", "Kdr", "Flt1",             #Endothelium
  "Wnt6", "Pdgfc",                   #Endoderm 2
  "Cdh6", "Krt8",                    #Endoderm 3
  "Ntn1", "Slit1",                   #Neural Tube 4
  "Postn", "Col1a2", "Col3a1",       #Mesenchyme
  "Epha4", "Sox9",                   #Endoderm 4
  "Hba-a1", "Hba-a2", "Hbb-y",       #Blood
  "Ahnak", "Podxl",                  #Endoderm 5
  "Ttr", "Afp", "Apoa1",             #Yolk Sac
  "T", "Shh"                         #Notochord
)
DP_E9 <- DotPlot(mef2c_v13_E9, features = features3) + RotatedAxis()
ggsave('FigS2F_DotPlot_E9_cluster_markers.pdf', plot = DP_E9, device = 'pdf', width = 19, height = 7, dpi = 300)

