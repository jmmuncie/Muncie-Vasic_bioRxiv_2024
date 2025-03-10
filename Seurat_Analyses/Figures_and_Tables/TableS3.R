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


#Run FindMarkers and Output DEG Lists


#E7.75

#Test for DEG between KO and WT cells in CMs/FHF cluster
Idents(mef2c_v13_E775_subset) <- "cell_type_subset"
markers_mef2c_E775_KOvWT_CMFHF <- FindMarkers(mef2c_v13_E775_subset, ident.1 = "KO", ident.2 = "WT", group.by = "Genotype", subset.ident = "CMs/FHF")
markers_mef2c_E775_KOvWT_CMFHF <- subset(markers_mef2c_E775_KOvWT_CMFHF, p_val_adj < 0.05)
markers_mef2c_E775_KOvWT_CMFHF <- markers_mef2c_E775_KOvWT_CMFHF[order(markers_mef2c_E775_KOvWT_CMFHF$avg_log2FC, decreasing = TRUE), ]
names(markers_mef2c_E775_KOvWT_CMFHF)[3] <- "pct.KO"
names(markers_mef2c_E775_KOvWT_CMFHF)[4] <- "pct.WT"
fwrite(markers_mef2c_E775_KOvWT_CMFHF, row.names = TRUE, file = "TableS3A_E775_CMFHF_KOvWT_DEGs.csv")

#Test for DEG between KO and WT cells in SHF clusters
Idents(mef2c_v13_E775_subset) <- "cell_type_subset"
markers_mef2c_E775_KOvWT_SHF <- FindMarkers(mef2c_v13_E775_subset, ident.1 = "KO", ident.2 = "WT", group.by = "Genotype", subset.ident = c("SHF1", "SHF2"))
markers_mef2c_E775_KOvWT_SHF <- subset(markers_mef2c_E775_KOvWT_SHF, p_val_adj < 0.05)
markers_mef2c_E775_KOvWT_SHF <- markers_mef2c_E775_KOvWT_SHF[order(markers_mef2c_E775_KOvWT_SHF$avg_log2FC, decreasing = TRUE), ]
names(markers_mef2c_E775_KOvWT_SHF)[3] <- "pct.KO"
names(markers_mef2c_E775_KOvWT_SHF)[4] <- "pct.WT"
fwrite(markers_mef2c_E775_KOvWT_SHF, row.names = TRUE, file = "TableS3B_E775_SHF_KOvWT_DEGs.csv")

#Test for DEG between KO and WT cells in JCF cluster
Idents(mef2c_v13_E775_subset) <- "cell_type_subset"
markers_mef2c_E775_KOvWT_JCF <- FindMarkers(mef2c_v13_E775_subset, ident.1 = "KO", ident.2 = "WT", group.by = "Genotype", subset.ident = "JCF")
markers_mef2c_E775_KOvWT_JCF <- subset(markers_mef2c_E775_KOvWT_JCF, p_val_adj < 0.05)
markers_mef2c_E775_KOvWT_JCF <- markers_mef2c_E775_KOvWT_JCF[order(markers_mef2c_E775_KOvWT_JCF$avg_log2FC, decreasing = TRUE), ]
names(markers_mef2c_E775_KOvWT_JCF)[3] <- "pct.KO"
names(markers_mef2c_E775_KOvWT_JCF)[4] <- "pct.WT"
fwrite(markers_mef2c_E775_KOvWT_JCF, row.names = TRUE, file = "TableS3C_E775_JCF_KOvWT_DEGs.csv")


#E8.5

#Test for DEG between KO and WT cells in IFT CMs clusters (3 clusters)
Idents(mef2c_v13_E85_v3_subset) <- "harmony_cell_type_subset"
markers_mef2c_E85_KOvWT_IFTCM <- FindMarkers(mef2c_v13_E85_v3_subset, asay = "SCT", ident.1 = "KO", ident.2 = "WT", group.by = "Genotype", subset.ident = c("CMs-IFT1", "CMs-IFT2", "CMs-IFT3"))
markers_mef2c_E85_KOvWT_IFTCM <- subset(markers_mef2c_E85_KOvWT_IFTCM, p_val_adj < 0.05)
markers_mef2c_E85_KOvWT_IFTCM <- markers_mef2c_E85_KOvWT_IFTCM[order(markers_mef2c_E85_KOvWT_IFTCM$avg_log2FC, decreasing = TRUE), ]
names(markers_mef2c_E85_KOvWT_IFTCM)[3] <- "pct.KO"
names(markers_mef2c_E85_KOvWT_IFTCM)[4] <- "pct.WT"
fwrite(markers_mef2c_E85_KOvWT_IFTCM, row.names = TRUE, file = "TableS3D_E85_IFT-CMs_KOvWT_DEGs.csv")
  
#Test for DEG between KO and WT cells in V CMs cluster
Idents(mef2c_v13_E85_v3_subset) <- "harmony_cell_type_subset"
markers_mef2c_E85_KOvWT_VCM <- FindMarkers(mef2c_v13_E85_v3_subset, assay = "SCT", ident.1 = "KO", ident.2 = "WT", group.by = "Genotype", subset.ident = "CMs-V")
markers_mef2c_E85_KOvWT_VCM <- subset(markers_mef2c_E85_KOvWT_VCM, p_val_adj < 0.05)
markers_mef2c_E85_KOvWT_VCM <- markers_mef2c_E85_KOvWT_VCM[order(markers_mef2c_E85_KOvWT_VCM$avg_log2FC, decreasing = TRUE), ]
names(markers_mef2c_E85_KOvWT_VCM)[3] <- "pct.KO"
names(markers_mef2c_E85_KOvWT_VCM)[4] <- "pct.WT"
fwrite(markers_mef2c_E85_KOvWT_VCM, row.names = TRUE, file = "TableS3E_E85_V-CMs_KOvWT_DEGs.csv")

#Test for DEG between KO and WT cells in OFT CMs cluster
Idents(mef2c_v13_E85_v3_subset) <- "harmony_cell_type_subset"
markers_mef2c_E85_KOvWT_OFTCM <- FindMarkers(mef2c_v13_E85_v3_subset, asay = "SCT", ident.1 = "KO", ident.2 = "WT", group.by = "Genotype", subset.ident = "CMs-OFT")
markers_mef2c_E85_KOvWT_OFTCM <- subset(markers_mef2c_E85_KOvWT_OFTCM, p_val_adj < 0.05)
markers_mef2c_E85_KOvWT_OFTCM <- markers_mef2c_E85_KOvWT_OFTCM[order(markers_mef2c_E85_KOvWT_OFTCM$avg_log2FC, decreasing = TRUE), ]
names(markers_mef2c_E85_KOvWT_OFTCM)[3] <- "pct.KO"
names(markers_mef2c_E85_KOvWT_OFTCM)[4] <- "pct.WT"
fwrite(markers_mef2c_E85_KOvWT_OFTCM, row.names = TRUE, file = "TableS3F_E85_OFT-CMs_KOvWT_DEGs.csv")


#E9

#Test for DEG between KO and WT cells in A CMs clusters
Idents(mef2c_v13_E9_subset) <- "cell_type_subset"
markers_mef2c_E9_KOvWT_ACM <- FindMarkers(mef2c_v13_E9_subset, ident.1 = "KO", ident.2 = "WT", group.by = "Genotype", subset.ident = c("CMs-A1", "CMs-A2"))
markers_mef2c_E9_KOvWT_ACM <- subset(markers_mef2c_E9_KOvWT_ACM, p_val_adj < 0.05)
markers_mef2c_E9_KOvWT_ACM <- markers_mef2c_E9_KOvWT_ACM[order(markers_mef2c_E9_KOvWT_ACM$avg_log2FC, decreasing = TRUE), ]
names(markers_mef2c_E9_KOvWT_ACM)[3] <- "pct.KO"
names(markers_mef2c_E9_KOvWT_ACM)[4] <- "pct.WT"
fwrite(markers_mef2c_E9_KOvWT_ACM, row.names = TRUE, file = "TableS3G_E9_A-CMs_KOvWT_DEGs.csv")

#Test for DEG between KO and WT cells in AVC CMs cluster
Idents(mef2c_v13_E9_subset) <- "cell_type_subset"
markers_mef2c_E9_KOvWT_AVCCM <- FindMarkers(mef2c_v13_E9_subset, ident.1 = "KO", ident.2 = "WT", group.by = "Genotype", subset.ident = "CMs-AVC")
markers_mef2c_E9_KOvWT_AVCCM <- subset(markers_mef2c_E9_KOvWT_AVCCM, p_val_adj < 0.05)
markers_mef2c_E9_KOvWT_AVCCM <- markers_mef2c_E9_KOvWT_AVCCM[order(markers_mef2c_E9_KOvWT_AVCCM$avg_log2FC, decreasing = TRUE), ]
names(markers_mef2c_E9_KOvWT_AVCCM)[3] <- "pct.KO"
names(markers_mef2c_E9_KOvWT_AVCCM)[4] <- "pct.WT"
fwrite(markers_mef2c_E9_KOvWT_AVCCM, row.names = TRUE, file = "TableS3H_E9_AVC-CMs_KOvWT_DEGs.csv")

#Test for DEG between KO and WT cells in V CMs cluster
Idents(mef2c_v13_E9_subset) <- "cell_type_subset"
markers_mef2c_E9_KOvWT_VCM <- FindMarkers(mef2c_v13_E9_subset, ident.1 = "KO", ident.2 = "WT", group.by = "Genotype", subset.ident = "CMs-V")
markers_mef2c_E9_KOvWT_VCM <- subset(markers_mef2c_E9_KOvWT_VCM, p_val_adj < 0.05)
markers_mef2c_E9_KOvWT_VCM <- markers_mef2c_E9_KOvWT_VCM[order(markers_mef2c_E9_KOvWT_VCM$avg_log2FC, decreasing = TRUE), ]
names(markers_mef2c_E9_KOvWT_VCM)[3] <- "pct.KO"
names(markers_mef2c_E9_KOvWT_VCM)[4] <- "pct.WT"
fwrite(markers_mef2c_E9_KOvWT_VCM, row.names = TRUE, file = "TableS3I_E9_V-CMs_KOvWT_DEGs.csv")

#Test for DEG between KO and WT cells in OFT CMs cluster
Idents(mef2c_v13_E9_subset) <- "cell_type_subset"
markers_mef2c_E9_KOvWT_OFTCM <- FindMarkers(mef2c_v13_E9_subset, ident.1 = "KO", ident.2 = "WT", group.by = "Genotype", subset.ident = "CMs-OFT")
markers_mef2c_E9_KOvWT_OFTCM <- subset(markers_mef2c_E9_KOvWT_OFTCM, p_val_adj < 0.05)
markers_mef2c_E9_KOvWT_OFTCM <- markers_mef2c_E9_KOvWT_OFTCM[order(markers_mef2c_E9_KOvWT_OFTCM$avg_log2FC, decreasing = TRUE), ]
names(markers_mef2c_E9_KOvWT_OFTCM)[3] <- "pct.KO"
names(markers_mef2c_E9_KOvWT_OFTCM)[4] <- "pct.WT"
fwrite(markers_mef2c_E9_KOvWT_OFTCM, row.names = TRUE, file = "TableS3J_E9_OFT-CMs_KOvWT_DEGs.csv")
