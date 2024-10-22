#This script is a starting point for analyzing scRNA-seq data with Seurat using
#scTransform for normalization and Harmony for batch correction. 
#Specifically, this script is for analyzing the 4 samples from the Mef2c Multiome
#E7.75 timepoint

#Script written by Jonathon M. Muncie-Vasic
#edited 10/22/2024


#load libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(magrittr)
library(XLConnect)
library(data.table)
library(harmony)
library(sctransform)


###---------Create Seurat Object----------------------------------------------

#read in 10X matrices and create Seurat objects

#Note: unlike pure scRNA-seq data, when importing Multiome data, all.data contains 
#a list with two entries: Gene Expression and Peaks so we need to refer to the 
#first entry of that list "Gene Expression" using the [[1]]
sample01.data <- Read10X(data.dir = "~/Desktop/Mef2c_v13_Seurat_scTrans_working/Inputs/sample01")
sample01 <- CreateSeuratObject(counts = sample01.data[[1]], project = "mef2c_v13_E775_sample01")
sample02.data <- Read10X(data.dir = "~/Desktop/Mef2c_v13_Seurat_scTrans_working/Inputs/sample02")
sample02 <- CreateSeuratObject(counts = sample02.data[[1]], project = "mef2c_v13_E775_sample02")
sample03.data <- Read10X(data.dir = "~/Desktop/Mef2c_v13_Seurat_scTrans_working/Inputs/sample03")
sample03 <- CreateSeuratObject(counts = sample03.data[[1]], project = "mef2c_v13_E775_sample03")
sample04.data <- Read10X(data.dir = "~/Desktop/Mef2c_v13_Seurat_scTrans_working/Inputs/sample04")
sample04 <- CreateSeuratObject(counts = sample04.data[[1]], project = "mef2c_v13_E775_sample04")

mef2c_v13_E775 <- merge(sample01, 
                        y = c(sample02, sample03, sample04), 
                        add.cell.ids = c("sample01", "sample02", "sample03", "sample04"), 
                        project = "mef2c_v13_E775")


###---------Add Meta Data Info----------------------------------------------

#Add sample number metadata column
class.vec1 <- gsub("_.*","", (colnames(mef2c_v13_E775@assays$RNA@data)))
names(class.vec1) <- colnames(mef2c_v13_E775@assays$RNA@data)
mef2c_v13_E775 <- AddMetaData(mef2c_v13_E775, class.vec1, "Sample_Number")

#Add sample name metadata column
class.vec1[class.vec1=="sample01"] <- "E7.75_WT1"
class.vec1[class.vec1=="sample02"] <- "E7.75_KO1"
class.vec1[class.vec1=="sample03"] <- "E7.75_KO2"
class.vec1[class.vec1=="sample04"] <- "E7.75_WT2"
mef2c_v13_E775 <- AddMetaData(mef2c_v13_E775, class.vec1, "Sample_Name")

#Add genotype metadata column
class.vec1[class.vec1 == "E7.75_WT1" | class.vec1 == "E7.75_WT2"] <- "WT"
class.vec1[class.vec1 == "E7.75_KO1" | class.vec1 == "E7.75_KO2"] <- "KO"
mef2c_v13_E775 <- AddMetaData(mef2c_v13_E775, class.vec1, "Genotype")

#Add percent mitochondria metadata column
mef2c_v13_E775 <- PercentageFeatureSet(mef2c_v13_E775, pattern = "^mt-", col.name = "percent.mt")

#Add perecent ribosomal reads to metadata column
#Initially, scTransform pipeline led to a cluster that was marked entirely by 
#ribosomal transcripts. Creating this metadata column so I can regress on them
mef2c_v13_E775 <- PercentageFeatureSet(mef2c_v13_E775, pattern = "^Rp[ls]", col.name = "percent.rib")

save(mef2c_v13_E775, file = "~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Objects/mef2c_v13_E775_scTrans.Robj")


###---------QC and Filtering----------------------------------------------

# Visualize QC metrics, grouped by sample, as a violin plot
# Note: the use of factor with set levels in the second line below will set the 
#order the groups are displayed any time that you use split.by Sample_Name
mef2c_v13_E775$Sample_Name <- factor(mef2c_v13_E775$Sample_Name, levels = c("E7.75_WT1", "E7.75_WT2", "E7.75_KO1", "E7.75_KO2"))
Idents(mef2c_v13_E775) <- "Sample_Name"
vlp1 <- VlnPlot(mef2c_v13_E775, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, split.by = "Sample_Name", pt.size = 0.1)
setwd("~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Outs/E775")
ggsave("QC_ViolinPlots_Mef2c_E775_prefilters.pdf", plot = vlp1, device = "pdf", width = 11, height = 5.5, dpi = 300)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
Idents(mef2c_v13_E775) <- "Sample_Name"
plot1 <- FeatureScatter(mef2c_v13_E775, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "Sample_Name", shuffle = T, pt.size = 0.5)
plot2 <- FeatureScatter(mef2c_v13_E775, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "Sample_Name", shuffle = T, pt.size = 0.5)
plot3 <- FeatureScatter(mef2c_v13_E775, feature1 = "percent.rib", feature2 = "percent.mt", group.by = "Sample_Name", shuffle = T, pt.size = 0.5)
fsc <- plot1+plot2+plot3
setwd("~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Outs/E775")
ggsave("QC_FeatureScatters_Mef2c_E775_prefilters.pdf", plot = fsc, device = "pdf", width = 15, height = 4, dpi = 300)

# Now filtering out GEMs with:
#  < 2000 features or > 8000 features
#  > 100000 counts
#  > 20% mt
#  > 20% ribosomal content
mef2c_v13_E775 <- subset(mef2c_v13_E775, subset = nFeature_RNA > 2000 & nFeature_RNA < 8000 & nCount_RNA < 100000 & percent.mt < 20 & percent.rib < 20)
#Removed ~7000 cells (28304 to 21237)

# Visualize QC metrics for filtered dataset, grouped by sample, as a violin plot
Idents(mef2c_v13_E775) <- "Sample_Name"
vlp2 <- VlnPlot(mef2c_v13_E775, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rib"), ncol = 4, split.by = "Sample_Name", pt.size = 0.1)
plot4 <- FeatureScatter(mef2c_v13_E775, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "Sample_Name", shuffle = T, pt.size = 0.5)
plot5 <- FeatureScatter(mef2c_v13_E775, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "Sample_Name", shuffle = T, pt.size = 0.5)
plot6 <- FeatureScatter(mef2c_v13_E775, feature1 = "percent.rib", feature2 = "percent.mt", group.by = "Sample_Name", shuffle = T, pt.size = 0.5)
fsc2 <- plot4+plot5+plot6
setwd("~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Outs/E775")
ggsave("QC_ViolinPlots_Mef2c_E775_postfilters.pdf", plot = vlp2, device = "pdf", width = 15, height = 5.5, dpi = 300)
ggsave("QC_FeatureScatters_Mef2c_E775_postfilters.pdf", plot = fsc2, device = "pdf", width = 15, height = 4, dpi = 300)


###------------------Run scTransform and Dim Reduction---------------------------
#Replaces/combines normalization, scaling, and finding variable features

# Note: glmGamPoi speeds up scTransform - if didn't load properly, just remove that argument
# Note: as in standard workflow, only the variable genes are saved to the scale.data
# matrices in the output assay 

mef2c_v13_E775 <- SCTransform(mef2c_v13_E775, 
                              method = "glmGamPoi", 
                              variable.features.n = 5000,
                              vars.to.regress = c("percent.mt", "percent.rib"), 
                              verbose = TRUE,
                              vst.flavor = "v2")

mef2c_v13_E775 <- RunPCA(mef2c_v13_E775, assay = "SCT")
elb <- ElbowPlot(mef2c_v13_E775, ndims = 40)
elb


###------Clustering and UMAPs with Harmony Batch Correction--------------
mef2c_v13_E775 <- RunHarmony(mef2c_v13_E775, group.by.vars = "Sample_Name", assay.use = "SCT", dims.use = 1:30)
mef2c_v13_E775 <- FindNeighbors(mef2c_v13_E775, reduction = "harmony", dims = 1:30)
mef2c_v13_E775 <- FindClusters(mef2c_v13_E775, resolution = 0.5, random.seed = 0)
mef2c_v13_E775 <- RunUMAP(mef2c_v13_E775, reduction = "harmony", assay = "SCT", dims = 1:30, seed.use = 42)

#Visualize UMAPs labeled by sample
UMAP_E775_Harmony_bysample <- DimPlot(object = mef2c_v13_E775, 
                                      reduction = "umap", 
                                      group.by = "Sample_Name", 
                                      cols = c("lightskyblue","turquoise1","blue","darkblue"), 
                                      shuffle = T, seed = 1, pt.size = 0.1
                                      ) + NoLegend()
UMAP_E775_Harmony_bysample_legend <- DimPlot(object = mef2c_v13_E775, 
                                             reduction = "umap", 
                                             group.by = "Sample_Name", 
                                             cols = c("lightskyblue","turquoise1","blue","darkblue"), 
                                             shuffle = T, seed = 1, pt.size = 0.1
                                             ) 
setwd("~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Outs/E775")
ggsave('UMAP_Mef2c_E775_Harmony_bysample.pdf', plot = UMAP_E775_Harmony_bysample + labs(title = 'Mef2c E7.75, Harmony corrected, 21237 cells'), device = 'pdf', width = 5, height = 5, dpi = 300)
ggsave('UMAP_Mef2c_E775_Harmony_bysample_legend.pdf', plot = UMAP_E775_Harmony_bysample_legend + labs(title = 'Mef2c E7.75, Harmony corrected, 21237 cells'), device = 'pdf', width = 5, height = 5, dpi = 300)

#Visualize UMAPs labeled by cluster
Idents(mef2c_v13_E775) <- "SCT_snn_res.0.5"
UMAP_E775_Harmony_res0.5 <- DimPlot(object = mef2c_v13_E775, 
                                    reduction = "umap", 
                                    label = T, raster = F, 
                                    label.size = 5, pt.size = 0.1
                                    ) + NoLegend()
setwd("~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Outs/E775")
ggsave('UMAP_Mef2c_E775_Harmony_res0.5.pdf', plot = UMAP_E775_Harmony_res0.5 + labs(title = 'Mef2c E7.75, Harmony corrected, 21237 cells, res 0.5'), device = 'pdf', width = 5, height = 5, dpi = 300)


###----------Identify Clusters---------------------------------------
#Make feature plots for genes used to identify clusters 
#Create list of of marker genes
features1 = c(
  "Pdgfra", "Lef1", "Cdh2",             #Posterior Mesoderm
  "Sox2", "Pax6", "Ror1", "Ror2",       #Ectoderm
  "T", "Wnt3a", "Cdx1", "Ncam1",        #Ectoderm & PS
  "Ttr","Cubn","Gata4",                 #Visceral Endoderm
  "Hand1", "Bmp4", "Ahnak",             #Extraembryonic Mesoderm
  "Foxp2", "Foxc2", "Tbx1",             #Anterior Mesoderm
  "Tfap2a", "Tfap2c",                   #Extraembryonic Ectoderm
  "Tbx4",                               #Allantois
  "Cdh5", "Kdr", "Flt1",                #Endothelium
  "Nfatc1",                             #Endocardium
  "Sox17", "Trh", "Epcam",              #Definitive Endoderm
  "Nkx2-5", "Smarcd3", "Tbx5",          #Cardiac Progenitors/Myocytes
  "Mef2c","Tnnt2",
  "Wnt6", "Flt1",                       #Extraembryonic Endoderm 
  "Nog", "Shh",                         #Axial Mesoderm
  "Gjb3", "Wnt7b",                      #Trophoblast/Placental
  "Hba-a1", "Hba-a2", "Gata1",          #Blood
  "Tfap2c", "Prdm1", "Ifitm3", "Dnd1"   #PGCs (Prdm1 = Blimp1, Ifitm3 = Fragilis)
)

#Use for loop to save plots
setwd("~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Outs/E775/Feature_Plots")
for (gene in features1) {
  gene_featureplot <- plot(
    FeaturePlot(mef2c_v13_E775, features = gene)
  )
  ggsave(file = paste0("FP_Mef2c_E775_", gene, ".pdf"), plot = gene_featureplot, device = "pdf", width = 5, height = 5, dpi = 300)
}

#Run find markers for additional unbiased cluster markers
#Should return markers for 17 clusters, corresponding to res 0.5
Idents(mef2c_v13_E775) <- "SCT_snn_res.0.5"
mef2c_E775.markers <- FindAllMarkers(mef2c_v13_E775, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
mef2c_E775.markers <- subset(mef2c_E775.markers, p_val_adj < 0.05)
mef2c_E775.markers <- mef2c_E775.markers[order(mef2c_E775.markers$avg_log2FC, decreasing = TRUE), ]
#Save to spreadsheet
setwd("~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Outs/E775")
fwrite(mef2c_E775.markers, row.names = TRUE, file = "Mef2c_v13_E775_markers_res0.5.csv")

#Determined cell types for each cluster:
#C0 = Posterior Mesoderm
#C1 = Ectoderm 1 
#C2 = Visceral Endoderm
#C3 = Ectoderm 2 and Primitive Streak
#C4 = Ectoderm 3
#C5 = Blood
#C6 = Allantois
#C7 = Anterior Mesoderm
#C8 = Extraembryonic Mesoderm
#C9 = Endothelium (includes endocardium)
#C10 = Extraembryonic Ectoderm
#C11 = Ectoderm 4
#C12 = Definitive Endoderm 
#C13 = CPs and CMs
#C14 = Axial Mesoderm
#C15 = Trophoblast/Placental
#C16 = C16 

#Now label clusters by cell type
Idents(mef2c_v13_E775) <- "SCT_snn_res.0.5"
new.cluster.ids <- c("PostM", "Ect1", "VE", "Ect2/PS",
                     "Ect3", "Blood", "Allantois", "AntM", "ExM",
                     "Endo", "ExEct", "Ect4", "DE",
                     "CPs/CMs", "AxM", "T/P", "C16"
                     )
names(new.cluster.ids) <- levels(mef2c_v13_E775$SCT_snn_res.0.5)
mef2c_v13_E775 <- RenameIdents(mef2c_v13_E775, new.cluster.ids)

#Add "cell_type" metadata column 
mef2c_v13_E775 <- AddMetaData(mef2c_v13_E775, mef2c_v13_E775@active.ident, "cell_type")

#Plot and save labeled UMAP
Idents(mef2c_v13_E775) <- "cell_type"
UMAP_E775_Harmony_res0.5_labeled <- DimPlot(object = mef2c_v13_E775, 
                                            reduction = "umap", 
                                            label = T, repel = F, 
                                            label.size = 3.5, pt.size = 0.1
                                            ) + NoLegend()
setwd("~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Outs/E775")
ggsave('UMAP_Mef2c_E775_Harmony_res0.5_labeled.pdf', plot = UMAP_E775_Harmony_res0.5_labeled + labs(title = 'Mef2c E7.75, Harmony corrected, 21237 cells, res 0.5'), device = 'pdf', width = 5, height = 5, dpi = 300)

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
features2 = c(
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
DP_E775 <- DotPlot(mef2c_v13_E775, features = features2) + RotatedAxis()
setwd("~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Outs/E775")
ggsave('DotPlot_mef2c_E775.pdf', plot = DP_E775, device = 'pdf', width = 15, height = 5, dpi = 300)


###-----------------Subset----------------------------------------------
#Subset clusters of interest, namely CPs/CMs, Anterior Mesoderm, Extraembryonic Mesoderm
Idents(mef2c_v13_E775) <- "cell_type"
mef2c_v13_E775_subset <- subset(mef2c_v13_E775, idents = c("CPs/CMs", "AntM", "ExM"))
mef2c_v13_E775_subset
#The subset contains 3356 cells

save(mef2c_v13_E775, file = "~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Objects/mef2c_v13_E775_scTrans.Robj")
save(mef2c_v13_E775_subset, file = "~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Objects/mef2c_v13_E775_subset_scTrans.Robj")


###-------------Subset - Run scTransform and Dim Reduction---------------------------
mef2c_v13_E775_subset <- SCTransform(mef2c_v13_E775_subset, 
                                     method = "glmGamPoi", 
                                     variable.features.n = 5000,
                                     vars.to.regress = c("percent.mt", "percent.rib"), 
                                     verbose = TRUE,
                                     vst.flavor = "v2")

mef2c_v13_E775_subset <- RunPCA(mef2c_v13_E775_subset, assay = "SCT")
elb <- ElbowPlot(mef2c_v13_E775_subset, ndims = 40)
elb


###------Subset - Clustering and UMAPs with Harmony Batch Correction--------------
mef2c_v13_E775_subset <- RunHarmony(mef2c_v13_E775_subset, group.by.vars = "Sample_Name", assay.use = "SCT", dims.use = 1:30)
mef2c_v13_E775_subset <- FindNeighbors(mef2c_v13_E775_subset, reduction = "harmony", dims = 1:30)
mef2c_v13_E775_subset <- FindClusters(mef2c_v13_E775_subset, resolution = 2, random.seed = 0)
mef2c_v13_E775_subset <- RunUMAP(mef2c_v13_E775_subset, reduction = "harmony", assay = "SCT", dims = 1:30, seed.use = 42)

#Visualize UMAP with sample labels
UMAP_E775_subset_Harmony_bysample <- DimPlot(object = mef2c_v13_E775_subset, 
                                             reduction = "umap", 
                                             group.by = "Sample_Name", 
                                             cols = c("lightskyblue","turquoise1","blue","darkblue"), 
                                             shuffle = T, seed = 1, pt.size = 0.5
                                             ) + NoLegend()
setwd("~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Outs/E775")
ggsave('UMAP_Mef2c_E775_subset_Harmony_bysample.pdf', plot = UMAP_E775_subset_Harmony_bysample + labs(title = 'Mef2c E7.75 Subset, Harmony corrected, 3356 cells'), device = 'pdf', width = 5, height = 5, dpi = 300)

#Visualize UMAP with cluster labels
Idents(mef2c_v13_E775_subset) <- "SCT_snn_res.2"
UMAP_E775_subset_Harmony_res2 <- DimPlot(object = mef2c_v13_E775_subset, 
                                         reduction = "umap", 
                                         label = T, raster = F, 
                                         label.size = 5, pt.size = 0.5
                                         ) + NoLegend()
setwd("~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Outs/E775")
ggsave('UMAP_Mef2c_E775_subset_Harmony_res2.pdf', plot = UMAP_E775_subset_Harmony_res2 + labs(title = 'Mef2c E7.75 Subset, Harmony corrected, 3356 cells, res 2.0'), device = 'pdf', width = 5, height = 5, dpi = 300)


###----------Subset - Identify Clusters---------------------------------------
#Make feature plots for genes used to identify clusters 
#Create list of of marker genes
features3 = c(
  "Hand1", "Bmp4",                      #ExMeso
  "Sobp", "Otx2", "Cyp26c1",            #Cranial Meso
  "Isl1", "Fgf8", "Fgf10", "Mef2c",     #SHF
  "Gata4", "Gata6",                     #Lateral Mesoderm
  "Fst", "Meox1",                       #Somitic Meso
  "Nkx1-2", "Sox2", "T",                #NMPs
  "Tbx1",                               #Paraxial Meso
  "Tbx5", "Nkx2-5", "Mab21l2",          #JCF
  "Runx1", "Itga4", "Cd44",             #HSCs
  "Tnnt2", "Ttn",                       #CMs
  "Slc7a8", "Slc2a2", "Cubn", "Amn"     #KPs
)

#Use for loop to save plots
setwd("~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Outs/E775/Feature_Plots_Subset")
for (gene in features3) {
  gene_featureplot <- plot(
    FeaturePlot(mef2c_v13_E775_subset, features = gene)
  )
  ggsave(file = paste0("FP_Mef2c_E775_subset_", gene, ".pdf"), plot = gene_featureplot, device = "pdf", width = 5, height = 5, dpi = 300)
}

#Run find markers for additional unbiased cluster markers
#Should return markers for 21 clusters, corresponding to res 2.0
Idents(mef2c_v13_E775_subset) <- "SCT_snn_res.2"
mef2c_E775_subset.markers <- FindAllMarkers(mef2c_v13_E775_subset, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
mef2c_E775_subset.markers <- subset(mef2c_E775_subset.markers, p_val_adj < 0.05)
mef2c_E775_subset.markers <- mef2c_E775_subset.markers[order(mef2c_E775_subset.markers$avg_log2FC, decreasing = TRUE), ]
#Save to spreadsheet
setwd("~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Outs/E775")
fwrite(mef2c_E775_subset.markers, row.names = TRUE, file = "Mef2c_v13_E775_subset_markers_res2.csv")

#Determined cell types for each cluster:
#C0 = Cranial Mesoderm 1
#C1 = Somitic Mesoderm 1
#C2 = Extraembryonic Mesoderm 1
#C3 = NMPs 
#C4 = Extraembryonic Mesoderm 2
#C5 = Paraxial Mesoderm
#C6 = Cranial Mesoderm 2
#C7 = Somitic Mesoderm 2
#C8 = Extraembryonic Mesoderm 3
#C9 = Extraembryonic Mesoderm 4
#C10 = Lateral Plate Mesoderm  
#C11 = SHF1
#C12 = Extraembryonic Mesoderm 5
#C13 = Extraembryonic Mesoderm 6
#C14 = JCF
#C15 = Somitic Mesoderm 3
#C16 = CMs/FHF
#C17 = HSCs
#C18 = SHF2
#C19 = Cranial Mesoderm 3
#C20 = KPs

#Now label clusters by cell type
Idents(mef2c_v13_E775_subset) <- "SCT_snn_res.2"
new.cluster.ids <- c("CrM1", "SoM1", "ExM1", "NMPs",
                     "ExM2", "PrxM", "CrM2", "SoM2",
                     "ExM3", "ExM4", "LPM", "SHF1",
                     "ExM5", "ExM6", "JCF", "SoM3",
                     "CMs/FHF", "HSCs", "SHF2", "CrM3", "KPs"
                     )
names(new.cluster.ids) <- levels(mef2c_v13_E775_subset$SCT_snn_res.2)
mef2c_v13_E775_subset <- RenameIdents(mef2c_v13_E775_subset, new.cluster.ids)

#Add "cell_type_subset" metadata column 
mef2c_v13_E775_subset <- AddMetaData(mef2c_v13_E775_subset, mef2c_v13_E775_subset@active.ident, "cell_type_subset")

#Plot and save labeled UMAP
Idents(mef2c_v13_E775_subset) <- "cell_type_subset"
UMAP_E775_subset_Harmony_res2_labeled <- DimPlot(object = mef2c_v13_E775_subset, 
                                                 reduction = "umap", 
                                                 label = T, repel = T, 
                                                 label.size = 4, pt.size = 0.5
                                                 ) + NoLegend()
setwd("~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Outs/E775")
ggsave('UMAP_Mef2c_E775_subset_Harmony_res2_labeled.pdf', plot = UMAP_E775_subset_Harmony_res2_labeled + labs(title = 'Mef2c E7.75 Subset, Harmony corrected, 3356 cells, res2'), device = 'pdf', width = 5, height = 5, dpi = 300)

#Plot and save UMAP labeled by full dataset labels to double-check that labels make sense
Idents(mef2c_v13_E775_subset) <- "cell_type"
UMAP_E775_subset_Harmony_labeled_full <- DimPlot(object = mef2c_v13_E775_subset, 
                                                 reduction = "umap", 
                                                 label = T, repel = F, 
                                                 label.size = 4, pt.size = 0.5
                                                 ) + NoLegend()
setwd("~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Outs/E775")
ggsave('UMAP_Mef2c_E775_subset_Harmony_labeled_fulldataset.pdf', 
       plot = UMAP_E775_subset_Harmony_labeled_full + 
         labs(title = 'Mef2c E7.75 Subset, Harmony corrected, 3356 cells') +
         theme(plot.title = element_text(size=10)),
       device = 'pdf', width = 5, height = 5, dpi = 300
)

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
features4 = c(
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
DP_E775_subset <- DotPlot(mef2c_v13_E775_subset, features = features4) + RotatedAxis()
setwd("~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Outs/E775")
ggsave('DotPlot_mef2c_E775_subset.pdf', plot = DP_E775_subset, device = 'pdf', width = 17, height = 7, dpi = 300)


###----------Subset - DEG Testing ------------------------------------------

#Test for DEG between KO and WT cells in CMs/FHF cluster
Idents(mef2c_v13_E775_subset) <- "cell_type_subset"
markers_mef2c_E775_KOvWT_CMFHF <- FindMarkers(mef2c_v13_E775_subset, ident.1 = "KO", ident.2 = "WT", group.by = "Genotype", subset.ident = "CMs/FHF")
markers_mef2c_E775_KOvWT_CMFHF <- subset(markers_mef2c_E775_KOvWT_CMFHF, p_val_adj < 0.05)
markers_mef2c_E775_KOvWT_CMFHF <- markers_mef2c_E775_KOvWT_CMFHF[order(markers_mef2c_E775_KOvWT_CMFHF$avg_log2FC, decreasing = TRUE), ]
names(markers_mef2c_E775_KOvWT_CMFHF)[3] <- "pct.KO"
names(markers_mef2c_E775_KOvWT_CMFHF)[4] <- "pct.WT"
setwd("~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Outs/E775")
fwrite(markers_mef2c_E775_KOvWT_CMFHF, row.names = TRUE, file = "Mef2c_v13_E775_subset_markers_KOvWT_CMFHF.csv")

#Test for DEG between KO and WT cells in SHF clusters
Idents(mef2c_v13_E775_subset) <- "cell_type_subset"
markers_mef2c_E775_KOvWT_SHF <- FindMarkers(mef2c_v13_E775_subset, ident.1 = "KO", ident.2 = "WT", group.by = "Genotype", subset.ident = c("SHF1", "SHF2"))
markers_mef2c_E775_KOvWT_SHF <- subset(markers_mef2c_E775_KOvWT_SHF, p_val_adj < 0.05)
markers_mef2c_E775_KOvWT_SHF <- markers_mef2c_E775_KOvWT_SHF[order(markers_mef2c_E775_KOvWT_SHF$avg_log2FC, decreasing = TRUE), ]
names(markers_mef2c_E775_KOvWT_SHF)[3] <- "pct.KO"
names(markers_mef2c_E775_KOvWT_SHF)[4] <- "pct.WT"
setwd("~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Outs/E775")
fwrite(markers_mef2c_E775_KOvWT_SHF, row.names = TRUE, file = "Mef2c_v13_E775_subset_markers_KOvWT_SHF.csv")

#Test for DEG between KO and WT cells in JCF cluster
Idents(mef2c_v13_E775_subset) <- "cell_type_subset"
markers_mef2c_E775_KOvWT_JCF <- FindMarkers(mef2c_v13_E775_subset, ident.1 = "KO", ident.2 = "WT", group.by = "Genotype", subset.ident = "JCF")
markers_mef2c_E775_KOvWT_JCF <- subset(markers_mef2c_E775_KOvWT_JCF, p_val_adj < 0.05)
markers_mef2c_E775_KOvWT_JCF <- markers_mef2c_E775_KOvWT_JCF[order(markers_mef2c_E775_KOvWT_JCF$avg_log2FC, decreasing = TRUE), ]
names(markers_mef2c_E775_KOvWT_JCF)[3] <- "pct.KO"
names(markers_mef2c_E775_KOvWT_JCF)[4] <- "pct.WT"
setwd("~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Outs/E775")
fwrite(markers_mef2c_E775_KOvWT_JCF, row.names = TRUE, file = "Mef2c_v13_E775_subset_markers_KOvWT_JCF.csv")


###----------Save Seurat Objects---------------------------------------
save(mef2c_v13_E775, file = "~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Objects/mef2c_v13_E775_scTrans.Robj")
save(mef2c_v13_E775_subset, file = "~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Objects/mef2c_v13_E775_subset_scTrans.Robj")

