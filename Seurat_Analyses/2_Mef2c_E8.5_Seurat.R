#This script is a starting point for analyzing scRNA-seq data with Seurat using
#scTransform for normalization and Harmony for batch correction. 
#Specifically, this script is for analyzing the 4 samples from the Mef2c Multiome
#E8.5 timepoint

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
sample21.data <- Read10X(data.dir = "~/Desktop/Mef2c_v13_Seurat_scTrans_working/Inputs/sample21")
sample21 <- CreateSeuratObject(counts = sample21.data[[1]], project = "mef2c_v13_E85_sample21")
sample22.data <- Read10X(data.dir = "~/Desktop/Mef2c_v13_Seurat_scTrans_working/Inputs/sample22")
sample22 <- CreateSeuratObject(counts = sample22.data[[1]], project = "mef2c_v13_E85_sample22")
sample23.data <- Read10X(data.dir = "~/Desktop/Mef2c_v13_Seurat_scTrans_working/Inputs/sample23")
sample23 <- CreateSeuratObject(counts = sample23.data[[1]], project = "mef2c_v13_E85_sample23")
sample24.data <- Read10X(data.dir = "~/Desktop/Mef2c_v13_Seurat_scTrans_working/Inputs/sample24")
sample24 <- CreateSeuratObject(counts = sample24.data[[1]], project = "mef2c_v13_E85_sample24")

mef2c_v13_E85_v3 <- merge(sample21, 
                          y = c(sample22, sample23, sample24), 
                          add.cell.ids = c("sample21", "sample22", "sample23", "sample24"), 
                          project = "mef2c_v13_E85_v3")


###---------Add Meta Data Info----------------------------------------------

#Add sample number metadata column
class.vec1 <- gsub("_.*","", (colnames(mef2c_v13_E85_v3@assays$RNA@data)))
names(class.vec1) <- colnames(mef2c_v13_E85_v3@assays$RNA@data)
mef2c_v13_E85_v3 <- AddMetaData(mef2c_v13_E85_v3, class.vec1, "Sample_Number")

#Add sample name metadata column
class.vec1[class.vec1=="sample21"] <- "E8.5_KO1"
class.vec1[class.vec1=="sample22"] <- "E8.5_WT1"
class.vec1[class.vec1=="sample23"] <- "E8.5_WT2"
class.vec1[class.vec1=="sample24"] <- "E8.5_KO2"
mef2c_v13_E85_v3 <- AddMetaData(mef2c_v13_E85_v3, class.vec1, "Sample_Name")

#Add genotype metadata column
class.vec1[class.vec1 == "E8.5_WT1" | class.vec1 == "E8.5_WT2"] <- "WT"
class.vec1[class.vec1 == "E8.5_KO1" | class.vec1 == "E8.5_KO2"] <- "KO"
mef2c_v13_E85_v3 <- AddMetaData(mef2c_v13_E85_v3, class.vec1, "Genotype")

#Add percent mitochondria metadata column
mef2c_v13_E85_v3 <- PercentageFeatureSet(mef2c_v13_E85_v3, pattern = "^mt-", col.name = "percent.mt")
#Add percent ribosomal reads to metadata column
#Initially, scTransform pipeline led to a cluster that was marked entirely by 
#ribosomal transcripts. Creating this metadata column so I can regress on them
mef2c_v13_E85_v3 <- PercentageFeatureSet(mef2c_v13_E85_v3, pattern = "^Rp[ls]", col.name = "percent.rib")

save(mef2c_v13_E85_v3, file = "~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Objects/mef2c_v13_E85_scTrans_v3_harmony.Robj")


###---------QC and Filtering----------------------------------------------

# Visualize QC metrics, grouped by sample, as a violin plot
# Note: the use of factor with set levels in the second line below will set the 
#order the groups are displayed any time that you use split.by Sample_Name
mef2c_v13_E85_v3$Sample_Name <- factor(mef2c_v13_E85_v3$Sample_Name, levels = c("E8.5_WT1", "E8.5_WT2", "E8.5_KO1", "E8.5_KO2"))
Idents(mef2c_v13_E85_v3) <- "Sample_Name"
vlp1 <- VlnPlot(mef2c_v13_E85_v3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, split.by = "Sample_Name", pt.size = 0.1)
setwd("~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Outs/E85_v3")
ggsave("QC_ViolinPlots_Mef2c_E85_prefilters.pdf", plot = vlp1, device = "pdf", width = 11, height = 5.5, dpi = 300)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
Idents(mef2c_v13_E85_v3) <- "Sample_Name"
plot1 <- FeatureScatter(mef2c_v13_E85_v3, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "Sample_Name", shuffle = T, pt.size = 0.5)
plot2 <- FeatureScatter(mef2c_v13_E85_v3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "Sample_Name", shuffle = T, pt.size = 0.5)
plot3 <- FeatureScatter(mef2c_v13_E85_v3, feature1 = "percent.rib", feature2 = "percent.mt", group.by = "Sample_Name", shuffle = T, pt.size = 0.5)
fsc <- plot1+plot2+plot3
setwd("~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Outs/E85_v3")
ggsave("QC_FeatureScatters_Mef2c_E85_prefilters.pdf", plot = fsc, device = "pdf", width = 15, height = 4, dpi = 300)

# Now filtering out GEMs with:
#  < 1750 features or >5500 features **Pretty aggressive/strict, intended to remove potential doublets
#  > 40000 counts, 
#  > 15% mt
#  > 20% rib
mef2c_v13_E85_v3
mef2c_v13_E85_v3 <- subset(mef2c_v13_E85_v3, subset = nFeature_RNA > 1750 & nFeature_RNA < 5500 & nCount_RNA < 40000 & percent.mt < 15 & percent.rib < 20)
mef2c_v13_E85_v3
#Removed ~6000 cells (34,354 to 28,373)

# Visualize QC metrics for filtered dataset, grouped by sample
Idents(mef2c_v13_E85_v3) <- "Sample_Name"
vlp2 <- VlnPlot(mef2c_v13_E85_v3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rib"), ncol = 4, split.by = "Sample_Name", pt.size = 0.1)
plot4 <- FeatureScatter(mef2c_v13_E85_v3, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "Sample_Name", shuffle = T, pt.size = 0.5)
plot5 <- FeatureScatter(mef2c_v13_E85_v3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "Sample_Name", shuffle = T, pt.size = 0.5)
plot6 <- FeatureScatter(mef2c_v13_E85_v3, feature1 = "percent.rib", feature2 = "percent.mt", group.by = "Sample_Name", shuffle = T, pt.size = 0.5)
fsc2 <- plot4+plot5+plot6
setwd("~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Outs/E85_v3")
ggsave("QC_ViolinPlots_Mef2c_E85_postfilters.pdf", plot = vlp2, device = "pdf", width = 15, height = 5.5, dpi = 300)
ggsave("QC_FeatureScatters_Mef2c_E85_postfilters.pdf", plot = fsc2, device = "pdf", width = 15, height = 4, dpi = 300)


###------------------Run scTransform and Dim Reduction---------------------------
#Replaces/combines normalization, scaling, and finding variable features

# Note: glmGamPoi speeds up scTransform - if didn't load properly, just remove that argument
# Note: as in standard workflow, only the variable genes are saved to the scale.data
# matrices in the output assay 

mef2c_v13_E85_v3 <- SCTransform(mef2c_v13_E85_v3, 
                             method = "glmGamPoi", 
                             variable.features.n = 5000,
                             vars.to.regress = c("percent.mt", "percent.rib"), 
                             verbose = TRUE,
                             vst.flavor = "v2")

mef2c_v13_E85_v3 <- RunPCA(mef2c_v13_E85_v3, assay = "SCT")
elb <- ElbowPlot(mef2c_v13_E85_v3, ndims = 40)
elb


###------Clustering and UMAPs with Harmony Batch Correction--------------
mef2c_v13_E85_v3 <- RunHarmony(mef2c_v13_E85_v3, group.by.vars = "Sample_Name", assay.use = "SCT", dims.use = 1:30)
mef2c_v13_E85_v3 <- FindNeighbors(mef2c_v13_E85_v3, reduction = "harmony", dims = 1:30)
mef2c_v13_E85_v3 <- FindClusters(mef2c_v13_E85_v3, resolution = 0.4, random.seed = 0)
mef2c_v13_E85_v3 <- RunUMAP(mef2c_v13_E85_v3, reduction = "harmony", assay = "SCT", dims = 1:30, seed.use = 42)

#Visualize UMAP with sample labels
UMAP_E85_Harmony_bysample <- DimPlot(object = mef2c_v13_E85_v3, 
                                     reduction = "umap", 
                                     group.by = "Sample_Name", 
                                     cols = c("#F3E5F5", "#D1C4E9", "#9C27B0", "#4527A0"), 
                                     shuffle = T, seed = 1, pt.size = 0.1
                                     ) + NoLegend()
UMAP_E85_Harmony_bysample_legend <- DimPlot(object = mef2c_v13_E85_v3, 
                                            reduction = "umap", 
                                            group.by = "Sample_Name", 
                                            cols = c("#F3E5F5", "#D1C4E9", "#9C27B0", "#4527A0"), 
                                            shuffle = T, seed = 1, pt.size = 0.1
                                            ) 
setwd("~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Outs/E85_v3")
ggsave('UMAP_Mef2c_E85_Harmony_bysample.pdf', plot = UMAP_E85_Harmony_bysample + labs(title = 'Mef2c E8.5 Harmony, 28373 cells'), device = 'pdf', width = 5, height = 5, dpi = 300)
ggsave('UMAP_Mef2c_E85_Harmony_bysample_legend.pdf', plot = UMAP_E85_Harmony_bysample_legend + labs(title = 'Mef2c E8.5 Harmony, 28373 cells'), device = 'pdf', width = 5, height = 5, dpi = 300)

#Visualize UMAP with cluster labels
Idents(mef2c_v13_E85_v3) <- "SCT_snn_res.0.4"
UMAP_E85_Harmony_res0.4 <- DimPlot(object = mef2c_v13_E85_v3, 
                                   reduction = "umap", 
                                   label = T, repel = F, raster = F, 
                                   label.size = 5, pt.size = 0.1
                                   ) + NoLegend()
setwd("~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Outs/E85_v3")
ggsave('UMAP_Mef2c_E85_Harmony_res0.4.pdf', plot = UMAP_E85_Harmony_res0.4 + labs(title = 'Mef2c E8.5 Harmony, 28373 cells, res 0.4'), device = 'pdf', width = 5, height = 5, dpi = 300)


###----------Identify Clusters---------------------------------------
#Make feature plots for genes used to identify clusters 
#Create list of of marker genes
features1 = c(
  "Pdgfra",                                       #Mesoderm
  "Gata6", "Alx1", "Dlk1",                        #Lateral Mesoderm
  "Foxf1", "Foxp2",                               #Splanchnic Mesoderm
  "Sox2", "Fgf14", "Nrcam",                       #Neural Tube
  "Meox1", "Fst", "Cdh11",                        #Somitic Mesoderm
  "Sox17","Epcam", "Trh", "Cdh1", "Cldn6",        #Definitive Endoderm
  "Isl1", "Tbx1", "Fgf8", "Fgf10", "Hand2",       #SHF
  "Nkx2-5", "Tbx5", "Mef2c", "Mab21l2",           #Cardiac markers
  "Ttn", "Tnnt2", "Myl7",                         #CMs
  "Tbx18", "Wt1",                                 #Proepicardium
  "Tfap2b", "Sox10", "Ets1",                      #Neural Crest
  "Cdh5", "Kdr", "Flt1",                          #Endothelium
  "Nfatc1",                                       #Endocardium
  "Hand1", "Bmp4", "Ahnak",                       #Extraembryonic Mesoderm
  "Hba-a1", "Hba-a2", "Hbb-y",                    #Blood
  "Ttr", "Afp", "Apoa1",                          #Yolk Sac
  "T", "Shh",                                     #Notochord
  "Postn", "Col1a2", "Col3a1"                     #Mesenchyme
)

#Use for loop to save plots
setwd("~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Outs/E85_v3/Feature_Plots_harmony")
for (gene in features1) {
  gene_featureplot <- plot(
    FeaturePlot(mef2c_v13_E85_v3, features = gene)
  )
  ggsave(file = paste0("FP_Mef2c_E85_", gene, ".pdf"), plot = gene_featureplot, device = "pdf", width = 5, height = 5, dpi = 300)
}

#Run find markers for additional unbiased cluster markers
#Should return markers for 18 clusters, corresponding to res 0.4
Idents(mef2c_v13_E85_v3) <- "SCT_snn_res.0.4"
mef2c_E85.markers <- FindAllMarkers(mef2c_v13_E85_v3, assay = "SCT", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.6)
mef2c_E85.markers <- subset(mef2c_E85.markers, p_val_adj < 0.05)
mef2c_E85.markers <- mef2c_E85.markers[order(mef2c_E85.markers$avg_log2FC, decreasing = TRUE), ]
#Save to spreadsheet
setwd("~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Outs/E85_v3")
fwrite(mef2c_E85.markers, row.names = TRUE, file = "mef2c_v13_E85_v3_markers_res0.4.csv")

#Determined cell type labels for each cluster: 
#C0 = Neural Tube 1
#C1 = Second Heart Field 1/Mixed Mesoderm 
#C2 = Endoderm 1
#C3 = Endoderm 2
#C4 = Somitic Mesoderm
#C5 = CMs1
#C6 = Second Heart Field 2/Pharyngel Mesoderm
#C7 = Neural Tube 2
#C8 = Neural Tube 3
#C9 = Second Heart Field 3/Lateral Plate Mesoderm
#C10 = CMs2/CPs
#C11 = Neural Crest
#C12 = Neural Tube 4
#C13 = Neural Tube 5
#C14 = Endothelium
#C15 = Endoderm 3
#C16 = Notochord
#C17 = Yolk Sac

#Now label clusters by cell type
Idents(mef2c_v13_E85_v3) <- "SCT_snn_res.0.4"
new.cluster.ids <- c("NT1", "SHF1/MixM", "En1", "En2",
                     "SoM", "CMs1", "SHF2/PhM", "NT2", 
                     "NT3", "SHF3/LPM", "CMs2/CPs", "NC",
                     "NT4", "NT5", "Endo", "En3",
                     "Noto", "YS")
names(new.cluster.ids) <- levels(mef2c_v13_E85_v3$SCT_snn_res.0.4)
mef2c_v13_E85_v3 <- RenameIdents(mef2c_v13_E85_v3, new.cluster.ids)

#Add "harmony_cell_type" metadata column 
mef2c_v13_E85_v3 <- AddMetaData(mef2c_v13_E85_v3, mef2c_v13_E85_v3@active.ident, "harmony_cell_type")

#Plot and save labeled UMAP
Idents(mef2c_v13_E85_v3) <- "harmony_cell_type"
UMAP_E85_harm_res0.4_labeled <- DimPlot(object = mef2c_v13_E85_v3, 
                                        reduction = "umap", 
                                        label = T, repel = T, 
                                        label.size = 4, pt.size = 0.1
                                        ) + NoLegend()
setwd("~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Outs/E85_v3")
ggsave('UMAP_Mef2c_E85_Harmony_res0.4_labeled.pdf', 
       plot = UMAP_E85_harm_res0.4_labeled + 
         labs(title = 'Mef2c E8.5 Harmony, 28373 cells, res 0.4') +
         theme(plot.title = element_text(size=12)),
       device = 'pdf', width = 5, height = 5, dpi = 300
)

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
setwd("~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Outs/E85_v3")
ggsave('DotPlot_mef2c_E85_harmony.pdf', plot = DP_E85, device = 'pdf', width = 18, height = 6, dpi = 300)


###-----------------Subset----------------------------------------------
#Subset clusters of interest, namely CMs1, CMs2/CPs, SHF1/MixM, SHF2/PhM, SHF3/LPM
Idents(mef2c_v13_E85_v3) <- "harmony_cell_type"
mef2c_v13_E85_v3_subset <- subset(mef2c_v13_E85_v3, idents = c("CMs1", "CMs2/CPs", "SHF1/MixM", "SHF2/PhM", "SHF3/LPM"))
mef2c_v13_E85_v3_subset
#The subset contains 9240 cells

save(mef2c_v13_E85_v3, file = "~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Objects/mef2c_v13_E85_scTrans_v3_harmony.Robj")
save(mef2c_v13_E85_v3_subset, file = "~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Objects/mef2c_v13_E85_scTrans_v3_harmony_subset.Robj")


###---------Subset -Run scTransform and Dim Reduction---------------------------
mef2c_v13_E85_v3_subset <- SCTransform(mef2c_v13_E85_v3_subset, 
                                        method = "glmGamPoi", 
                                        variable.features.n = 5000,
                                        vars.to.regress = c("percent.mt", "percent.rib"), 
                                        verbose = TRUE,
                                        vst.flavor = "v2")

mef2c_v13_E85_v3_subset <- RunPCA(mef2c_v13_E85_v3_subset, assay = "SCT")
elb <- ElbowPlot(mef2c_v13_E85_v3_subset, ndims = 40)
elb


###------Subset - Clustering and UMAPs with Harmony Batch Correction--------------
mef2c_v13_E85_v3_subset <- RunHarmony(mef2c_v13_E85_v3_subset, group.by.vars = "Sample_Name", assay.use = "SCT", dims.use = 1:30)
mef2c_v13_E85_v3_subset <- FindNeighbors(mef2c_v13_E85_v3_subset, reduction = "harmony", dims = 1:30)
mef2c_v13_E85_v3_subset <- FindClusters(mef2c_v13_E85_v3_subset, resolution = 1, random.seed = 0)
mef2c_v13_E85_v3_subset <- RunUMAP(mef2c_v13_E85_v3_subset, reduction = "harmony", assay = "SCT", dims = 1:30, seed.use = 42)

#Visualize UMAP with sample labels
UMAP_E85_sub_Harmony_bysample <- DimPlot(object = mef2c_v13_E85_v3_subset, 
                                         reduction = "umap", 
                                         group.by = "Sample_Name", 
                                         cols = c("#F3E5F5", "#D1C4E9", "#9C27B0", "#4527A0"), 
                                         shuffle = T, seed = 1, pt.size = 0.1
                                         ) + NoLegend()
UMAP_E85_sub_Harmony_bysample_legend <- DimPlot(object = mef2c_v13_E85_v3_subset, 
                                                reduction = "umap", 
                                                group.by = "Sample_Name", 
                                                cols = c("#F3E5F5", "#D1C4E9", "#9C27B0", "#4527A0"), 
                                                shuffle = T, seed = 1, pt.size = 0.1
                                                ) 
setwd("~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Outs/E85_v3")
ggsave('UMAP_Mef2c_E85_Harmony_subset_bysample.pdf', plot = UMAP_E85_sub_Harmony_bysample + labs(title = 'Mef2c E8.5 Subset Harmony, 9240 cells'), device = 'pdf', width = 5, height = 5, dpi = 300)
ggsave('UMAP_Mef2c_E85_Harmony_subset_bysample_legend.pdf', plot = UMAP_E85_sub_Harmony_bysample_legend + labs(title = 'Mef2c E8.5 Subset Harmony, 9240 cells'), device = 'pdf', width = 5, height = 5, dpi = 300)

#Visualize UMAP with cluster labels
Idents(mef2c_v13_E85_v3_subset) <- "SCT_snn_res.1"
UMAP_E85_sub_Harmony_res1 <- DimPlot(object = mef2c_v13_E85_v3_subset, 
                                     reduction = "umap", 
                                     label = T, repel = F, raster = F, 
                                     label.size = 5, pt.size = 0.1
                                     ) + NoLegend()
setwd("~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Outs/E85_v3")
ggsave('UMAP_Mef2c_E85_Harmony_subset_res1.pdf', 
       plot = UMAP_E85_sub_Harmony_res1 + 
       labs(title = 'Mef2c E8.5 Subset Harmony, 9240 cells, res 1') +
       theme(plot.title = element_text(size=10)), 
       device = 'pdf', width = 5, height = 5, dpi = 300
       )


###-------Subset - Identify Clusters---------------------------------------
#Make feature plots for genes used to identify clusters 
#Create list of of marker genes
features3 = c(
  "Ttn", "Tnnt2", "Myl7",                         #CMs
  "Nppa", "Myl2", "Irx4",                         #Ventricles
  "Hand1", "Cited1",                              #LV
  "Cck", "Irx2", "Pln",                           #RV
  "Tbx5", "Angpt1", "Vsnl1", "Nr2f2",             #Atria
  "Rspo3", "Msx2", "Bmp2", "Tbx2",                #AVC
  "Wnt11", "Rgs5", "Sema3c", "Nrp2", "Tdgf1",     #OFT
  "Sfrp5", "Wnt2",                                #IFT
  "Wt1", "Tbx18",                                 #Proepicardium
  "Isl1", "Tbx1", "Fgf8", "Fgf10",                #aSHF
  "Aldh1a2", "Foxf1", "Gata6",                    #pSHF
  "Pdgfra", "Dlk1", "Bmp4", "Pbx1",               #Lateral Meso
  "Lef1"                                          #Posterior Meso
)

#Use for loop to save plots
setwd("~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Outs/E85_v3/Feature_Plots_Subset_harmony")
for (gene in features3) {
  gene_featureplot <- plot(
    FeaturePlot(mef2c_v13_E85_v3_subset, features = gene)
  )
  ggsave(file = paste0("FP_Mef2c_E85_subset_", gene, ".pdf"), plot = gene_featureplot, device = "pdf", width = 5, height = 5, dpi = 300)
}

#Run find markers for additional unbiased cluster markers
#Should return markers for 17 clusters, corresponding to res 1
Idents(mef2c_v13_E85_v3_subset) <- "SCT_snn_res.1"
mef2c_E85_sub.markers <- FindAllMarkers(mef2c_v13_E85_v3_subset, assay = "SCT", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.6)
mef2c_E85_sub.markers <- subset(mef2c_E85_sub.markers, p_val_adj < 0.05)
mef2c_E85_sub.markers <- mef2c_E85_sub.markers[order(mef2c_E85_sub.markers$avg_log2FC, decreasing = TRUE), ]
#Save to spreadsheet
setwd("~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Outs/E85_v3")
fwrite(mef2c_E85_sub.markers, row.names = TRUE, file = "mef2c_v13_E85_v3_markers_subset_res1.csv")

#Determined cell type labels for each cluster:
#C0 = CMs-V
#C1 = pSHF
#C2 = MixM1
#C3 = PhM1
#C4 = MixM2
#C5 = CMs-IFT1
#C6 = PhM2
#C7 = CMs-OFT
#C8 = LPM1
#C9 = PhM3
#C10 = CMs-IFT2
#C11 = aSHF
#C12 = LPM2
#C13 = CMs-IFT3
#C14 = Posterior Mesoderm 
#C15 = MixM3
#C16 = C16 - likley doublets

#Now label clusters by cell type
Idents(mef2c_v13_E85_v3_subset) <- "SCT_snn_res.1"
new.cluster.ids <- c("CMs-V", "pSHF", "MixM1", "PhM1",
                     "MixM2", "CMs-IFT1", "PhM2", "CMs-OFT",
                     "LPM1", "PhM3", "CMs-IFT2", "aSHF",
                     "LPM2", "CMs-IFT3", "PostM", "MixM3",
                     "C16")
names(new.cluster.ids) <- levels(mef2c_v13_E85_v3_subset$SCT_snn_res.1)
mef2c_v13_E85_v3_subset <- RenameIdents(mef2c_v13_E85_v3_subset, new.cluster.ids)

#Add "harmony_cell_type_subset" metadata column 
mef2c_v13_E85_v3_subset <- AddMetaData(mef2c_v13_E85_v3_subset, mef2c_v13_E85_v3_subset@active.ident, "harmony_cell_type_subset")

#Plot and save labeled UMAP
Idents(mef2c_v13_E85_v3_subset) <- "harmony_cell_type_subset"
UMAP_E85_harm_sub_res1_labeled <- DimPlot(object = mef2c_v13_E85_v3_subset, 
                                          reduction = "umap", 
                                          label = T, repel = T, 
                                          label.size = 4, pt.size = 0.1
                                          ) + NoLegend()
setwd("~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Outs/E85_v3")
ggsave('UMAP_Mef2c_E85_Harmony_subset_res1_labeled.pdf', 
       plot = UMAP_E85_harm_sub_res1_labeled + 
         labs(title = 'Mef2c E8.5 Harmony Subset, 9240 cells, res 1') +
         theme(plot.title = element_text(size=12)),
       device = 'pdf', width = 5, height = 5, dpi = 300
)

#Plot and save UMAP with cells labeled by full dataset cell types to double-check that cell type labels make sense
Idents(mef2c_v13_E85_v3_subset) <- "harmony_cell_type"
UMAP_E85_harm_sub_labeled2 <- DimPlot(mef2c_v13_E85_v3_subset, reduction = "umap", label = T, repel = F, label.size = 4, pt.size = 0.1) + NoLegend()
setwd("~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Outs/E85_v3")
ggsave('UMAP_Mef2c_E85_Harmony_subset_labeled_fulldatalabels.pdf', 
       plot = UMAP_E85_harm_sub_labeled2 + 
         labs(title = 'Mef2c E8.5 Harmony Subset, 9240 cells') +
         theme(plot.title = element_text(size=12)),
       device = 'pdf', width = 5, height = 5, dpi = 300
)

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
DP_E85_sub <- DotPlot(mef2c_v13_E85_v3_subset, features = features2) + RotatedAxis()
setwd("~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Outs/E85_v3")
ggsave('DotPlot_mef2c_E85_harmony_subset.pdf', plot = DP_E85_sub, device = 'pdf', width = 18, height = 6, dpi = 300)


###----------Subset - DEG Testing ------------------------------------------

#Test for DEG between KO and WT cells in IFT CMs clusters (3 clusters)
Idents(mef2c_v13_E85_v3_subset) <- "harmony_cell_type_subset"
markers_mef2c_E85_KOvWT_IFTCM <- FindMarkers(mef2c_v13_E85_v3_subset, asay = "SCT", ident.1 = "KO", ident.2 = "WT", group.by = "Genotype", subset.ident = c("CMs-IFT1", "CMs-IFT2", "CMs-IFT3"))
markers_mef2c_E85_KOvWT_IFTCM <- subset(markers_mef2c_E85_KOvWT_IFTCM, p_val_adj < 0.05)
markers_mef2c_E85_KOvWT_IFTCM <- markers_mef2c_E85_KOvWT_IFTCM[order(markers_mef2c_E85_KOvWT_IFTCM$avg_log2FC, decreasing = TRUE), ]
names(markers_mef2c_E85_KOvWT_IFTCM)[3] <- "pct.KO"
names(markers_mef2c_E85_KOvWT_IFTCM)[4] <- "pct.WT"
setwd("~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Outs/E85_v3")
fwrite(markers_mef2c_E85_KOvWT_IFTCM, row.names = TRUE, file = "Mef2c_v13_E85_harmony_subset_markers_KOvWT_IFTCM.csv")

#Test for DEG between KO and WT cells in V CMs cluster
Idents(mef2c_v13_E85_v3_subset) <- "harmony_cell_type_subset"
markers_mef2c_E85_KOvWT_VCM <- FindMarkers(mef2c_v13_E85_v3_subset, assay = "SCT", ident.1 = "KO", ident.2 = "WT", group.by = "Genotype", subset.ident = "CMs-V")
markers_mef2c_E85_KOvWT_VCM <- subset(markers_mef2c_E85_KOvWT_VCM, p_val_adj < 0.05)
markers_mef2c_E85_KOvWT_VCM <- markers_mef2c_E85_KOvWT_VCM[order(markers_mef2c_E85_KOvWT_VCM$avg_log2FC, decreasing = TRUE), ]
names(markers_mef2c_E85_KOvWT_VCM)[3] <- "pct.KO"
names(markers_mef2c_E85_KOvWT_VCM)[4] <- "pct.WT"
setwd("~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Outs/E85_v3")
fwrite(markers_mef2c_E85_KOvWT_VCM, row.names = TRUE, file = "Mef2c_v13_E85_harmony_subset_markers_KOvWT_VCM.csv")

#Test for DEG between KO and WT cells in OFT CMs cluster
Idents(mef2c_v13_E85_v3_subset) <- "harmony_cell_type_subset"
markers_mef2c_E85_KOvWT_OFTCM <- FindMarkers(mef2c_v13_E85_v3_subset, asay = "SCT", ident.1 = "KO", ident.2 = "WT", group.by = "Genotype", subset.ident = "CMs-OFT")
markers_mef2c_E85_KOvWT_OFTCM <- subset(markers_mef2c_E85_KOvWT_OFTCM, p_val_adj < 0.05)
markers_mef2c_E85_KOvWT_OFTCM <- markers_mef2c_E85_KOvWT_OFTCM[order(markers_mef2c_E85_KOvWT_OFTCM$avg_log2FC, decreasing = TRUE), ]
names(markers_mef2c_E85_KOvWT_OFTCM)[3] <- "pct.KO"
names(markers_mef2c_E85_KOvWT_OFTCM)[4] <- "pct.WT"
setwd("~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Outs/E85_v3")
fwrite(markers_mef2c_E85_KOvWT_OFTCM, row.names = TRUE, file = "Mef2c_v13_E85_harmony_subset_markers_KOvWT_OFTCM.csv")


###----------------Save Robjects-----------------------------------------
save(mef2c_v13_E85_v3, file = "~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Objects/mef2c_v13_E85_scTrans_v3_harmony.Robj")
save(mef2c_v13_E85_v3_subset, file = "~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Objects/mef2c_v13_E85_scTrans_v3_harmony_subset.Robj")
