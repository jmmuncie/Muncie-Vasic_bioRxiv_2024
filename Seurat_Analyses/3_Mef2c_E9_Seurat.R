#This script is a starting point for analyzing scRNA-seq data with Seurat using
#scTransform for normalization and Harmony for batch correction. 
#Specifically, this script is for analyzing the 4 samples from the Mef2c Multiome
#E9 timepoint

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
sample05.data <- Read10X(data.dir = "~/Desktop/Mef2c_v13_Seurat_scTrans_working/Inputs/sample05")
sample05 <- CreateSeuratObject(counts = sample05.data[[1]], project = "mef2c_v13_E775_sample05")
sample06.data <- Read10X(data.dir = "~/Desktop/Mef2c_v13_Seurat_scTrans_working/Inputs/sample06")
sample06 <- CreateSeuratObject(counts = sample06.data[[1]], project = "mef2c_v13_E775_sample06")
sample07.data <- Read10X(data.dir = "~/Desktop/Mef2c_v13_Seurat_scTrans_working/Inputs/sample07")
sample07 <- CreateSeuratObject(counts = sample07.data[[1]], project = "mef2c_v13_E775_sample07")
sample08.data <- Read10X(data.dir = "~/Desktop/Mef2c_v13_Seurat_scTrans_working/Inputs/sample08")
sample08 <- CreateSeuratObject(counts = sample08.data[[1]], project = "mef2c_v13_E775_sample08")

mef2c_v13_E9 <- merge(sample05, 
                        y = c(sample06, sample07, sample08), 
                        add.cell.ids = c("sample05", "sample06", "sample07", "sample08"), 
                        project = "mef2c_v13_E9")

save(mef2c_v13_E9, file = "~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Objects/mef2c_v13_E9_scTrans.Robj")


###---------Add Meta Data Info----------------------------------------------

#Add sample number metadata column
class.vec1 <- gsub("_.*","", (colnames(mef2c_v13_E9@assays$RNA@data)))
names(class.vec1) <- colnames(mef2c_v13_E9@assays$RNA@data)
mef2c_v13_E9 <- AddMetaData(mef2c_v13_E9, class.vec1, "Sample_Number")

#Add sample name metadata column
class.vec1[class.vec1=="sample05"] <- "E9.0_WT1"
class.vec1[class.vec1=="sample06"] <- "E9.0_KO1"
class.vec1[class.vec1=="sample07"] <- "E9.0_KO2"
class.vec1[class.vec1=="sample08"] <- "E9.0_WT2"
mef2c_v13_E9 <- AddMetaData(mef2c_v13_E9, class.vec1, "Sample_Name")

#Add genotype metadata column
class.vec1[class.vec1 == "E9.0_WT1" | class.vec1 == "E9.0_WT2"] <- "WT"
class.vec1[class.vec1 == "E9.0_KO1" | class.vec1 == "E9.0_KO2"] <- "KO"
mef2c_v13_E9 <- AddMetaData(mef2c_v13_E9, class.vec1, "Genotype")

#Add percent mitochondria metadata column
mef2c_v13_E9 <- PercentageFeatureSet(mef2c_v13_E9, pattern = "^mt-", col.name = "percent.mt")

#Add perecent ribosomal reads to metadata column
#Initially, scTransform pipeline led to a cluster that was marked entirely by 
#ribosomal transcripts. Creating this metadata column so I can regress on them
mef2c_v13_E9 <- PercentageFeatureSet(mef2c_v13_E9, pattern = "^Rp[ls]", col.name = "percent.rib")

save(mef2c_v13_E9, file = "~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Objects/mef2c_v13_E9_scTrans.Robj")


###---------QC and Filtering----------------------------------------------

#Visualize QC metrics, grouped by sample, as a violin plot
# Note: the use of factor with set levels in the second line below will set the 
#order the groups are displayed any time that you use split.by Sample_Name
mef2c_v13_E9$Sample_Name <- factor(mef2c_v13_E9$Sample_Name, levels = c("E9.0_WT1", "E9.0_WT2", "E9.0_KO1", "E9.0_KO2"))
Idents(mef2c_v13_E9) <- "Sample_Name"
vlp1 <- VlnPlot(mef2c_v13_E9, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, split.by = "Sample_Name", pt.size = 0.1)
setwd("~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Outs/E9")
ggsave("QC_ViolinPlots_Mef2c_E9_prefilters.pdf", plot = vlp1, device = "pdf", width = 11, height = 5.5, dpi = 300)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
Idents(mef2c_v13_E9) <- "Sample_Name"
plot1 <- FeatureScatter(mef2c_v13_E9, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "Sample_Name", shuffle = T, pt.size = 0.5)
plot2 <- FeatureScatter(mef2c_v13_E9, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "Sample_Name", shuffle = T, pt.size = 0.5)
plot3 <- FeatureScatter(mef2c_v13_E9, feature1 = "percent.rib", feature2 = "percent.mt", group.by = "Sample_Name", shuffle = T, pt.size = 0.5)
fsc <- plot1+plot2+plot3
setwd("~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Outs/E9")
ggsave("QC_FeatureScatters_Mef2c_E9_prefilters.pdf", plot = fsc, device = "pdf", width = 15, height = 4, dpi = 300)

# Now filtering out GEMs with:
#  < 1000 features or > 5000 features
#  > 100000 counts
#  > 20% mt
#  > 20% ribosomal content
mef2c_v13_E9 <- subset(mef2c_v13_E9, subset = nFeature_RNA > 1000 & nFeature_RNA < 5000 & nCount_RNA < 100000 & percent.mt < 20 & percent.rib < 20)
#Removed ~5000 cells (52768 to 47482)

# Visualize QC metrics for filtered dataset, grouped by sample, as a violin plot
Idents(mef2c_v13_E9) <- "Sample_Name"
vlp2 <- VlnPlot(mef2c_v13_E9, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rib"), ncol = 4, split.by = "Sample_Name", pt.size = 0.1)
plot4 <- FeatureScatter(mef2c_v13_E9, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "Sample_Name", shuffle = T, pt.size = 0.5)
plot5 <- FeatureScatter(mef2c_v13_E9, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "Sample_Name", shuffle = T, pt.size = 0.5)
plot6 <- FeatureScatter(mef2c_v13_E9, feature1 = "percent.rib", feature2 = "percent.mt", group.by = "Sample_Name", shuffle = T, pt.size = 0.5)
fsc2 <- plot4+plot5+plot6
setwd("~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Outs/E9")
ggsave("QC_ViolinPlots_Mef2c_E9_postfilters.pdf", plot = vlp2, device = "pdf", width = 15, height = 5.5, dpi = 300)
ggsave("QC_FeatureScatters_Mef2c_E9_postfilters.pdf", plot = fsc2, device = "pdf", width = 15, height = 4, dpi = 300)


###------------------Run scTransform and Dim Reduction---------------------------
#Replaces/combines normalization, scaling, and finding variable features

# Note: glmGamPoi speeds up scTransform - if didn't load properly, just remove that argument
# Note: as in standard workflow, only the variable genes are saved to the scale.data
# matrices in the output assay 

mef2c_v13_E9 <- SCTransform(mef2c_v13_E9, 
                              method = "glmGamPoi", 
                              variable.features.n = 5000,
                              vars.to.regress = c("percent.mt", "percent.rib"), 
                              verbose = TRUE,
                              vst.flavor = "v2")

mef2c_v13_E9 <- RunPCA(mef2c_v13_E9, assay = "SCT")
elb <- ElbowPlot(mef2c_v13_E9, ndims = 40)
elb


###------Clustering and UMAPs with Harmony Batch Correction--------------
mef2c_v13_E9 <- RunHarmony(mef2c_v13_E9, group.by.vars = "Sample_Name", assay.use = "SCT", dims.use = 1:30)
mef2c_v13_E9 <- FindNeighbors(mef2c_v13_E9, reduction = "harmony", dims = 1:30)
mef2c_v13_E9 <- FindClusters(mef2c_v13_E9, resolution = 0.6, random.seed = 0)
mef2c_v13_E9 <- RunUMAP(mef2c_v13_E9, reduction = "harmony", assay = "SCT", dims = 1:30, seed.use = 42)

#Visualize UMAP with sample labels
UMAP_E9_Harmony_bysample <- DimPlot(object = mef2c_v13_E9, 
                                    reduction = "umap", 
                                    group.by = "Sample_Name", 
                                    cols = c("pink", "plum1", "red", "red4"), 
                                    shuffle = T, seed = 1, pt.size = 0.1
                                    ) + NoLegend()
UMAP_E9_Harmony_bysample_legend <- DimPlot(object = mef2c_v13_E9, 
                                           reduction = "umap", 
                                           group.by = "Sample_Name", 
                                           cols = c("pink", "plum1", "red", "red4"), 
                                           shuffle = T, seed = 1, pt.size = 0.1
                                           )
setwd("~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Outs/E9")
ggsave('UMAP_Mef2c_E9_Harmony_bysample.pdf', plot = UMAP_E9_Harmony_bysample + labs(title = 'Mef2c E9.0, Harmony corrected, 47482 cells'), device = 'pdf', width = 5, height = 5, dpi = 300)
ggsave('UMAP_Mef2c_E9_Harmony_bysample_legend.pdf', plot = UMAP_E9_Harmony_bysample_legend + labs(title = 'Mef2c E9.0, Harmony corrected, 47482 cells'), device = 'pdf', width = 5, height = 5, dpi = 300)

#Visualize UMAP with cluster labels
Idents(mef2c_v13_E9) <- "SCT_snn_res.0.6"
UMAP_E9_Harmony_res0.6 <- DimPlot(mef2c_v13_E9, reduction = "umap", label = T, repel = F, raster = F, label.size = 5, pt.size = 0.1) + NoLegend()
setwd("~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Outs/E9")
ggsave('UMAP_Mef2c_E9_Harmony_res0.6.pdf', plot = UMAP_E9_Harmony_res0.6 + labs(title = 'Mef2c E9.0, Harmony corrected, 47482 cells, res 0.6'), device = 'pdf', width = 5, height = 5, dpi = 300)


###----------Identify Clusters---------------------------------------
#Make feature plots for genes used to identify clusters 
#Create list of of marker genes
features1 = c(
  "Pdgfra", "Foxf1", "Alx1", "Gata6",      #Lateral Plate Mesoderm
  "Dlk1", "Bmp5",                          #Mixed Mesoderm
  "Foxf1", "Foxp2",                        #Splanchnic Mesoderm
  "Sox2", "Fgf14", "Nrcam",                #Neural Tube
  "Meox1", "Fst", "Cdh11",                 #Somitic Mesoderm
  "Sox17","Epcam", "Cdh1", "Cldn6",        #Definitive Endoderm
  "Isl1", "Tbx1", "Fgf8", "Fgf10",         #SHF
  "Nkx2-5", "Tbx5", "Mef2c",               #Cardiac markers
  "Ttn", "Tnnt2", "Myl7",                  #CMs
  "Tbx18", "Wt1",                          #Proepicardium
  "Tfap2b", "Sox10", "Ets1",               #Neural Crest
  "Cdh5", "Kdr", "Flt1",                   #Endothelium
  "Nfatc1",                                #Endocardium
  "Postn", "Col1a2", "Col3a1",             #Mesenchyme
  "Hba-a1", "Hba-a2", "Hbb-y",             #Blood
  "Ahnak", "Podxl",                        #ExMeso
  "Ttr", "Afp", "Apoa1",                   #Yolk Sac
  "T", "Shh"                               #Notochord
)

#Use for loop to save plots
setwd("~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Outs/E9/Feature_Plots")
for (gene in features1) {
  gene_featureplot <- plot(
    FeaturePlot(mef2c_v13_E9, features = gene)
  )
  ggsave(file = paste0("FP_Mef2c_E9_", gene, ".pdf"), plot = gene_featureplot, device = "pdf", width = 5, height = 5, dpi = 300)
}

#Run find markers for additional unbiased cluster markers
#Should return markers for 22 clusters, corresponding to res 0.6
Idents(mef2c_v13_E9) <- "SCT_snn_res.0.6"
mef2c_E9.markers <- FindAllMarkers(mef2c_v13_E9, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.6)
mef2c_E9.markers <- subset(mef2c_E9.markers, p_val_adj < 0.05)
mef2c_E9.markers <- mef2c_E9.markers[order(mef2c_E9.markers$avg_log2FC, decreasing = TRUE), ]
#Save to spreadsheet
setwd("~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Outs/E9")
fwrite(mef2c_E9.markers, row.names = TRUE, file = "Mef2c_v13_E9_markers_res0.6.csv")

#Determined cell type labels for each cluster:
#C0 = Endoderm 1
#C1 = MixM
#C2 = Neural Tube 1
#C3 = Somitic Mesoderm 1
#C4 = CMs
#C5 = Neural Tube 2
#C6 = SHF 
#C7 = Neural Tube 3
#C8 = Proepicardium
#C9 = Neural Crest
#C10 = Splanchnic Mesoderm (Foxf1, Foxp2) 
#C11 = Somitic Mesoderm 2
#C12 = Endothelium (includes Nfatc1+ epicardium)
#C13 = Endoderm 2
#C14 = Endoderm 3
#C15 = Neural Tube 4
#C16 = Mesenchyme (Postn, Col1a2, Col3a1)
#C17 = Endoderm 4
#C18 = Blood
#C19 = Endoderm 5 
#C20 = Yolk Sac
#C21 = Notochord

#Now label clusters by cell type
Idents(mef2c_v13_E9) <- "SCT_snn_res.0.6"
new.cluster.ids <- c("En1", "MixM", "NT1", "SoM1", 
                     "CMs", "NT2", "SHF", "NT3",
                     "Pe", "NC", "SpM", "SoM2",
                     "Endo", "En2", "En3", "NT4",
                     "Me", "En4", "Blood", "En5",
                     "YS", "Noto"
                     )
names(new.cluster.ids) <- levels(mef2c_v13_E9$SCT_snn_res.0.6)
mef2c_v13_E9 <- RenameIdents(mef2c_v13_E9, new.cluster.ids)

#Add "cell_type" metadata column 
mef2c_v13_E9 <- AddMetaData(mef2c_v13_E9, mef2c_v13_E9@active.ident, "cell_type")

#Plot and save labeled UMAP
Idents(mef2c_v13_E9) <- "cell_type"
UMAP_E9_Harmony_res0.6_labeled <- DimPlot(object = mef2c_v13_E9, 
                                          reduction = "umap", 
                                          label = T, repel = F, 
                                          label.size = 4, pt.size = 0.1
                                          ) + NoLegend()
setwd("~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Outs/E9")
ggsave('UMAP_Mef2c_E9_Harmony_res0.6_labeled.pdf', 
       plot = UMAP_E9_Harmony_res0.6_labeled + 
       labs(title = 'Mef2c E9.0, Harmony corrected, 47482 cells, res 0.6') +
       theme(plot.title = element_text(size=12)),
       device = 'pdf', width = 5, height = 5, dpi = 300
       )

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
features2 = c(
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
DP_E9 <- DotPlot(mef2c_v13_E9, features = features2) + RotatedAxis()
setwd("~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Outs/E9")
ggsave('DotPlot_mef2c_E9.pdf', plot = DP_E9, device = 'pdf', width = 19, height = 7, dpi = 300)


###-----------------Subset----------------------------------------------
#Subset clusters of interest, namely CMs, Pe, SHF 
Idents(mef2c_v13_E9) <- "cell_type"
mef2c_v13_E9_subset <- subset(mef2c_v13_E9, idents = c("CMs", "Pe", "SHF"))
#The subset contains 9259 cells

save(mef2c_v13_E9, file = "~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Objects/mef2c_v13_E9_scTrans.Robj")
save(mef2c_v13_E9_subset, file = "~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Objects/mef2c_v13_E9_subset_scTrans.Robj")


###-------------Subset - Run scTransform and Dim Reduction---------------------------
mef2c_v13_E9_subset <- SCTransform(mef2c_v13_E9_subset, 
                                     method = "glmGamPoi", 
                                     variable.features.n = 5000,
                                     vars.to.regress = c("percent.mt", "percent.rib"), 
                                     verbose = TRUE,
                                     vst.flavor = "v2")

mef2c_v13_E9_subset <- RunPCA(mef2c_v13_E9_subset, assay = "SCT")
elb <- ElbowPlot(mef2c_v13_E9_subset, ndims = 40)
elb


###------Subset - Clustering and UMAPs with Harmony Batch Correction--------------
mef2c_v13_E9_subset <- RunHarmony(mef2c_v13_E9_subset, group.by.vars = "Sample_Name", assay.use = "SCT", dims.use = 1:30)
mef2c_v13_E9_subset <- FindNeighbors(mef2c_v13_E9_subset, reduction = "harmony", dims = 1:30)
mef2c_v13_E9_subset <- FindClusters(mef2c_v13_E9_subset, resolution = 0.6, random.seed = 0)
mef2c_v13_E9_subset <- RunUMAP(mef2c_v13_E9_subset, reduction = "harmony", assay = "SCT", dims = 1:30, seed.use = 42)

#Visualize UMAP with sample labels
UMAP_E9_subset_Harmony_bysample <- DimPlot(object = mef2c_v13_E9_subset, 
                                           reduction = "umap", 
                                           group.by = "Sample_Name", 
                                           cols = c("pink", "plum1", "red", "red4"), 
                                           shuffle = T, seed = 1, pt.size = 0.5
                                           ) + NoLegend()
setwd("~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Outs/E9")
ggsave('UMAP_Mef2c_E9_subset_Harmony_bysample.pdf', 
       plot = UMAP_E9_subset_Harmony_bysample + 
       labs(title = 'Mef2c E9.0 Subset, Harmony corrected, 9259 cells') +
       theme(plot.title = element_text(size=12)),
       device = 'pdf', width = 5, height = 5, dpi = 300
       )

#Visualize UMAP with cluster labels
Idents(mef2c_v13_E9_subset) <- "SCT_snn_res.0.6"
UMAP_E9_subset_Harmony_res0.6 <- DimPlot(object = mef2c_v13_E9_subset, 
                                         reduction = "umap", 
                                         label = T, raster = F, 
                                         label.size = 5, pt.size = 0.5
                                         ) + NoLegend()
setwd("~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Outs/E9")
ggsave('UMAP_Mef2c_E9_subset_Harmony_res0.6.pdf', 
       plot = UMAP_E9_subset_Harmony_res0.6 + 
       labs(title = 'Mef2c E9.0 Subset, Harmony corrected, 9259 cells, res 0.6') +
       theme(plot.title = element_text(size=12)),
       device = 'pdf', width = 5, height = 5, dpi = 300
       )


###----------Subset - Identify Clusters---------------------------------------

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
  "Sfrp5", "Wnt2",                                #Venous Pole (VP)
  "Wt1", "Tbx18",                                 #Proepicardium
  "Isl1", "Fgf8", "Fgf10",                        #SHF
  "Pdgfra", "Dlk1", "Bmp5", "Pbx1",               #Lateral Meso
  "Ebf1", "Tbx1", "Cdh11",                        #PhM
  "Fgf14", "Nrcam", "Dcc"                         #Neural adhesion genes that show up in some clusters
)

#Use for loop to save plots
setwd("~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Outs/E9/Feature_Plots_Subset")
for (gene in features3) {
  gene_featureplot <- plot(
    FeaturePlot(mef2c_v13_E9_subset, features = gene)
  )
  ggsave(file = paste0("FP_Mef2c_E9_subset_", gene, ".pdf"), plot = gene_featureplot, device = "pdf", width = 5, height = 5, dpi = 300)
}

#Run find markers for additional unbiased cluster markers
#Should return markers for 12 clusters, corresponding to res 0.6
Idents(mef2c_v13_E9_subset) <- "SCT_snn_res.0.6"
mef2c_E9_subset.markers <- FindAllMarkers(mef2c_v13_E9_subset, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.7)
mef2c_E9_subset.markers <- subset(mef2c_E9_subset.markers, p_val_adj < 0.05)
mef2c_E9_subset.markers <- mef2c_E9_subset.markers[order(mef2c_E9_subset.markers$avg_log2FC, decreasing = TRUE), ]
#Save to spreadsheet
setwd("~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Outs/E9")
fwrite(mef2c_E9_subset.markers, row.names = TRUE, file = "Mef2c_v13_E9_subset_markers_res0.6.csv")

#Determined cell type labels for each cluster: 
#C0 = SHF1
#C1 = CMs - Atrial 1
#C2 = CMs - Ventricular
#C3 = Venous Pole 1
#C4 = CMs - AVC
#C5 = SHF2
#C6 = CMs - OFT
#C7 = PhM
#C8 = CMs - Atrial 2
#C9 = Pe
#C10 = Venous Pole 2
#C11 = C11 (likely doublets)

#Now label clusters by cell type
Idents(mef2c_v13_E9_subset) <- "SCT_snn_res.0.6"
new.cluster.ids <- c("SHF1", "CMs-A1", "CMs-V", "VP1", "CMs-AVC",
                     "SHF2", "CMs-OFT", "PhM", "CMs-A2", "Pe",
                     "VP2", "C11"
)
names(new.cluster.ids) <- levels(mef2c_v13_E9_subset$SCT_snn_res.0.6)
mef2c_v13_E9_subset <- RenameIdents(mef2c_v13_E9_subset, new.cluster.ids)

#Add "cell_type_subset" metadata column 
mef2c_v13_E9_subset <- AddMetaData(mef2c_v13_E9_subset, mef2c_v13_E9_subset@active.ident, "cell_type_subset")

#Plot and save labeled UMAP
Idents(mef2c_v13_E9_subset) <- "cell_type_subset"
UMAP_E9_subset_Harmony_res0.6_labeled <- DimPlot(object = mef2c_v13_E9_subset, 
                                                 reduction = "umap", 
                                                 label = T, repel = T, 
                                                 label.size = 4, pt.size = 0.5
                                                 ) + NoLegend()
setwd("~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Outs/E9")
ggsave('UMAP_Mef2c_E9_subset_Harmony_res0.6_labeled.pdf', 
       plot = UMAP_E9_subset_Harmony_res0.6_labeled + 
       labs(title = 'Mef2c E9.0 Subset, Harmony corrected, 9259 cells, res 0.6') +
       theme(plot.title = element_text(size=10)),
       device = 'pdf', width = 5, height = 5, dpi = 300
       )

#Plot and save UMAP labeled by full dataset labels to double check that cell type annotations make sense
Idents(mef2c_v13_E9_subset) <- "cell_type"
UMAP_E9_subset_Harmony_labeled_full <- DimPlot(object = mef2c_v13_E9_subset, 
                                               reduction = "umap", 
                                               label = T, repel = T, 
                                               label.size = 4, pt.size = 0.5
                                               ) + NoLegend()
setwd("~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Outs/E9")
ggsave('UMAP_Mef2c_E9_subset_Harmony_labeled_fulldataset.pdf', 
       plot = UMAP_E9_subset_Harmony_labeled_full + 
         labs(title = 'Mef2c E9.0 Subset, Harmony corrected, 9259 cells') +
         theme(plot.title = element_text(size=10)),
       device = 'pdf', width = 5, height = 5, dpi = 300
)

#Create cluster labels for dot plot
Idents(mef2c_v13_E9_subset) <- "SCT_snn_res.0.6"
new.cluster.ids <- c("C0-SHF1", "C1-CMs-A1", "C2-CMs-V", "C3-VP1", "C4-CMs-AVC",
                     "C5-SHF2", "C6-CMs-OFT", "C7-PhM", "C8-CMs-A2", "C9-Pe",
                     "C10-VP2", "C11"
)
names(new.cluster.ids) <- levels(mef2c_v13_E9_subset$SCT_snn_res.0.6)
mef2c_v13_E9_subset <- RenameIdents(mef2c_v13_E9_subset, new.cluster.ids)

#Plot and save dot plot of key markers 
features4 = c(
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
DP_E9_subset <- DotPlot(mef2c_v13_E9_subset, features = features4) + RotatedAxis()
setwd("~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Outs/E9")
ggsave('DotPlot_mef2c_E9_subset.pdf', plot = DP_E9_subset, device = 'pdf', width = 15, height = 6, dpi = 300)


###----------Subset - DEG Testing ------------------------------------------

#Test for DEG between KO and WT cells in A CMs clusters
Idents(mef2c_v13_E9_subset) <- "cell_type_subset"
markers_mef2c_E9_KOvWT_ACM <- FindMarkers(mef2c_v13_E9_subset, ident.1 = "KO", ident.2 = "WT", group.by = "Genotype", subset.ident = c("CMs-A1", "CMs-A2"))
markers_mef2c_E9_KOvWT_ACM <- subset(markers_mef2c_E9_KOvWT_ACM, p_val_adj < 0.05)
markers_mef2c_E9_KOvWT_ACM <- markers_mef2c_E9_KOvWT_ACM[order(markers_mef2c_E9_KOvWT_ACM$avg_log2FC, decreasing = TRUE), ]
names(markers_mef2c_E9_KOvWT_ACM)[3] <- "pct.KO"
names(markers_mef2c_E9_KOvWT_ACM)[4] <- "pct.WT"
setwd("~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Outs/E9")
fwrite(markers_mef2c_E9_KOvWT_ACM, row.names = TRUE, file = "Mef2c_v13_E9_subset_markers_KOvWT_ACM.csv")

#Test for DEG between KO and WT cells in AVC CMs cluster
Idents(mef2c_v13_E9_subset) <- "cell_type_subset"
markers_mef2c_E9_KOvWT_AVCCM <- FindMarkers(mef2c_v13_E9_subset, ident.1 = "KO", ident.2 = "WT", group.by = "Genotype", subset.ident = "CMs-AVC")
markers_mef2c_E9_KOvWT_AVCCM <- subset(markers_mef2c_E9_KOvWT_AVCCM, p_val_adj < 0.05)
markers_mef2c_E9_KOvWT_AVCCM <- markers_mef2c_E9_KOvWT_AVCCM[order(markers_mef2c_E9_KOvWT_AVCCM$avg_log2FC, decreasing = TRUE), ]
names(markers_mef2c_E9_KOvWT_AVCCM)[3] <- "pct.KO"
names(markers_mef2c_E9_KOvWT_AVCCM)[4] <- "pct.WT"
setwd("~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Outs/E9")
fwrite(markers_mef2c_E9_KOvWT_AVCCM, row.names = TRUE, file = "Mef2c_v13_E9_subset_markers_KOvWT_AVCCM.csv")

#Test for DEG between KO and WT cells in V CMs cluster
Idents(mef2c_v13_E9_subset) <- "cell_type_subset"
markers_mef2c_E9_KOvWT_VCM <- FindMarkers(mef2c_v13_E9_subset, ident.1 = "KO", ident.2 = "WT", group.by = "Genotype", subset.ident = "CMs-V")
markers_mef2c_E9_KOvWT_VCM <- subset(markers_mef2c_E9_KOvWT_VCM, p_val_adj < 0.05)
markers_mef2c_E9_KOvWT_VCM <- markers_mef2c_E9_KOvWT_VCM[order(markers_mef2c_E9_KOvWT_VCM$avg_log2FC, decreasing = TRUE), ]
names(markers_mef2c_E9_KOvWT_VCM)[3] <- "pct.KO"
names(markers_mef2c_E9_KOvWT_VCM)[4] <- "pct.WT"
setwd("~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Outs/E9")
fwrite(markers_mef2c_E9_KOvWT_VCM, row.names = TRUE, file = "Mef2c_v13_E9_subset_markers_KOvWT_VCM.csv")

#Test for DEG between KO and WT cells in OFT CMs cluster
Idents(mef2c_v13_E9_subset) <- "cell_type_subset"
markers_mef2c_E9_KOvWT_OFTCM <- FindMarkers(mef2c_v13_E9_subset, ident.1 = "KO", ident.2 = "WT", group.by = "Genotype", subset.ident = "CMs-OFT")
markers_mef2c_E9_KOvWT_OFTCM <- subset(markers_mef2c_E9_KOvWT_OFTCM, p_val_adj < 0.05)
markers_mef2c_E9_KOvWT_OFTCM <- markers_mef2c_E9_KOvWT_OFTCM[order(markers_mef2c_E9_KOvWT_OFTCM$avg_log2FC, decreasing = TRUE), ]
names(markers_mef2c_E9_KOvWT_OFTCM)[3] <- "pct.KO"
names(markers_mef2c_E9_KOvWT_OFTCM)[4] <- "pct.WT"
setwd("~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Outs/E9")
fwrite(markers_mef2c_E9_KOvWT_OFTCM, row.names = TRUE, file = "Mef2c_v13_E9_subset_markers_KOvWT_OFTCM.csv")


###----------Save Seurat Objects---------------------------------------
save(mef2c_v13_E9, file = "~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Objects/mef2c_v13_E9_scTrans.Robj")
save(mef2c_v13_E9_subset, file = "~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Objects/mef2c_v13_E9_subset_scTrans.Robj")

