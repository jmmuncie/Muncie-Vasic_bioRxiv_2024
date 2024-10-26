###----------------Description-------------------------------------

#This script is used to perform integrated analyses of snRNA-seq and snATAC-seq
#data in ArchR. Due to low numbers of MEF2C KO OFT-CMs, this analysis combines
#the eight samples from E8.5 and E9 to perform MEF2C KO-vs-WT comparisons for
#OFT-CMs

#Script written by Jonathon M. Muncie-Vasic
#edited 10/25/2024


###---------------Load Libraries and Set Threads------------------------------

#To install ArchR release_1.0.2 branch:
devtools::install_github("GreenleafLab/ArchR", ref="release_1.0.2", repos = BiocManager::repositories())

library(ArchR); library(Seurat); library(dplyr)
library(GenomeInfoDb); library(ensembldb)
library(ggplot2); library(patchwork)
set.seed(1234)
library(BSgenome) 
library(hexbin)

#Set ArchR Threads to 4:
addArchRThreads(threads = 4, force = FALSE)


###-------------Load Genome and Gene Annotations------------------------------------
#These are created in script "1_Mef2c_ArchR_create_arrows.R"
load("~/Genome_Gene_Annotations/genomeAnnotation_mm10v13.Robj")
load("~/Genome_Gene_Annotations/geneAnnotation_mm10v13.Robj")
setwd("~/Mef2c_ArchR_working")


###--------------Create ArchR Project----------------------
#Arrow Files previously created in script "1_Mef2c_ArchR_create_arrows.R"
#Move/copy arrow files to "~/Mef2c_ArchR_working/Original_Arrows"

#Create ArchR Project
proj_Mef2c_v13_E85_E9_ATAC_only <- ArchRProject(
  ArrowFiles = c("Original_Arrows/sample05.arrow", 
                 "Original_Arrows/sample06.arrow", 
                 "Original_Arrows/sample07.arrow", 
                 "Original_Arrows/sample08.arrow",
                 "Original_Arrows/sample21.arrow", 
                 "Original_Arrows/sample22.arrow", 
                 "Original_Arrows/sample23.arrow", 
                 "Original_Arrows/sample24.arrow"
                 ),
  outputDirectory = "Mef2c_v13_E85_E9_ATAC_only",
  copyArrows = TRUE,
  geneAnnotation = geneAnnotation_mm10v13,
  genomeAnnotation = genomeAnnotation_mm10v13,
  showLogo = FALSE
)


###-------------Add GEX data to ArchR project--------------------------

#Import scRNA data
#Filtered feature matrices in .h5 format that are output from Cellranger-arc counts
#should be located in "~/Mef2c_ArchR_working/Mef2c_inputs"

#Note: to do this on laptop, needed to set Renviron max memory to 64 Gb (includes virtual mem)
#$library(usethis) 
#$usethis::edit_r_environ()
#When tab opens, enter "R_MAX_VSIZE=64Gb" on first line, save, and restart R

#Initially had trouble importing 10x matrices from multiple samples
#This was resoled on Github Issue #507: https://github.com/GreenleafLab/ArchR/issues/507
#In brief: Cellranger seems to occasionally use mismatched in Ensembl IDs and ranges 
#for the same gene in different samples. This is fixed with the strictMatch
#input. When set to TRUE, it discards entries with non-matching metadata (seem to be rare)
#When set to FALSE, it re-writes mis-matched metadata to match that for the first sample
#verbose=TRUE will read out the mismatches 

seRNA_Mef2c_v13_E85_E9 <- import10xFeatureMatrix(
  input = c("Mef2c_inputs/sample05_filtered_feature_bc_matrix.h5",
            "Mef2c_inputs/sample06_filtered_feature_bc_matrix.h5",
            "Mef2c_inputs/sample07_filtered_feature_bc_matrix.h5",
            "Mef2c_inputs/sample08_filtered_feature_bc_matrix.h5",
            "Mef2c_inputs/sample21_filtered_feature_bc_matrix.h5", 
            "Mef2c_inputs/sample22_filtered_feature_bc_matrix.h5", 
            "Mef2c_inputs/sample23_filtered_feature_bc_matrix.h5", 
            "Mef2c_inputs/sample24_filtered_feature_bc_matrix.h5"
            ),
  names = c("sample05", "sample06", "sample07", "sample08",
            "sample21", "sample22", "sample23", "sample24"
            ),
  strictMatch = FALSE,
  verbose = TRUE,
  featureType = "Gene Expression"
)

#Note: at this point there are 83,373 cells in the ArchR project that passed
#QC based on the ATAC data and  87,122 cells in the seRNA ojbect with gene 
#expression data. We will need to subset ArchR project to contain only cells 
#that have info in both.
overlap <- intersect(proj_Mef2c_v13_E85_E9_ATAC_only$cellNames, seRNA_Mef2c_v13_E85_E9@colData@rownames)
#There are 79,123 cells that exist in both the ArchR project and the seRNA object

#Subset ArchR project to only contain cells that exist in seRNA

#IMPORTANT: Need to install and use "dev_emptyChr" branch 
detach("package:ArchR", unload = TRUE)
devtools::install_github("GreenleafLab/ArchR", ref="dev_emptyChr", repos = BiocManager::repositories())
library(ArchR)

#NOTE: Request either merging dev_empthyChr branch with current release or 
#adding excludeChr 
#UPDATE: I think this has been merged into release_1.0.3, but for now I am
#sticking with using release_1.0.2 + this dev_emptyChr branch at this step

proj_Mef2c_v13_E85_E9 <- subsetArchRProject(
  ArchRProj = proj_Mef2c_v13_E85_E9_ATAC_only,
  cells = overlap,
  outputDirectory = "Mef2c_v13_E85_E9_ATAC_and_GEX",
  dropCells = TRUE,
  force = TRUE
)

#Inspect project - should be 79,123 cells
proj_Mef2c_v13_E85_E9

#IMPORTANT: Re-install "release_1.0.2" branch before proceeding
detach("package:ArchR", unload = TRUE)
devtools::install_github("GreenleafLab/ArchR", ref="release_1.0.2", repos = BiocManager::repositories())
library(ArchR)

#Add the seRNA data to ArchR project
proj_Mef2c_v13_E85_E9 <- addGeneExpressionMatrix(
  input = proj_Mef2c_v13_E85_E9, 
  seRNA = seRNA_Mef2c_v13_E85_E9, 
  excludeChr = c("chrM", "chrY"), 
  strictMatch = FALSE,
  force = TRUE
)


###---------QC and additional TSSEnrichment and nFrag Filtering----------------

#Based on QC metrics, I think I want to further refine cutoff for TSSEnrichment
#Let's build the QC plot TSSEnrichment vs. log10 Unique Frags now that cells have
#been filtered

#Get nFrags and TSSEnrichment scores
df <- getCellColData(proj_Mef2c_v13_E85_E9, select = c("log10(nFrags)", "TSSEnrichment"))
df
#Build plot
p <- ggPoint(
  x = df[,1], 
  y = df[,2], 
  colorDensity = TRUE,
  continuousSet = "sambaNight",
  xlabel = "Log10 Unique Fragments",
  ylabel = "TSS Enrichment",
  xlim = c(log10(500), quantile(df[,1], probs = 1.0)+0.5),
  ylim = c(0, quantile(df[,2], probs = 1.0)+5)
) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")
#Save plot
ggsave("Mef2c_v13_E85_E9_TSSbyLog10nFrags_initial.pdf", plot = p, device = "pdf", path = "~/Mef2c_ArchR_working/QualityControl", width = 5, height = 5, dpi = 300)

#Also want to filter out cells with high/low GEX UMI/Gene counts - let's look at GEX data
df2 <- getCellColData(proj_Mef2c_v13_E85_E9, select = c("Gex_nUMI", "Gex_nGenes"))
#Build plot
p2 <- ggPoint(
  x = df2[,1], 
  y = df2[,2], 
  colorDensity = TRUE,
  continuousSet = "sambaNight",
  xlabel = "Gex_nUMI",
  ylabel = "Gex_nGenes",
)
#Save plot
ggsave("Mef2c_v13_E85_E9_nUMIbynGenes_initial.pdf", plot = p2, device = "pdf", path = "~/Mef2c_ArchR_working/QualityControl", width = 5, height = 5, dpi = 300)

#Before filtering, visualize subset of cells with new cutoffs
#TSS Enrichiment >= 13, TSS Enrichiment <= 30, nFrags >= 4000 (log10 3.6), Gex_nGenes >= 500, Gex_nGenes <= 9000
#Note: using less aggressive cutoffs from E9 analysis for nFrags and Gex filtering
idxPass <- which(proj_Mef2c_v13_E85_E9$TSSEnrichment >= 13 & 
                   proj_Mef2c_v13_E85_E9$TSSEnrichment <= 30 & 
                   proj_Mef2c_v13_E85_E9$nFrags >=4000 & 
                   proj_Mef2c_v13_E85_E9$Gex_nGenes >= 500 &    
                   proj_Mef2c_v13_E85_E9$Gex_nGenes <= 9000)     
cellsPass <- proj_Mef2c_v13_E85_E9$cellNames[idxPass]
df3 <- getCellColData(proj_Mef2c_v13_E85_E9[cellsPass, ], select = c("log10(nFrags)", "TSSEnrichment"))
p3 <- ggPoint(
  x = df3[,1], 
  y = df3[,2], 
  colorDensity = TRUE,
  continuousSet = "sambaNight",
  xlabel = "Log10 Unique Fragments",
  ylabel = "TSS Enrichment",
  xlim = c(log10(500), quantile(df3[,1], probs = 1.0)+0.5),
  ylim = c(0, quantile(df3[,2], probs = 1.0)+5)
) + geom_hline(yintercept = 13, lty = "dashed") + geom_hline(yintercept = 30, lty = "dashed") + geom_vline(xintercept = 3.6, lty = "dashed")
#Save plot
ggsave("Mef2c_v13_E85_E9_TSSbyLog10nFrags_filtered.pdf", plot = p3, device = "pdf", path = "~/Mef2c_ArchR_working/QualityControl", width = 5, height = 5, dpi = 300)

#And the GEX plot
df4 <- getCellColData(proj_Mef2c_v13_E85_E9[cellsPass, ], select = c("Gex_nUMI", "Gex_nGenes"))
#Build plot
p4 <- ggPoint(
  x = df4[,1], 
  y = df4[,2], 
  colorDensity = TRUE,
  continuousSet = "sambaNight",
  xlabel = "Gex_nUMI",
  ylabel = "Gex_nGenes",
)
#Save plot
ggsave("Mef2c_v13_E85_E9_nUMIbynGenes_filtered.pdf", plot = p4, device = "pdf", path = "~/Mef2c_ArchR_working/QualityControl", width = 5, height = 5, dpi = 300)

#Looks good, now filter these cells from project 
proj_Mef2c_v13_E85_E9 <- proj_Mef2c_v13_E85_E9[cellsPass, ]
proj_Mef2c_v13_E85_E9
#Now at 66,597 cells


###----------------Doublet Filtering-----------------------------------------

#First add doublet scores
dir.create("~/Mef2c_ArchR_working/Mef2c_v13_E85_E9_ATAC_and_GEX/DoubletSummaries")
proj_Mef2c_v13_E85_E9 <- addDoubletScores(proj_Mef2c_v13_E85_E9, outDir = paste(getOutputDirectory(proj_Mef2c_v13_E85_E9), "DoubletSummaries", sep = "/"))

#Filter doublets
proj_Mef2c_v13_E85_E9 <- filterDoublets(proj_Mef2c_v13_E85_E9, filterRatio = 1)

#Filtering 6002 cells from ArchRProject!
  #sample05 : 991 of 9955 (10%)
  #sample06 : 408 of 6393 (6.4%)
  #sample07 : 1264 of 11243 (11.2%)
  #sample08 : 1298 of 11394 (11.4%)
  #sample21 : 816 of 9034 (9%)
  #sample22 : 533 of 7301 (7.3%)
  #sample23 : 537 of 7334 (7.3%)
  #sample24 : 155 of 3943 (3.9%)

#Now at  60,595 cells
saveArchRProject(ArchRProj = proj_Mef2c_v13_E85_E9, overwrite = TRUE, dropCells = TRUE, load = TRUE)


###--------------Dim Reduction - LSI ATAC -----------------------------------
proj_Mef2c_v13_E85_E9 <- addIterativeLSI(
  ArchRProj = proj_Mef2c_v13_E85_E9, 
  name = "LSI_ATAC",
  useMatrix = "TileMatrix", 
  depthCol = "nFrags",
  firstSelection = "top",
  binarize = TRUE,
  iterations = 2, 
  varFeatures = 25000,
  dimsToUse = 1:30,
  clusterParams = list(
    resolution = 0.2, 
    sampleCells = 10000,
    n.start = 10
  ),
  saveIterations = FALSE,
  excludeChr = c("chrM", "chrY"),
  force = TRUE, 
  seed = 1
)

###--------------Dim Reduction - LSI RNA -----------------------------------
proj_Mef2c_v13_E85_E9 <- addIterativeLSI(
  ArchRProj = proj_Mef2c_v13_E85_E9, 
  name = "LSI_RNA",
  useMatrix = "GeneExpressionMatrix", 
  depthCol = "Gex_nUMI",
  firstSelection = "var",
  binarize = FALSE,
  iterations = 2,
  varFeatures = 5000,
  dimsToUse = 1:30,
  clusterParams = list(
    resolution = 0.2, 
    sampleCells = 10000,
    n.start = 10
  ),
  saveIterations = FALSE,
  excludeChr = c("chrM", "chrY"), 
  force = TRUE,
  seed = 1
)

###------------Combined Dims-------------------------------------------------
proj_Mef2c_v13_E85_E9 <- addCombinedDims(
  ArchRProj = proj_Mef2c_v13_E85_E9, 
  reducedDims = c("LSI_ATAC", "LSI_RNA"), 
  name =  "LSI_Combined"
)

###------------Harmony Correction for Combined Dims-------------------------
proj_Mef2c_v13_E85_E9 <- addHarmony(
  ArchRProj = proj_Mef2c_v13_E85_E9, 
  reducedDims = "LSI_Combined", 
  name = "Harmony_Combined", 
  groupBy = "Sample", 
  force = T
)

###-------------UMAPs--------------------------------------------------------
#ATAC UMAP
proj_Mef2c_v13_E85_E9 <- addUMAP(
  ArchRProj = proj_Mef2c_v13_E85_E9, 
  reducedDims = "LSI_ATAC", 
  name = "UMAP_ATAC", 
  minDist = 0.5, 
  force = TRUE,
  seed = 1
)

#RNA UMAP
proj_Mef2c_v13_E85_E9 <- addUMAP(
  ArchRProj = proj_Mef2c_v13_E85_E9, 
  reducedDims = "LSI_RNA", 
  name = "UMAP_RNA", 
  minDist = 0.5, 
  force = TRUE,
  seed = 1
)

#Combined UMAP
proj_Mef2c_v13_E85_E9 <- addUMAP(
  ArchRProj = proj_Mef2c_v13_E85_E9, 
  reducedDims = "LSI_Combined", 
  name = "UMAP_Combined", 
  minDist = 0.5, 
  force = TRUE,
  seed = 1
)

#Harmony Combined UMAP
proj_Mef2c_v13_E85_E9 <- addUMAP(
  ArchRProj = proj_Mef2c_v13_E85_E9, 
  reducedDims = "Harmony_Combined", 
  name = "UMAP_Harmony_Combined", 
  minDist = 0.5, 
  force = TRUE,
  seed = 1
)

#Add clusters from Harmony Combined dims
proj_Mef2c_v13_E85_E9 <- addClusters(
  input = proj_Mef2c_v13_E85_E9, 
  reducedDims = "Harmony_Combined", 
  maxClusters = 30,
  name = "Clus_Harm_Comb_res0.5", 
  resolution = 0.5, 
  force = TRUE,
  seed = 1
)

#Plot embeddings
u1 <- plotEmbedding(
  ArchRProj = proj_Mef2c_v13_E85_E9, 
  name = "Clus_Harm_Comb_res0.5", 
  embedding = "UMAP_ATAC", 
  size = 0.5, 
  labelAsFactors=F, 
  labelMeans=F,
  colorTitle = "",
  baseSize = 16,
) + theme(legend.text = element_text(size=16)) + guides(colour = guide_legend(override.aes = list(size=3)))

u2 <- plotEmbedding(
  ArchRProj = proj_Mef2c_v13_E85_E9, 
  name = "Clus_Harm_Comb_res0.5", 
  embedding = "UMAP_RNA", 
  size = 0.5, 
  labelAsFactors=F, 
  labelMeans=F,
  colorTitle = "",
  baseSize = 16,
) + theme(legend.text = element_text(size=16)) + guides(colour = guide_legend(override.aes = list(size=3)))

u3 <- plotEmbedding(
  ArchRProj = proj_Mef2c_v13_E85_E9, 
  name = "Clus_Harm_Comb_res0.5", 
  embedding = "UMAP_Combined", 
  size = 0.5, 
  labelAsFactors=F, 
  labelMeans=F,
  colorTitle = "",
  baseSize = 16,
) + theme(legend.text = element_text(size=16)) + guides(colour = guide_legend(override.aes = list(size=3)))

u4 <- plotEmbedding(
  ArchRProj = proj_Mef2c_v13_E85_E9, 
  name = "Clus_Harm_Comb_res0.5", 
  embedding = "UMAP_Harmony_Combined", 
  size = 0.5, 
  labelAsFactors=F, 
  labelMeans=F,
  colorTitle = "",
  baseSize = 16,
) + theme(legend.text = element_text(size=16)) + guides(colour = guide_legend(override.aes = list(size=3)))

plotPDF(u1, u2, u3, u4, name = "UMAPs_ATAC_RNA_Combined_HarmonyComb_res0.5", width = 10, height = 10, ArchRProj = proj_Mef2c_v13_E85_E9, addDOC = F)

#UMAPs colored by sample
#First, label cells by sample name - samples 22 and 23 are WT, samples 21 and 24 are KO
sn1 <- gsub("sample21", "E8.5_KO1", proj_Mef2c_v13_E85_E9$Sample)
sn2 <- gsub("sample22", "E8.5_WT1", sn1)
sn3 <- gsub("sample23", "E8.5_WT2", sn2)
sn4 <- gsub("sample24", "E8.5_KO2", sn3)
sn5 <- gsub("sample05", "E9.0_WT1", sn4)
sn6 <- gsub("sample06", "E9.0_KO1", sn5)
sn7 <- gsub("sample07", "E9.0_KO2", sn6)
sn8 <- gsub("sample08", "E9.0_WT2", sn7)
proj_Mef2c_v13_E85_E9$SampleNames <- sn8

jon_cols = c("E8.5_KO1" = "#9C27B0",
             "E8.5_WT1" = "#F3E5F5", 
             "E8.5_WT2" = "#D1C4E9", 
             "E8.5_KO2" = "#4527A0",
             "E9.0_WT1" = "pink",
             "E9.0_KO1" = "red", 
             "E9.0_KO2" = "red4", 
             "E9.0_WT2" = "plum1"
             )

u5 <- plotEmbedding(
  proj_Mef2c_v13_E85_E9, 
  name = "SampleNames", 
  embedding = "UMAP_ATAC", 
  size = 0.5, 
  labelAsFactors=F, 
  labelMeans=F, 
  pal = jon_cols,
  colorTitle = ""
) + theme(legend.text = element_text(size=10)) + guides(colour = guide_legend(override.aes = list(size=3)))

u6 <- plotEmbedding(
  proj_Mef2c_v13_E85_E9, 
  name = "SampleNames", 
  embedding = "UMAP_RNA", 
  size = 0.5, 
  labelAsFactors=F, 
  labelMeans=F, 
  pal = jon_cols,
  colorTitle = ""
) + theme(legend.text = element_text(size=10)) + guides(colour = guide_legend(override.aes = list(size=3)))

u7 <- plotEmbedding(
  proj_Mef2c_v13_E85_E9, 
  name = "SampleNames", 
  embedding = "UMAP_Combined", 
  size = 0.5, 
  labelAsFactors=F, 
  labelMeans=F, 
  pal = jon_cols,
  colorTitle = ""
) + theme(legend.text = element_text(size=10)) + guides(colour = guide_legend(override.aes = list(size=3)))

u8 <- plotEmbedding(
  proj_Mef2c_v13_E85_E9, 
  name = "SampleNames", 
  embedding = "UMAP_Harmony_Combined", 
  size = 0.5, 
  labelAsFactors=F, 
  labelMeans=F, 
  pal = jon_cols,
  colorTitle = ""
) + theme(legend.text = element_text(size=10)) + guides(colour = guide_legend(override.aes = list(size=3)))

plotPDF(u5, u6, u7, u8, name = "UMAPs_ATAC_RNA_Combined_HarmonyComb_bySample", ArchRProj = proj_Mef2c_v13_E85_E9, addDOC = F)

#Check doublet scores and other QC metrics
u9 <- plotEmbedding(proj_Mef2c_v13_E85_E9, name = "DoubletScore", embedding = "UMAP_Harmony_Combined", size = 0.5, labelAsFactors=F, labelMeans=F)
u10 <- plotEmbedding(proj_Mef2c_v13_E85_E9, name = "DoubletEnrichment", embedding = "UMAP_Harmony_Combined", size = 0.5, labelAsFactors=F, labelMeans=F)
u11 <- plotEmbedding(proj_Mef2c_v13_E85_E9, name = "nFrags", embedding = "UMAP_Harmony_Combined", size = 0.5, labelAsFactors=F, labelMeans=F)
u12 <- plotEmbedding(proj_Mef2c_v13_E85_E9, name = "TSSEnrichment", embedding = "UMAP_Harmony_Combined", size = 0.5, labelAsFactors=F, labelMeans=F)
u13 <- plotEmbedding(proj_Mef2c_v13_E85_E9, name = "NucleosomeRatio", embedding = "UMAP_Harmony_Combined", size = 0.5, labelAsFactors=F, labelMeans=F)
u14 <- plotEmbedding(proj_Mef2c_v13_E85_E9, name = "Gex_nUMI", embedding = "UMAP_Harmony_Combined", size = 0.5, labelAsFactors=F, labelMeans=F)
u15 <- plotEmbedding(proj_Mef2c_v13_E85_E9, name = "Gex_nGenes", embedding = "UMAP_Harmony_Combined", size = 0.5, labelAsFactors=F, labelMeans=F)
plotPDF(u9, u10, u11, u12, u13, u14, u15, name = "additional_QC_Combined", ArchRProj = proj_Mef2c_v13_E85_E9, addDOC = F)


###--------------Identifying Clusters----------------------------------------

#First check Seurat Labels on ArchR UMAPs

library(Seurat)
#Load in Seurat objects: "mef2c_v13_E85_v3" and "mef2c_v13_E9" 

#Get barcodes from ArchR project
archr.barcodes.list <- as.data.frame(proj_Mef2c_v13_E85_E9$cellNames)[,1]

#Get barcodes and cell type labels from Seurat project
seurat.df.E85 <- as.data.frame(mef2c_v13_E85_v3_harmony$harmony_cell_type) 
names(seurat.df.E85) <- "Label"
seurat.df.E85$Label <- paste("E8.5_", seurat.df.E85$Label, sep = "")
seurat.df.E9 <- as.data.frame(mef2c_v13_E9$cell_type)
names(seurat.df.E9) <- "Label"
seurat.df.E9$Label <- paste("E9_", seurat.df.E9$Label, sep = "")
seurat.df <- rbind(seurat.df.E85, seurat.df.E9)
head(seurat.df)
tail(seurat.df)

#Make Seurat barcodes match ArchR formatting
corrected_barcodes <- gsub("_","#",rownames(seurat.df))
rownames(seurat.df) <- corrected_barcodes
seurat.barcodes.list <- rownames(seurat.df)
#Find matching barcodes
matches <- intersect(archr.barcodes.list, seurat.barcodes.list)
#For matched barcodes, add Seurat cluster labels to ArchR metadata
proj_Mef2c_v13_E85_E9$seurat_labels <- NA
proj_Mef2c_v13_E85_E9@cellColData[matches,"seurat_labels"] <- seurat.df[matches,1]

u50 <- plotEmbedding(
  proj_Mef2c_v13_E85_E9, 
  name = "seurat_labels", 
  embedding = "UMAP_Harmony_Combined", 
  size = 0.5,
  labelAsFactors=F, 
  labelMeans=T,
  labelSize = 3,
  baseSize = 10,
  colorTitle = ""
) + theme(legend.text = element_text(size=10)) + guides(colour = guide_legend(override.aes = list(size=3)))

plotPDF(u50, 
        name = "UMAP_HarmonyComb_E85_E9_seurat_labels", 
        ArchRProj = proj_Mef2c_v13_E85_E9, 
        addDOC = F
)

#Now use Seurat labels to label ArchR clusters ------------------------------

#C1 = Blood
#C2 = Yolk Sac
#C3 = Endothelium
#C4 = CMs1
#C5 = En1
#C6 = En2
#C7 = En3
#C8 = SoM
#C9 = PhM
#C10 = Pe/Me
#C11 = CMs2/SHF
#C12 = MixM
#C13 = SpM
#C14 = C14
#C15 = NT1
#C16 = NT2
#C17 = NT3
#C18 = NT4
#C19 = NT5
#C20 = En4
#C21 = NC

#Highlight a few clusters at a time
idxs <- which(proj_Mef2c_v13_E85_E9$Clus_Harm_Comb_res0.5 == "C20" 
              #  proj_Mef2c_v13_E85_E9$Clus_Harm_Comb_res0.5 == "C11" | 
              #  proj_Mef2c_v13_E85_E9$Clus_Harm_Comb_res0.5 == "C13" |
              #  proj_Mef2c_v13_E85_E9$Clus_Harm_Comb_res0.5 == "C17" 
              )
cells <- proj_Mef2c_v13_E85_E9$cellNames[idxs]
plotEmbedding(
  ArchRProj = proj_Mef2c_v13_E85_E9, 
  name = "Clus_Harm_Comb_res0.5", 
  embedding = "UMAP_Harmony_Combined",
  highlightCells = cells,
  size = 0.5, labelAsFactors=F, labelMeans=F
)

oldlabels <- c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10",
               "C11", "C12", "C13", "C14", "C15", "C16", "C17", "C18", "C19",
               "C20", "C21"
               )

newlabels <- c("Blood", "YS", "Endo", "CMs1", "En1", "En2", "En3",
               "SoM", "PhM", "Pe/Me", "CMs2/SHF", "MixM", "SpM",
               "C14", "NT1", "NT2", "NT3", "NT4", "NT5", "En4", "NC"
               )

proj_Mef2c_v13_E85_E9$Cluster_Labels <- mapLabels(proj_Mef2c_v13_E85_E9$Clus_Harm_Comb_res0.5, newLabels = newlabels, oldLabels = oldlabels)

u21 <- plotEmbedding(
  proj_Mef2c_v13_E85_E9, 
  name = "Cluster_Labels", 
  embedding = "UMAP_Harmony_Combined", 
  size = 0.5, 
  labelAsFactors=F, 
  labelMeans=T,
  labelSize = 3,
  baseSize = 10,
  colorTitle = ""
  ) + theme(legend.text = element_text(size=10)) + guides(colour = guide_legend(override.aes = list(size=3)))

plotPDF(u21, 
        name = "UMAP_HarmonyComb_res0.5_Labeled", 
        ArchRProj = proj_Mef2c_v13_E85_E9, 
        addDOC = F
        )


###--------------Subset--------------------------------------------------

#Subset cells of interest, in this case CMs, SHF, and Pe/ME

subset_ids <- which(proj_Mef2c_v13_E85_E9$Cluster_Labels == "CMs1" | 
                    proj_Mef2c_v13_E85_E9$Cluster_Labels == "CMs2/SHF" |
                    proj_Mef2c_v13_E85_E9$Cluster_Labels == "Pe/Me"  
                    )
subset_cells <- proj_Mef2c_v13_E85_E9$cellNames[subset_ids]
proj_Mef2c_v13_E85_E9_subset <- subsetArchRProject(
                            ArchRProj = proj_Mef2c_v13_E85_E9,
                            cells = subset_cells,
                            outputDirectory = "~/Mef2c_ArchR_working/Mef2c_v13_E85_E9_subset",
                            force = TRUE,
                            dropCells = TRUE
                            )

proj_Mef2c_v13_E85_E9_subset
#Subset project contains 11891 cells


###--------------Subset Dim Reduction - LSI ATAC -----------------------------------
proj_Mef2c_v13_E85_E9_subset <- addIterativeLSI(
  ArchRProj = proj_Mef2c_v13_E85_E9_subset, 
  name = "LSI_ATAC_sub",
  useMatrix = "TileMatrix", 
  depthCol = "nFrags",
  firstSelection = "top",
  binarize = TRUE,
  iterations = 2, 
  varFeatures = 25000,
  dimsToUse = 1:30,
  clusterParams = list(
    resolution = 0.2, 
    sampleCells = 10000,
    n.start = 10
  ),
  saveIterations = FALSE,
  excludeChr = c("chrM", 
                 "chrY"
                 ),
  seed = 1,
  force = TRUE
)

###--------------Subset Dim Reduction - LSI RNA -----------------------------------
proj_Mef2c_v13_E85_E9_subset <- addIterativeLSI(
  ArchRProj = proj_Mef2c_v13_E85_E9_subset, 
  name = "LSI_RNA_sub",
  useMatrix = "GeneExpressionMatrix", 
  depthCol = "Gex_nUMI",
  firstSelection = "var",
  binarize = FALSE,
  iterations = 2,
  varFeatures = 5000,
  dimsToUse = 1:30,
  clusterParams = list(
    resolution = 0.2, 
    sampleCells = 10000,
    n.start = 10
  ),
  saveIterations = FALSE,
  excludeChr = c("chrM", 
                 "chrY"
                 ), 
  seed = 1,
  force = TRUE
)

###------------Subset Combined Dims-------------------------------------------------
proj_Mef2c_v13_E85_E9_subset <- addCombinedDims(
  ArchRProj = proj_Mef2c_v13_E85_E9_subset, 
  reducedDims = c("LSI_ATAC_sub", "LSI_RNA_sub"), 
  name =  "LSI_Combined_sub"
)

###------------Harmony Correction for Combined Dims-------------------------
proj_Mef2c_v13_E85_E9_subset <- addHarmony(
  ArchRProj = proj_Mef2c_v13_E85_E9_subset, 
  reducedDims = "LSI_Combined_sub", 
  name = "Harmony_Combined_sub", 
  groupBy = "Sample", 
  force = T
)

###-------------Subset UMAPs--------------------------------------------------------

#ATAC UMAP
proj_Mef2c_v13_E85_E9_subset <- addUMAP(
  ArchRProj = proj_Mef2c_v13_E85_E9_subset, 
  reducedDims = "LSI_ATAC_sub", 
  name = "UMAP_ATAC_sub", 
  minDist = 0.5,
  seed = 1,
  force = TRUE
)

#RNA UMAP
proj_Mef2c_v13_E85_E9_subset <- addUMAP(
  ArchRProj = proj_Mef2c_v13_E85_E9_subset, 
  reducedDims = "LSI_RNA_sub", 
  name = "UMAP_RNA_sub", 
  minDist = 0.5, 
  seed = 1,
  force = TRUE
)

#Combined UMAP
proj_Mef2c_v13_E85_E9_subset <- addUMAP(
  ArchRProj = proj_Mef2c_v13_E85_E9_subset, 
  reducedDims = "LSI_Combined_sub", 
  name = "UMAP_Combined_sub", 
  minDist = 0.5,
  seed = 1,
  force = TRUE
)

#Harmony Combined UMAP
proj_Mef2c_v13_E85_E9_subset <- addUMAP(
  ArchRProj = proj_Mef2c_v13_E85_E9_subset, 
  reducedDims = "Harmony_Combined_sub", 
  name = "UMAP_Harmony_Combined_sub", 
  minDist = 0.5, 
  seed = 1,
  force = TRUE
)

#Add clusters from Harmony Combined dims
proj_Mef2c_v13_E85_E9_subset <- addClusters(
  input = proj_Mef2c_v13_E85_E9_subset, 
  reducedDims = "Harmony_Combined_sub", 
  maxClusters = 30,
  name = "Subset_Clus_Harm_Comb_res1", 
  resolution = 1, 
  seed = 1,
  force = TRUE
)

#Plot embeddings
u31 <- plotEmbedding(
  ArchRProj = proj_Mef2c_v13_E85_E9_subset, 
  name = "Subset_Clus_Harm_Comb_res1", 
  embedding = "UMAP_ATAC_sub", 
  size = 1, 
  labelAsFactors=F, 
  labelMeans=F,
  colorTitle = "",
  baseSize = 16,
) + theme(legend.text = element_text(size=16)) + guides(colour = guide_legend(override.aes = list(size=3)))

u32 <- plotEmbedding(
  ArchRProj = proj_Mef2c_v13_E85_E9_subset, 
  name = "Subset_Clus_Harm_Comb_res1", 
  embedding = "UMAP_RNA_sub", 
  size = 1, 
  labelAsFactors=F, 
  labelMeans=F,
  colorTitle = "",
  baseSize = 16,
) + theme(legend.text = element_text(size=16)) + guides(colour = guide_legend(override.aes = list(size=3)))

u33 <- plotEmbedding(
  ArchRProj = proj_Mef2c_v13_E85_E9_subset, 
  name = "Subset_Clus_Harm_Comb_res1", 
  embedding = "UMAP_Combined_sub", 
  size = 1, 
  labelAsFactors=F, 
  labelMeans=F,
  colorTitle = "",
  baseSize = 16,
) + theme(legend.text = element_text(size=16)) + guides(colour = guide_legend(override.aes = list(size=3)))

u34 <- plotEmbedding(
  ArchRProj = proj_Mef2c_v13_E85_E9_subset, 
  name = "Subset_Clus_Harm_Comb_res1", 
  embedding = "UMAP_Harmony_Combined_sub", 
  size = 1, 
  labelAsFactors=F, 
  labelMeans=F,
  colorTitle = "",
  baseSize = 16,
) + theme(legend.text = element_text(size=16)) + guides(colour = guide_legend(override.aes = list(size=3)))

plotPDF(u31, u32, u33, u34, name = "UMAPs_ATAC_RNA_Combined_HarmonyComb_E85_E9_Subset_res1", width = 10, height = 10, ArchRProj = proj_Mef2c_v13_E85_E9_subset, addDOC = F)

jon_cols = c("E8.5_KO1" = "#9C27B0",
             "E8.5_WT1" = "#F3E5F5", 
             "E8.5_WT2" = "#D1C4E9", 
             "E8.5_KO2" = "#4527A0",
             "E9.0_WT1" = "pink",
             "E9.0_KO1" = "red", 
             "E9.0_KO2" = "red4", 
             "E9.0_WT2" = "plum1"
)

u35 <- plotEmbedding(
  proj_Mef2c_v13_E85_E9_subset, 
  name = "SampleNames", 
  embedding = "UMAP_ATAC_sub", 
  size = 1, 
  labelAsFactors=F, 
  labelMeans=F, 
  pal = jon_cols,
  colorTitle = ""
) + theme(legend.text = element_text(size=10)) + guides(colour = guide_legend(override.aes = list(size=3)))

u36 <- plotEmbedding(
  proj_Mef2c_v13_E85_E9_subset, 
  name = "SampleNames", 
  embedding = "UMAP_RNA_sub", 
  size = 1, 
  labelAsFactors=F, 
  labelMeans=F, 
  pal = jon_cols,
  colorTitle = ""
) + theme(legend.text = element_text(size=10)) + guides(colour = guide_legend(override.aes = list(size=3)))

u37 <- plotEmbedding(
  proj_Mef2c_v13_E85_E9_subset, 
  name = "SampleNames", 
  embedding = "UMAP_Combined_sub", 
  size = 1, 
  labelAsFactors=F, 
  labelMeans=F, 
  pal = jon_cols,
  colorTitle = ""
) + theme(legend.text = element_text(size=10)) + guides(colour = guide_legend(override.aes = list(size=3)))

u38 <- plotEmbedding(
  proj_Mef2c_v13_E85_E9_subset, 
  name = "SampleNames", 
  embedding = "UMAP_Harmony_Combined_sub", 
  size = 1, 
  labelAsFactors=F, 
  labelMeans=F, 
  pal = jon_cols,
  colorTitle = ""
) + theme(legend.text = element_text(size=10)) + guides(colour = guide_legend(override.aes = list(size=3)))

plotPDF(u35, u36, u37, u38, name = "UMAPs_ATAC_RNA_Combined_HarmonyComb_E85_E9_Subset_bySample", ArchRProj = proj_Mef2c_v13_E85_E9_subset, addDOC = F)


###---------------Identifying Subset Clusters-------------------------------

#First check Seurat Labels on ArchR UMAPs------------------------------------

library(Seurat)
#Load in Seurat objects: "mef2c_v13_E85_v3_subset" and "mef2c_v13_E9_subset"

#Get barcodes from ArchR project
archr.barcodes.list.sub <- as.data.frame(proj_Mef2c_v13_E85_E9_subset$cellNames)[,1]

#Get barcodes and cell type labels from Seurat project
seurat.df.E85.sub <- as.data.frame(mef2c_v13_E85_v3_harmony_subset$harmony_cell_type_subset_pool) 
names(seurat.df.E85.sub) <- "Subset_Label"
seurat.df.E85.sub$Subset_Label <- paste("E8.5_", seurat.df.E85.sub$Subset_Label, sep = "")
seurat.df.E9.sub <- as.data.frame(mef2c_v13_E9_subset$cell_type_subset_pooled)
names(seurat.df.E9.sub) <- "Subset_Label"
seurat.df.E9.sub$Subset_Label <- paste("E9_", seurat.df.E9.sub$Subset_Label, sep = "")
seurat.df.sub <- rbind(seurat.df.E85.sub, seurat.df.E9.sub)
head(seurat.df.sub)
tail(seurat.df.sub)

#Make Seurat barcodes match ArchR formatting
corrected_barcodes_sub <- gsub("_","#",rownames(seurat.df.sub))
rownames(seurat.df.sub) <- corrected_barcodes_sub
seurat.barcodes.list.sub <- rownames(seurat.df.sub)
#Find matching barcodes
matches_sub <- intersect(archr.barcodes.list.sub, seurat.barcodes.list.sub)
#For matched barcodes, add Seurat cluster labels to ArchR metadata
proj_Mef2c_v13_E85_E9_subset$sub_seurat_labels <- NA
proj_Mef2c_v13_E85_E9_subset@cellColData[matches_sub,"sub_seurat_labels"] <- seurat.df.sub[matches_sub,1]

u51 <- plotEmbedding(
  proj_Mef2c_v13_E85_E9_subset, 
  name = "sub_seurat_labels", 
  embedding = "UMAP_Harmony_Combined_sub", 
  size = 0.5,
  labelAsFactors=F, 
  labelMeans=T,
  labelSize = 3,
  baseSize = 10,
  colorTitle = ""
) + theme(legend.text = element_text(size=10)) + guides(colour = guide_legend(override.aes = list(size=3)))

plotPDF(u51, 
        name = "UMAP_HarmonyComb_E85_E9_subset_seurat_labels", 
        ArchRProj = proj_Mef2c_v13_E85_E9_subset, 
        addDOC = F
)

#Now use Seurat Labels to label Subset ArchR Cluster Labels ------------------------------

#C1 = CMs-V1
#C2 = CMs-V2/CMs-A
#C3 = CMs-V3/CMs-AVC
#C4 = CMs-IFT1
#C5 = CMs-IFT2
#C6 = C6
#C7 = CMs-OFT
#C8 = SHF1
#C9 = SHF2
#C10 = SHF3/PostM1
#C11 = SHF4
#C12 = LPM1
#C13 = LPM2
#C14 = LPM3
#C15 = VP1/PostM2
#C16 = Pe
#C17 = VP2/CMs-IFT3

oldlabels2 <- c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10",
               "C11", "C12", "C13", "C14", "C15", "C16", "C17"
               )

newlabels2 <- c("CMs_V1", "CMs_V2/CMs_A", "CMs_V3/CMs_AVC", "CMs_IFT1", 
               "CMs_IFT2", "C6", "CMs_OFT", "SHF1", "SHF2", "SHF3/PostM1",
               "SHF4", "LPM1", "LPM2", "LPM3", "VP1/PostM2", "Pe", "VP2/CMs_IFT3"
               )

proj_Mef2c_v13_E85_E9_subset$Cluster_Labels_Subset <- mapLabels(proj_Mef2c_v13_E85_E9_subset$Subset_Clus_Harm_Comb_res1, newLabels = newlabels2, oldLabels = oldlabels2)

u40 <- plotEmbedding(
  proj_Mef2c_v13_E85_E9_subset, 
  name = "Cluster_Labels_Subset", 
  embedding = "UMAP_Harmony_Combined_sub", 
  size = 1,
  labelAsFactors=F, 
  labelMeans=T,
  labelSize = 3,
  baseSize = 10,
  colorTitle = ""
) + theme(legend.text = element_text(size=10)) + guides(colour = guide_legend(override.aes = list(size=3)))

plotPDF(u40, 
        name = "UMAP_HarmonyComb_E85_E9_Subset_res1_labeled", 
        ArchRProj = proj_Mef2c_v13_E85_E9_subset, 
        addDOC = F
)


###----------Prepare to Create Pseudobulk Replicates--------------------------------

#Label cells by genotype - samples 5, 8, 22, 23 are WT, samples 6, 7, 21, 24 are KO
geno_WT <- gsub("sample22|sample23|sample05|sample08", "WT", proj_Mef2c_v13_E85_E9_subset$Sample)
geno <- gsub("sample21|sample24|sample06|sample07", "KO", geno_WT)
proj_Mef2c_v13_E85_E9_subset$Genotype <- geno

#Label cells by timepoint
time_E85 <- gsub("sample21|sample22|sample23|sample24", "E8.5", proj_Mef2c_v13_E85_E9_subset$Sample)
time <- gsub("sample05|sample06|sample07|sample08", "E9", time_E85)
proj_Mef2c_v13_E85_E9_subset$Timepoint <- time

#Next, create CellTypeByGenotype column for grouping
proj_Mef2c_v13_E85_E9_subset <- addCellColData(ArchRProj = proj_Mef2c_v13_E85_E9_subset, 
                                              data = paste0(proj_Mef2c_v13_E85_E9_subset$Cluster_Labels_Subset,"_x_",
                                                            proj_Mef2c_v13_E85_E9_subset$Genotype), 
                                              name = "CellTypeByGenotype", 
                                              cells = getCellNames(proj_Mef2c_v13_E85_E9_subset), 
                                              force = TRUE
                                              )

#Create CellTypeBySample column to allow easy viewing of number of cells per pseudobulk rep
proj_Mef2c_v13_E85_E9_subset <- addCellColData(ArchRProj = proj_Mef2c_v13_E85_E9_subset, 
                                               data = paste0(proj_Mef2c_v13_E85_E9_subset$Cluster_Labels_Subset,"_x_",
                                                             proj_Mef2c_v13_E85_E9_subset$Sample), 
                                               name = "CellTypeBySample", 
                                               cells = getCellNames(proj_Mef2c_v13_E85_E9_subset), 
                                               force = TRUE
)


###------------Add Group Coverages aka Create Pseudobulk Reps------------------

#Check cell count tables to inform pseudobulking params
table(proj_Mef2c_v13_E85_E9_subset$CellTypeByGenotype)
table(proj_Mef2c_v13_E85_E9_subset$CellTypeBySample)

#Note - because only 17 CM_OFT cells in sample 24, allow min reps = 3 to essentially exclude that sample

#Now addGroupCoverages
proj_Mef2c_v13_E85_E9_subset <- addGroupCoverages(
  ArchRProj = proj_Mef2c_v13_E85_E9_subset,
  groupBy = "CellTypeByGenotype",
  useLabels = TRUE,
  sampleLabels = "Sample",
  minCells = 40,
  maxCells = 700,
  maxFragments = 25 * 10^6,
  minReplicates = 3,
  maxReplicates = 4,
  sampleRatio = 0.8,
  kmerLength = 6,
  threads = getArchRThreads(), 
  returnGroups = FALSE,
  parallelParam = NULL,
  force = TRUE,
  verbose = TRUE,
  logFile = createLogFile("addGroupCoverages")
)


###-------------Add Peak Set and Peak Matrix to Project----------------###

proj_Mef2c_v13_E85_E9_subset <- addReproduciblePeakSet(
  ArchRProj = proj_Mef2c_v13_E85_E9_subset,
  groupBy = "CellTypeByGenotype",
  #Note: with reproducibility = 2, at least 2 pseudobulk rep for each group will
  #be required to contain a peak call at the locus 
  reproducibility = "2",
  peaksPerCell = 500,
  maxPeaks = 150000,
  minCells = 20,
  excludeChr = c("chrM",
                 "chrY"
  ),
  pathToMacs2 = "/Library/Frameworks/Python.framework/Versions/3.6/lib/python3.6/site-packages/MACS2-2.2.7.1-py3.6-macosx-10.9-x86_64.egg/EGG-INFO/scripts/macs2",
  shift = -75,
  extsize = 150,
  cutOff = 0.1,
  extendSummits = 250,
  promoterRegion = c(2000, 100),
  plot = TRUE,
  force = TRUE,
  verbose = TRUE,
  logFile = createLogFile("addReproduciblePeakSet")
)

#Group info used for peak calling

#                                  Group nCells nCellsUsed nReplicates nMin nMax maxPeaks
#C6_x_KO                         C6_x_KO    232        232           3   59  109   116000
#C6_x_WT                         C6_x_WT    238        238           3   70   85   119000
#CMs_IFT1_x_KO             CMs_IFT1_x_KO    493        493           4   70  193   150000
#CMs_IFT1_x_WT             CMs_IFT1_x_WT    456        456           4   64  136   150000
#CMs_IFT2_x_KO             CMs_IFT2_x_KO    387        387           4   55  141   150000
#CMs_IFT2_x_WT             CMs_IFT2_x_WT    178        161           3   49   56    80500
#CMs_OFT_x_KO               CMs_OFT_x_KO    314        297           3   46  127   148500
#CMs_OFT_x_WT               CMs_OFT_x_WT    773        773           4  165  246   150000
#CMs_V1_x_KO                 CMs_V1_x_KO     36         36           3   23   26    18000
#CMs_V1_x_WT                 CMs_V1_x_WT    720        720           4  132  253   150000
#CMs_V2/CMs_A_x_KO     CMs_V2/CMs_A_x_KO    261        261           3   40  146   130500
#CMs_V2/CMs_A_x_WT     CMs_V2/CMs_A_x_WT    575        575           4  119  213   150000
#CMs_V3/CMs_AVC_x_KO CMs_V3/CMs_AVC_x_KO   1015       1015           4  119  405   150000
#CMs_V3/CMs_AVC_x_WT CMs_V3/CMs_AVC_x_WT     30         30           3   18   24    15000
#LPM1_x_KO                     LPM1_x_KO    181        181           3   40   89    90500
#LPM1_x_WT                     LPM1_x_WT    236        236           3   40  112   118000
#LPM2_x_KO                     LPM2_x_KO    231        231           3   52   99   115500
#LPM2_x_WT                     LPM2_x_WT    324        298           3   84  109   149000
#LPM3_x_KO                     LPM3_x_KO    413        413           4   43  149   150000
#LPM3_x_WT                     LPM3_x_WT    338        338           4   63  100   150000
#Pe_x_KO                         Pe_x_KO    325        325           3   40  159   150000
#Pe_x_WT                         Pe_x_WT    155        155           3   40   59    77500
#SHF1_x_KO                     SHF1_x_KO    225        225           3   40  116   112500
#SHF1_x_WT                     SHF1_x_WT    486        457           3   53  257   150000
#SHF2_x_KO                     SHF2_x_KO    267        267           3   42  163   133500
#SHF2_x_WT                     SHF2_x_WT    574        574           4   72  191   150000
#SHF3/PostM1_x_KO       SHF3/PostM1_x_KO    391        376           3   74  192   150000
#SHF3/PostM1_x_WT       SHF3/PostM1_x_WT    406        406           4   54  179   150000
#SHF4_x_KO                     SHF4_x_KO    216        216           3   40  109   108000
#SHF4_x_WT                     SHF4_x_WT      3          3           3    3    3     1500
#VP1/PostM2_x_KO         VP1/PostM2_x_KO    323        284           3   75  126   142000
#VP1/PostM2_x_WT         VP1/PostM2_x_WT    282        282           4   40  116   141000
#VP2/CMs_IFT3_x_KO     VP2/CMs_IFT3_x_KO    487        487           4   68  206   150000
#VP2/CMs_IFT3_x_WT     VP2/CMs_IFT3_x_WT    320        320           3   71  153   150000

#addPeakMatrix
proj_Mef2c_v13_E85_E9_subset <- addPeakMatrix(proj_Mef2c_v13_E85_E9_subset)

#Check that peak matrix appears in available matrices
getAvailableMatrices(proj_Mef2c_v13_E85_E9_subset)


###------------Export Group Bigwigs (optional)-----------------------------------

getGroupBW(
  ArchRProj = proj_Mef2c_v13_E85_E9_subset,
  groupBy = "CellTypeByGenotype",
  normMethod = "ReadsInTSS",
  tileSize = 100,
  maxCells = 1000,
  ceiling = 4,
  verbose = TRUE,
  threads = getArchRThreads(),
  logFile = createLogFile("getGroupBW")
)


###--------------Add Motif Annotations----------------------------------------
#Using new Vierstra motif sets (in Beta) which attempts to de-duplicate similar motifs
#Described here: https://github.com/GreenleafLab/ArchR/discussions/1386
#As of 10/4/2022, Vierstra motifs are version beta2.1: https://resources.altius.org/~jvierstra/projects/motif-clustering-v2.1beta/
proj_Mef2c_v13_E85_E9_subset <- addMotifAnnotations(ArchRProj = proj_Mef2c_v13_E85_E9_subset, 
                                                motifSet = "vierstra", 
                                                collection = "archetype", 
                                                annoName = "Motif", 
                                                force = T
)


###-------------Add Peak2Gene Linkages----------------------------------------
proj_Mef2c_v13_E85_E9_subset <- addPeak2GeneLinks(
  ArchRProj = proj_Mef2c_v13_E85_E9_subset,
  reducedDims = "Harmony_Combined_sub",
  useMatrix = "GeneExpressionMatrix"
)


###--------------Save and Load ArchR Project------------------------------------
saveArchRProject(ArchRProj = proj_Mef2c_v13_E85_E9, overwrite = TRUE, dropCells = FALSE, load = TRUE)
saveArchRProject(ArchRProj = proj_Mef2c_v13_E85_E9_subset, overwrite = TRUE, dropCells = FALSE, load = TRUE)

#to load
proj_Mef2c_v13_E85_E9 <- loadArchRProject("~/Mef2c_ArchR_working/Mef2c_v13_E85_E9_ATAC_and_GEX")
proj_Mef2c_v13_E85_E9_subset <- loadArchRProject("~/Mef2c_ArchR_working/Mef2c_v13_E85_E9_subset")
setwd("~/Mef2c_ArchR_working")

