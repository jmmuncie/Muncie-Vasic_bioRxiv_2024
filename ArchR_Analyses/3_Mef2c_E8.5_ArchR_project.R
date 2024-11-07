###----------------Description-------------------------------------

#This script is used to perform integrated analyses of snRNA-seq and snATAC-seq
#data in ArchR for the four E8.5 Mef2c Multiome samples 

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
proj_Mef2c_v13_E85_ATAC_only <- ArchRProject(
  ArrowFiles = c("Original_Arrows/sample21.arrow", 
                 "Original_Arrows/sample22.arrow", 
                 "Original_Arrows/sample23.arrow", 
                 "Original_Arrows/sample24.arrow"
                 ),
  outputDirectory = "Mef2c_E85_ATAC_only",
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

seRNA_Mef2c_v13_E85 <- import10xFeatureMatrix(
  input = c("Mef2c_inputs/sample21_filtered_feature_bc_matrix.h5", 
            "Mef2c_inputs/sample22_filtered_feature_bc_matrix.h5", 
            "Mef2c_inputs/sample23_filtered_feature_bc_matrix.h5", 
            "Mef2c_inputs/sample24_filtered_feature_bc_matrix.h5"
            ),
  names = c("sample21", "sample22", "sample23", "sample24"),
  strictMatch = FALSE,
  verbose = TRUE,
  featureType = "Gene Expression"
)

#Note: at this point there are 32,868 cells in the ArchR project that passed
#QC based on the ATAC data and  34,354 cells in the seRNA ojbect with gene 
#expression data. We will need to subset ArchR project to contain only cells 
#that have info in both.
overlap <- intersect(proj_Mef2c_v13_E85_ATAC_only$cellNames, seRNA_Mef2c_v13_E85@colData@rownames)
#There are 31,278 cells that exist in both the ArchR project and the seRNA object

#Add the seRNA data to ArchR project
#Note: will get warning that not all cells in input exist in seRNA - that's 
#okay, we will subset the overlapping cells next
proj_Mef2c_v13_E85_ATAC_only <- addGeneExpressionMatrix(
  input = proj_Mef2c_v13_E85_ATAC_only, 
  seRNA = seRNA_Mef2c_v13_E85, 
  excludeChr = c("chrM", "chrY"), 
  strictMatch = FALSE,
  force = TRUE
)

#Subset ArchR project to only contain cells that exist in seRNA

#IMPORTANT: Need to install and use "dev_emptyChr" branch 
detach("package:ArchR", unload = TRUE)
devtools::install_github("GreenleafLab/ArchR", ref="dev_emptyChr", repos = BiocManager::repositories())
library(ArchR)

  #NOTE: Request either merging dev_empthyChr branch with current release or 
  #adding excludeChr 
  #UPDATE: I think this has been merged into release_1.0.3, but for now I am
  #sticking with using release_1.0.2 + this dev_emptyChr branch at this step

proj_Mef2c_v13_E85 <- subsetArchRProject(
  ArchRProj = proj_Mef2c_v13_E85_ATAC_only,
  cells = overlap,
  outputDirectory = "Mef2c_v13_E85_ATAC_and_GEX",
  dropCells = TRUE,
  force = TRUE
)

#Inspect project - should be 31,278 cells
proj_Mef2c_v13_E85

#IMPORTANT: Re-install "release_1.0.2" branch before proceeding
devtools::install_github("GreenleafLab/ArchR", ref="release_1.0.2", repos = BiocManager::repositories())
detach("package:ArchR", unload = TRUE)
library(ArchR)


###---------QC and additional TSSEnrichment and nFrag Filtering----------------

#Based on QC metrics, I think I want to further refine cutoff for TSSEnrichment
#Let's build the QC plot TSSEnrichment vs. log10 Unique Frags now that cells have
#been filtered

#Get nFrags and TSSEnrichment scores
df <- getCellColData(proj_Mef2c_v13_E85, select = c("log10(nFrags)", "TSSEnrichment"))
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
ggsave("Mef2c_v13_E85_TSSbyLog10nFrags_initial.pdf", plot = p, device = "pdf", path = "~/Mef2c_ArchR_working/QualityControl", width = 5, height = 5, dpi = 300)

#Also want to filter out cells with high/low GEX UMI/Gene counts - let's look at GEX data
df2 <- getCellColData(proj_Mef2c_v13_E85, select = c("Gex_nUMI", "Gex_nGenes"))
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
ggsave("Mef2c_v13_E85_nUMIbynGenes_initial.pdf", plot = p2, device = "pdf", path = "~/Mef2c_ArchR_working/QualityControl", width = 5, height = 5, dpi = 300)

#Before filtering, visualize subset of cells with new cutoffs
#TSS Enrichiment >= 13, TSS Enrichiment <= 30 nFrags >= 6300 (log10 3.8), Gex_nGenes >= 1750, Gex_nGenes <= 5500
idxPass <- which(proj_Mef2c_v13_E85$TSSEnrichment >= 13 & 
                   proj_Mef2c_v13_E85$TSSEnrichment <= 30 & 
                   proj_Mef2c_v13_E85$nFrags >=6300 & 
                   proj_Mef2c_v13_E85$Gex_nGenes >= 1750 &    #Same cutoffs as Seurat analysis
                   proj_Mef2c_v13_E85$Gex_nGenes <= 5500)     #Same cutoffs as Seurat analysis
cellsPass <- proj_Mef2c_v13_E85$cellNames[idxPass]
df3 <- getCellColData(proj_Mef2c_v13_E85[cellsPass, ], select = c("log10(nFrags)", "TSSEnrichment"))
p3 <- ggPoint(
  x = df3[,1], 
  y = df3[,2], 
  colorDensity = TRUE,
  continuousSet = "sambaNight",
  xlabel = "Log10 Unique Fragments",
  ylabel = "TSS Enrichment",
  xlim = c(log10(500), quantile(df3[,1], probs = 1.0)+0.5),
  ylim = c(0, quantile(df3[,2], probs = 1.0)+5)
) + geom_hline(yintercept = 13, lty = "dashed") + geom_hline(yintercept = 30, lty = "dashed") + geom_vline(xintercept = 3.8, lty = "dashed")
#Save plot
ggsave("Mef2c_v13_E85_TSSbyLog10nFrags_filtered.pdf", plot = p3, device = "pdf", path = "~/Mef2c_ArchR_working/QualityControl", width = 5, height = 5, dpi = 300)

#And the GEX plot
df4 <- getCellColData(proj_Mef2c_v13_E85[cellsPass, ], select = c("Gex_nUMI", "Gex_nGenes"))
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
ggsave("Mef2c_v13_E85_nUMIbynGenes_filtered.pdf", plot = p4, device = "pdf", path = "~/Mef2c_ArchR_working/QualityControl", width = 5, height = 5, dpi = 300)

#Looks good, now filter these cells from project 
proj_Mef2c_v13_E85 <- proj_Mef2c_v13_E85[cellsPass, ]
proj_Mef2c_v13_E85
#Now at 23,913 cells


###----------------Doublet Filtering-----------------------------------------

#First add doublet scores
dir.create("~/Mef2c_ArchR_working/Mef2c_v13_E85_ATAC_and_GEX/DoubletSummaries")
proj_Mef2c_v13_E85 <- addDoubletScores(proj_Mef2c_v13_E85, outDir = paste(getOutputDirectory(proj_Mef2c_v13_E85), "DoubletSummaries", sep = "/"))

#Filter doublets
proj_Mef2c_v13_E85 <- filterDoublets(proj_Mef2c_v13_E85, filterRatio = 1)

#Filtering 1534 cells from ArchRProject!
  #sample21 : 571 of 7558 (7.6%)
  #sample22 : 384 of 6200 (6.2%)
  #sample23 : 472 of 6875 (6.9%)
  #sample24 : 107 of 3280 (3.3%)

#Now at  22,379 cells
saveArchRProject(ArchRProj = proj_Mef2c_v13_E85, overwrite = TRUE, dropCells = TRUE, load = TRUE)


###--------------Dim Reduction - LSI ATAC -----------------------------------
proj_Mef2c_v13_E85 <- addIterativeLSI(
  ArchRProj = proj_Mef2c_v13_E85, 
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
proj_Mef2c_v13_E85 <- addIterativeLSI(
  ArchRProj = proj_Mef2c_v13_E85, 
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
proj_Mef2c_v13_E85 <- addCombinedDims(
  ArchRProj = proj_Mef2c_v13_E85, 
  reducedDims = c("LSI_ATAC", "LSI_RNA"), 
  name =  "LSI_Combined"
)

###------------Harmony Correction for Combined Dims-------------------------
proj_Mef2c_v13_E85 <- addHarmony(
  ArchRProj = proj_Mef2c_v13_E85, 
  reducedDims = "LSI_Combined", 
  name = "Harmony_Combined", 
  groupBy = "Sample", 
  force = T
)

###-------------UMAPs--------------------------------------------------------
#ATAC UMAP
proj_Mef2c_v13_E85 <- addUMAP(
  ArchRProj = proj_Mef2c_v13_E85, 
  reducedDims = "LSI_ATAC", 
  name = "UMAP_ATAC", 
  minDist = 0.5, 
  force = TRUE,
  seed = 1
)

#RNA UMAP
proj_Mef2c_v13_E85 <- addUMAP(
  ArchRProj = proj_Mef2c_v13_E85, 
  reducedDims = "LSI_RNA", 
  name = "UMAP_RNA", 
  minDist = 0.5, 
  force = TRUE,
  seed = 1
)

#Combined UMAP
proj_Mef2c_v13_E85 <- addUMAP(
  ArchRProj = proj_Mef2c_v13_E85, 
  reducedDims = "LSI_Combined", 
  name = "UMAP_Combined", 
  minDist = 0.5, 
  force = TRUE,
  seed = 1
)

#Harmony Combined UMAP
proj_Mef2c_v13_E85 <- addUMAP(
  ArchRProj = proj_Mef2c_v13_E85, 
  reducedDims = "Harmony_Combined", 
  name = "UMAP_Harmony_Combined", 
  minDist = 0.5, 
  force = TRUE,
  seed = 1
)

#Add clusters from Harmony Combined dims
proj_Mef2c_v13_E85 <- addClusters(
  input = proj_Mef2c_v13_E85, 
  reducedDims = "Harmony_Combined", 
  maxClusters = 30,
  name = "Clus_Harm_Comb_res0.5", 
  resolution = 0.5, 
  force = TRUE,
  seed = 1
)

#Plot embeddings
u1 <- plotEmbedding(
  ArchRProj = proj_Mef2c_v13_E85, 
  name = "Clus_Harm_Comb_res0.5", 
  embedding = "UMAP_ATAC", 
  size = 0.5, 
  labelAsFactors=F, 
  labelMeans=F,
  colorTitle = "",
  baseSize = 16,
) + theme(legend.text = element_text(size=16)) + guides(colour = guide_legend(override.aes = list(size=3)))

u2 <- plotEmbedding(
  ArchRProj = proj_Mef2c_v13_E85, 
  name = "Clus_Harm_Comb_res0.5", 
  embedding = "UMAP_RNA", 
  size = 0.5, 
  labelAsFactors=F, 
  labelMeans=F,
  colorTitle = "",
  baseSize = 16,
) + theme(legend.text = element_text(size=16)) + guides(colour = guide_legend(override.aes = list(size=3)))

u3 <- plotEmbedding(
  ArchRProj = proj_Mef2c_v13_E85, 
  name = "Clus_Harm_Comb_res0.5", 
  embedding = "UMAP_Combined", 
  size = 0.5, 
  labelAsFactors=F, 
  labelMeans=F,
  colorTitle = "",
  baseSize = 16,
) + theme(legend.text = element_text(size=16)) + guides(colour = guide_legend(override.aes = list(size=3)))

u4 <- plotEmbedding(
  ArchRProj = proj_Mef2c_v13_E85, 
  name = "Clus_Harm_Comb_res0.5", 
  embedding = "UMAP_Harmony_Combined", 
  size = 0.5, 
  labelAsFactors=F, 
  labelMeans=F,
  colorTitle = "",
  baseSize = 16,
) + theme(legend.text = element_text(size=16)) + guides(colour = guide_legend(override.aes = list(size=3)))

plotPDF(u1, u2, u3, u4, name = "UMAPs_ATAC_RNA_Combined_HarmonyComb_res0.5", width = 10, height = 10, ArchRProj = proj_Mef2c_v13_E85, addDOC = F)

#UMAPs colored by sample
#First, label cells by sample name - samples 22 and 23 are WT, samples 21 and 24 are KO
sn1 <- gsub("sample21", "E8.5_KO1", proj_Mef2c_v13_E85$Sample)
sn2 <- gsub("sample22", "E8.5_WT1", sn1)
sn3 <- gsub("sample23", "E8.5_WT2", sn2)
sn4 <- gsub("sample24", "E8.5_KO2", sn3)
proj_Mef2c_v13_E85$SampleNames <- sn4

jon_cols = c("E8.5_KO1" = "#9C27B0",
             "E8.5_WT1" = "#F3E5F5", 
             "E8.5_WT2" = "#D1C4E9", 
             "E8.5_KO2" = "#4527A0"
             )

u5 <- plotEmbedding(
  proj_Mef2c_v13_E85, 
  name = "SampleNames", 
  embedding = "UMAP_ATAC", 
  size = 0.5, 
  labelAsFactors=F, 
  labelMeans=F, 
  pal = jon_cols,
  colorTitle = ""
) + theme(legend.text = element_text(size=10)) + guides(colour = guide_legend(override.aes = list(size=3)))

u6 <- plotEmbedding(
  proj_Mef2c_v13_E85, 
  name = "SampleNames", 
  embedding = "UMAP_RNA", 
  size = 0.5, 
  labelAsFactors=F, 
  labelMeans=F, 
  pal = jon_cols,
  colorTitle = ""
) + theme(legend.text = element_text(size=10)) + guides(colour = guide_legend(override.aes = list(size=3)))

u7 <- plotEmbedding(
  proj_Mef2c_v13_E85, 
  name = "SampleNames", 
  embedding = "UMAP_Combined", 
  size = 0.5, 
  labelAsFactors=F, 
  labelMeans=F, 
  pal = jon_cols,
  colorTitle = ""
) + theme(legend.text = element_text(size=10)) + guides(colour = guide_legend(override.aes = list(size=3)))

u8 <- plotEmbedding(
  proj_Mef2c_v13_E85, 
  name = "SampleNames", 
  embedding = "UMAP_Harmony_Combined", 
  size = 0.5, 
  labelAsFactors=F, 
  labelMeans=F, 
  pal = jon_cols,
  colorTitle = ""
) + theme(legend.text = element_text(size=10)) + guides(colour = guide_legend(override.aes = list(size=3)))

plotPDF(u5, u6, u7, u8, name = "UMAPs_ATAC_RNA_Combined_HarmonyComb_bySample", ArchRProj = proj_Mef2c_v13_E85, addDOC = F)

#Check doublet scores and other QC metrics

#QC Plots
u9 <- plotEmbedding(proj_Mef2c_v13_E85, name = "DoubletScore", embedding = "UMAP_Harmony_Combined", size = 0.5, labelAsFactors=F, labelMeans=F)
u10 <- plotEmbedding(proj_Mef2c_v13_E85, name = "DoubletEnrichment", embedding = "UMAP_Harmony_Combined", size = 0.5, labelAsFactors=F, labelMeans=F)
u11 <- plotEmbedding(proj_Mef2c_v13_E85, name = "nFrags", embedding = "UMAP_Harmony_Combined", size = 0.5, labelAsFactors=F, labelMeans=F)
u12 <- plotEmbedding(proj_Mef2c_v13_E85, name = "TSSEnrichment", embedding = "UMAP_Harmony_Combined", size = 0.5, labelAsFactors=F, labelMeans=F)
u13 <- plotEmbedding(proj_Mef2c_v13_E85, name = "NucleosomeRatio", embedding = "UMAP_Harmony_Combined", size = 0.5, labelAsFactors=F, labelMeans=F)
u14 <- plotEmbedding(proj_Mef2c_v13_E85, name = "Gex_nUMI", embedding = "UMAP_Harmony_Combined", size = 0.5, labelAsFactors=F, labelMeans=F)
u15 <- plotEmbedding(proj_Mef2c_v13_E85, name = "Gex_nGenes", embedding = "UMAP_Harmony_Combined", size = 0.5, labelAsFactors=F, labelMeans=F)
plotPDF(u9, u10, u11, u12, u13, u14, u15, name = "additional_QC_Combined", ArchRProj = proj_Mef2c_v13_E85, addDOC = F)


###--------------Identifying Clusters----------------------------------------

#First check Seurat Labels on ArchR UMAPs------------------------------------

library(Seurat)
#Load in Seurat object, in this case "mef2c_v13_E85_v3"

#Get barcodes from ArchR project
archr.barcodes.list <- as.data.frame(proj_Mef2c_v13_E85$cellNames)[,1]
#Get barcodes and cell type labels from Seurat project
seurat.df <- as.data.frame(mef2c_v13_E85_v3_harmony$harmony_cell_type) 
#Make Seurat barcodes match ArchR formatting
corrected_barcodes <- gsub("_","#",rownames(seurat.df))
rownames(seurat.df) <- corrected_barcodes
seurat.barcodes.list <- rownames(seurat.df)
#Find matching barcodes
matches <- intersect(archr.barcodes.list, seurat.barcodes.list)
#For matched barcodes, add Seurat cluster labels to ArchR metadata
proj_Mef2c_v13_E85$seurat_labels <- NA
proj_Mef2c_v13_E85@cellColData[matches,"seurat_labels"] <- seurat.df[matches,1]

u50 <- plotEmbedding(
  proj_Mef2c_v13_E85, 
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
        name = "UMAP_HarmonyComb_E85_seurat_labels", 
        ArchRProj = proj_Mef2c_v13_E85, 
        addDOC = F
)

#Then look at ArchR feature plots on ArchR UMAPs-----------------------------
markerGenes <- c(
  "Ttn", "Tnnt2", "Myl7",                  #CMs/CPs
  "Tbx5", "Mef2c",
  "Isl1", "Fgf8", "Fgf10",                 #SHF
  "Tbx1", "Foxp2", "Ebf1",                 #PhM
  "Bmp4", "Hand1", "Dlk1",                 #LPM
  "Meox1", "Fst", "Foxd1",                 #SoM
  "Pbx1", "Pdgfra", "Dach1",               #MixM
  "Nrcam", "Fgf14", "Pax6",                #Neural Tube
  "Rmst", "Nrp2", "Rfx4",         
  "Sox10", "Tfap2b", "Ets1",               #Neural Crest
  "Epcam", "Cldn6", "Bcam", "Wnt6",        #Endoderm
  "Flt1", "Kdr", "Cdh5", "Nfatc1",         #Endothelium
  "T", "Shh",                              #Notochord
  "Ttr", "Afp", "Cubn"                     #YS
)

FP <- plotEmbedding(
  ArchRProj = proj_Mef2c_v13_E85,
  colorBy = "GeneExpressionMatrix",
  name = markerGenes,
  embedding = "UMAP_Harmony_Combined"
)

Fplot1 <- lapply(FP[1:5], function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})
FP1 <- do.call(cowplot::plot_grid, c(list(ncol = 3),Fplot1))

Fplot2 <- lapply(FP[6:8], function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})
FP2 <- do.call(cowplot::plot_grid, c(list(ncol = 2),Fplot2))

Fplot3 <- lapply(FP[9:11], function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})
FP3 <- do.call(cowplot::plot_grid, c(list(ncol = 2),Fplot3))

Fplot4 <- lapply(FP[12:14], function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})
FP4 <- do.call(cowplot::plot_grid, c(list(ncol = 2),Fplot4))

Fplot5 <- lapply(FP[15:17], function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})
FP5 <- do.call(cowplot::plot_grid, c(list(ncol = 2),Fplot5))

Fplot6 <- lapply(FP[18:20], function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})
FP6 <- do.call(cowplot::plot_grid, c(list(ncol = 2),Fplot6))

Fplot7 <- lapply(FP[21:26], function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})
FP7 <- do.call(cowplot::plot_grid, c(list(ncol = 2),Fplot7))

Fplot8 <- lapply(FP[27:29], function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})
FP8 <- do.call(cowplot::plot_grid, c(list(ncol = 2),Fplot8))

Fplot9 <- lapply(FP[30:33], function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})
FP9 <- do.call(cowplot::plot_grid, c(list(ncol = 2),Fplot9))

Fplot10 <- lapply(FP[34:37], function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})
FP10 <- do.call(cowplot::plot_grid, c(list(ncol = 2),Fplot10))

Fplot11 <- lapply(FP[38:39], function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})
FP11 <- do.call(cowplot::plot_grid, c(list(ncol = 2),Fplot11))

Fplot12 <- lapply(FP[40:42], function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})
FP12 <- do.call(cowplot::plot_grid, c(list(ncol = 2),Fplot12))

plotPDF(FP1, FP2, FP3, FP4, FP5, FP6, FP7, FP8, FP9, FP10, FP11, FP12,  
        name = "FeaturePlots_E85_full_HarmonyComb_UMAP", 
        ArchRProj = proj_Mef2c_v13_E85, 
        addDOC = F
)

#Now label E8.5 - Full Dataset Cluster Labels ------------------------------

#C1 = Yolk Sac
#C2 = Neural Crest
#C3 = Neural Tube 1
#C4 = Neural Tube 2
#C5 = Neural Tube 3
#C6 = Neural Tube 4
#C7 = Neural Tube 5
#C8 = Endothelium
#C9 = Somitic Mesoderm
#C10 = SHF1/PhM
#C11 = Notochord
#C12 = En1
#C13 = En2
#C14 = CPs1
#C15 = CMs/CPs2
#C16 = SHF2/MixM
#C17 = LPM
#C18 = SHF3

oldlabels <- c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10",
               "C11", "C12", "C13", "C14", "C15", "C16", "C17", "C18"
               )
newlabels <- c("YS", "NC", "NT1", "NT2", "NT3", "NT4", "NT5",
               "Endo", "SoM", "SHF1/PhM", "Noto", "En1", "En2",
               "CPs1", "CMs/CPs2", "SHF2/MixM", "LPM", "SHF3"
               )

proj_Mef2c_v13_E85$Cluster_Labels <- mapLabels(proj_Mef2c_v13_E85$Clus_Harm_Comb_res0.5, newLabels = newlabels, oldLabels = oldlabels)

u21 <- plotEmbedding(
  proj_Mef2c_v13_E85, 
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
        ArchRProj = proj_Mef2c_v13_E85, 
        addDOC = F
        )


###--------------Subset E8.5--------------------------------------------------
#Subset clusters of interest, namely CPs/CMs, and clusters containing SHF

subset_ids <- which(proj_Mef2c_v13_E85$Cluster_Labels == "CPs1" | 
                    proj_Mef2c_v13_E85$Cluster_Labels == "CMs/CPs2" |
                    proj_Mef2c_v13_E85$Cluster_Labels == "SHF1/PhM" |
                    proj_Mef2c_v13_E85$Cluster_Labels == "SHF2/MixM" |
                    proj_Mef2c_v13_E85$Cluster_Labels == "SHF3" 
                    )
subset_cells <- proj_Mef2c_v13_E85$cellNames[subset_ids]
proj_Mef2c_v13_E85_subset <- subsetArchRProject(
                            ArchRProj = proj_Mef2c_v13_E85,
                            cells = subset_cells,
                            outputDirectory = "~/Mef2c_ArchR_working/Mef2c_v13_E85_subset",
                            force = TRUE,
                            dropCells = TRUE
                            )

proj_Mef2c_v13_E85_subset
#Subset project contains 5988 cells


###--------------Subset Dim Reduction - LSI ATAC -----------------------------------
proj_Mef2c_v13_E85_subset <- addIterativeLSI(
  ArchRProj = proj_Mef2c_v13_E85_subset, 
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
proj_Mef2c_v13_E85_subset <- addIterativeLSI(
  ArchRProj = proj_Mef2c_v13_E85_subset, 
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
proj_Mef2c_v13_E85_subset <- addCombinedDims(
  ArchRProj = proj_Mef2c_v13_E85_subset, 
  reducedDims = c("LSI_ATAC_sub", "LSI_RNA_sub"), 
  name =  "LSI_Combined_sub"
)

###------------Harmony Correction for Combined Dims-------------------------
proj_Mef2c_v13_E85_subset <- addHarmony(
  ArchRProj = proj_Mef2c_v13_E85_subset, 
  reducedDims = "LSI_Combined_sub", 
  name = "Harmony_Combined_sub", 
  groupBy = "Sample", 
  force = T
)

###-------------Subset UMAPs--------------------------------------------------------

#ATAC UMAP
proj_Mef2c_v13_E85_subset <- addUMAP(
  ArchRProj = proj_Mef2c_v13_E85_subset, 
  reducedDims = "LSI_ATAC_sub", 
  name = "UMAP_ATAC_sub", 
  minDist = 0.5,
  seed = 1,
  force = TRUE
)

#RNA UMAP
proj_Mef2c_v13_E85_subset <- addUMAP(
  ArchRProj = proj_Mef2c_v13_E85_subset, 
  reducedDims = "LSI_RNA_sub", 
  name = "UMAP_RNA_sub", 
  minDist = 0.5, 
  seed = 1,
  force = TRUE
)

#Combined UMAP
proj_Mef2c_v13_E85_subset <- addUMAP(
  ArchRProj = proj_Mef2c_v13_E85_subset, 
  reducedDims = "LSI_Combined_sub", 
  name = "UMAP_Combined_sub", 
  minDist = 0.5,
  seed = 1,
  force = TRUE
)

#Harmony Combined UMAP
proj_Mef2c_v13_E85_subset <- addUMAP(
  ArchRProj = proj_Mef2c_v13_E85_subset, 
  reducedDims = "Harmony_Combined_sub", 
  name = "UMAP_Harmony_Combined_sub", 
  minDist = 0.5, 
  seed = 1,
  force = TRUE
)

#Add clusters from Harmony Combined dims
proj_Mef2c_v13_E85_subset <- addClusters(
  input = proj_Mef2c_v13_E85_subset, 
  reducedDims = "Harmony_Combined_sub", 
  maxClusters = 30,
  name = "Subset_Clus_Harm_Comb_res1", 
  resolution = 1, 
  seed = 1,
  force = TRUE
)

#Plot embeddings
u31 <- plotEmbedding(
  ArchRProj = proj_Mef2c_v13_E85_subset, 
  name = "Subset_Clus_Harm_Comb_res1", 
  embedding = "UMAP_ATAC_sub", 
  size = 1, 
  labelAsFactors=F, 
  labelMeans=F,
  colorTitle = "",
  baseSize = 16,
) + theme(legend.text = element_text(size=16)) + guides(colour = guide_legend(override.aes = list(size=3)))

u32 <- plotEmbedding(
  ArchRProj = proj_Mef2c_v13_E85_subset, 
  name = "Subset_Clus_Harm_Comb_res1", 
  embedding = "UMAP_RNA_sub", 
  size = 1, 
  labelAsFactors=F, 
  labelMeans=F,
  colorTitle = "",
  baseSize = 16,
) + theme(legend.text = element_text(size=16)) + guides(colour = guide_legend(override.aes = list(size=3)))

u33 <- plotEmbedding(
  ArchRProj = proj_Mef2c_v13_E85_subset, 
  name = "Subset_Clus_Harm_Comb_res1", 
  embedding = "UMAP_Combined_sub", 
  size = 1, 
  labelAsFactors=F, 
  labelMeans=F,
  colorTitle = "",
  baseSize = 16,
) + theme(legend.text = element_text(size=16)) + guides(colour = guide_legend(override.aes = list(size=3)))

u34 <- plotEmbedding(
  ArchRProj = proj_Mef2c_v13_E85_subset, 
  name = "Subset_Clus_Harm_Comb_res1", 
  embedding = "UMAP_Harmony_Combined_sub", 
  size = 1, 
  labelAsFactors=F, 
  labelMeans=F,
  colorTitle = "",
  baseSize = 16,
) + theme(legend.text = element_text(size=16)) + guides(colour = guide_legend(override.aes = list(size=3)))

plotPDF(u31, u32, u33, u34, name = "UMAPs_ATAC_RNA_Combined_HarmonyComb_E85_Subset_res1", width = 10, height = 10, ArchRProj = proj_Mef2c_v13_E85_subset, addDOC = F)

jon_cols = c("E8.5_KO1" = "#9C27B0",
             "E8.5_WT1" = "#F3E5F5", 
             "E8.5_WT2" = "#D1C4E9", 
             "E8.5_KO2" = "#4527A0"
             )

u35 <- plotEmbedding(
  proj_Mef2c_v13_E85_subset, 
  name = "SampleNames", 
  embedding = "UMAP_ATAC_sub", 
  size = 1, 
  labelAsFactors=F, 
  labelMeans=F, 
  pal = jon_cols,
  colorTitle = ""
) + theme(legend.text = element_text(size=10)) + guides(colour = guide_legend(override.aes = list(size=3)))

u36 <- plotEmbedding(
  proj_Mef2c_v13_E85_subset, 
  name = "SampleNames", 
  embedding = "UMAP_RNA_sub", 
  size = 1, 
  labelAsFactors=F, 
  labelMeans=F, 
  pal = jon_cols,
  colorTitle = ""
) + theme(legend.text = element_text(size=10)) + guides(colour = guide_legend(override.aes = list(size=3)))

u37 <- plotEmbedding(
  proj_Mef2c_v13_E85_subset, 
  name = "SampleNames", 
  embedding = "UMAP_Combined_sub", 
  size = 1, 
  labelAsFactors=F, 
  labelMeans=F, 
  pal = jon_cols,
  colorTitle = ""
) + theme(legend.text = element_text(size=10)) + guides(colour = guide_legend(override.aes = list(size=3)))

u38 <- plotEmbedding(
  proj_Mef2c_v13_E85_subset, 
  name = "SampleNames", 
  embedding = "UMAP_Harmony_Combined_sub", 
  size = 1, 
  labelAsFactors=F, 
  labelMeans=F, 
  pal = jon_cols,
  colorTitle = ""
) + theme(legend.text = element_text(size=10)) + guides(colour = guide_legend(override.aes = list(size=3)))

plotPDF(u35, u36, u37, u38, name = "UMAPs_ATAC_RNA_Combined_HarmonyComb_E85_Subset_bySample", ArchRProj = proj_Mef2c_v13_E85_subset, addDOC = F)


###---------------Identifying Subset Clusters-------------------------------

#First check Seurat Labels on ArchR UMAPs------------------------------------

library(Seurat)
#Load in Seurat object, in this case "mef2c_v13_E85_v3_subset"

#Get barcodes from ArchR project
archr.barcodes.list <- as.data.frame(proj_Mef2c_v13_E85_subset$cellNames)[,1]
#Get barcodes and cell type labels from Seurat project
seurat.df <- as.data.frame(mef2c_v13_E85_v3_harmony_subset$harmony_cell_type_subset) 
#Make Seurat barcodes match ArchR formatting
corrected_barcodes <- gsub("_","#",rownames(seurat.df))
rownames(seurat.df) <- corrected_barcodes
seurat.barcodes.list <- rownames(seurat.df)
#Find matching barcodes
matches <- intersect(archr.barcodes.list, seurat.barcodes.list)
#For matched barcodes, add Seurat cluster labels to ArchR metadata
proj_Mef2c_v13_E85_subset$seurat_labels_subset <- NA
proj_Mef2c_v13_E85_subset@cellColData[matches,"seurat_labels_subset"] <- seurat.df[matches,1]

u51 <- plotEmbedding(
  proj_Mef2c_v13_E85_subset, 
  name = "seurat_labels_subset", 
  embedding = "UMAP_Harmony_Combined_sub", 
  size = 1,
  labelAsFactors=F, 
  labelMeans=T,
  labelSize = 3,
  baseSize = 10,
  colorTitle = ""
) + theme(legend.text = element_text(size=10)) + guides(colour = guide_legend(override.aes = list(size=3)))

plotPDF(u51, 
        name = "UMAP_HarmonyComb_E85_subset_seurat_labels", 
        ArchRProj = proj_Mef2c_v13_E85_subset, 
        addDOC = F
)

#Then look at ArchR feature plots on ArchR UMAPs-----------------------------
markerGenes_sub <- c(
  "Ttn", "Nkx2-5", "Myl2", "Irx4",              #CMs-V
  "Tdgf1", "Rgs5", "Nrp2",                      #CMs-OFT
  "Tbx5", "Wnt2", "Mef2c", "Gata4", "Tbx20",    #CMs-IFT
  "Isl1", "Fgf8", "Fgf10",                      #aSHF1
  "Nrg1", "Irx5", "Lef1", "Bmp5",               #aSHF2/3
  "Isl1", "Aldh1a2", "Foxp2", "Foxf1",          #pSHF
  "Pdgfra", "Vegfc", "Epha7",                   #MixM
  "Nrip1", "Arg1",                              
  "Fst", "Ebf1", "Tbx1",                        #PhM
  "Cdh11", "Foxd1"                             
)

FPs <- plotEmbedding(
  ArchRProj = proj_Mef2c_v13_E85_subset,
  colorBy = "GeneExpressionMatrix",
  name = markerGenes_sub,
  embedding = "UMAP_Harmony_Combined_sub",
  plotAs = "points",
  size = 1
)

Fplot1s <- lapply(FPs[1:4], function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})
FP1s <- do.call(cowplot::plot_grid, c(list(ncol = 2),Fplot1s))

Fplot2s <- lapply(FPs[5:7], function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})
FP2s <- do.call(cowplot::plot_grid, c(list(ncol = 2),Fplot2s))

Fplot3s <- lapply(FPs[8:12], function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})
FP3s <- do.call(cowplot::plot_grid, c(list(ncol = 2),Fplot3s))

Fplot4s <- lapply(FPs[13:15], function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})
FP4s <- do.call(cowplot::plot_grid, c(list(ncol = 2),Fplot4s))

Fplot5s <- lapply(FPs[16:19], function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})
FP5s <- do.call(cowplot::plot_grid, c(list(ncol = 2),Fplot5s))

Fplot6s <- lapply(FPs[20:23], function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})
FP6s <- do.call(cowplot::plot_grid, c(list(ncol = 2),Fplot6s))

Fplot7s <- lapply(FPs[24:26], function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})
FP7s <- do.call(cowplot::plot_grid, c(list(ncol = 2),Fplot7s))

Fplot8s <- lapply(FPs[27:28], function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})
FP8s <- do.call(cowplot::plot_grid, c(list(ncol = 2),Fplot8s))

Fplot9s <- lapply(FPs[29:31], function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})
FP9s <- do.call(cowplot::plot_grid, c(list(ncol = 2),Fplot9s))

Fplot10s <- lapply(FPs[32:33], function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})
FP10s <- do.call(cowplot::plot_grid, c(list(ncol = 2),Fplot10s))

plotPDF(FP1s, FP2s, FP3s, FP4s, FP5s, FP6s, FP7s, FP8s, FP9s, FP10s,
        name = "FeaturePlots_E85_subset_HarmonyComb_UMAP", 
        ArchRProj = proj_Mef2c_v13_E85_subset, 
        addDOC = F
)

#Now label E8.5 - Subset Cluster Labels ------------------------------

#C1 = CMs-KO
#C2 = CMs-IFT
#C3 = CMs-OFT
#C4 = CMs-V
#C5 = PhM1 
#C6 = PhM2
#C7 = PhM3
#C8 = C8
#C9 = PostM
#C10 = aSHF
#C11 = LPM1
#C12 = LPM2/pSHF

oldlabels2 <- c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10",
               "C11", "C12"
               )
newlabels2 <- c("CMs_KO", "CMs_IFT", "CMs_OFT", "CMs_V", 
               "PhM1", "PhM2", "PhM3", "C8",
               "PostM", "aSHF", "LPM1", "LPM2/pSHF"  
      
               )

proj_Mef2c_v13_E85_subset$Cluster_Labels_Subset <- mapLabels(proj_Mef2c_v13_E85_subset$Subset_Clus_Harm_Comb_res1, newLabels = newlabels2, oldLabels = oldlabels2)

u40 <- plotEmbedding(
  proj_Mef2c_v13_E85_subset, 
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
        name = "UMAP_HarmonyComb_E85_Subset_res1_labeled", 
        ArchRProj = proj_Mef2c_v13_E85_subset, 
        addDOC = F
)


###----------Prepare to Create Pseudobulk Replicates--------------------------------

#First, create identity labels with all cell type clusters grouped together\
#Grouped cells together based on Seurat labels

oldlabels3 <- c("aSHF", "C16", "CMs-IFT1", "CMs-IFT2", 
                "CMs-IFT3", "CMs-OFT", "CMs-V", "LPM1",
                "LPM2", "MixM1", "MixM2", "MixM3",  
                "PhM1", "PhM2", "PhM3", "PostM", "pSHF"
)

newlabels3 <- c("aSHF", "C16", "CMs_IFT", "CMs_IFT", 
                "CMs_IFT", "CMs_OFT", "CMs_V", "LPM",
                "LPM", "LPM", "LPM", "LPM",  
                "PhM", "PhM", "PhM", "PostM", "pSHF"
)

proj_Mef2c_v13_E85_subset$Cluster_Labels_Subset_celltypegroup <- mapLabels(proj_Mef2c_v13_E85_subset$seurat_labels_subset, newLabels = newlabels3, oldLabels = oldlabels3)

#Custom palette for figure panel
pal = c("pSHF" = "#FEE500", "PhM" = "#C06CAB", "CMs_IFT" = "#20A39E", "LPM" = "#D8A767", "NA" = "#999999", "PostM" = "#208A42", "C16" = "#89288F", "CMs_OFT" = "#D1495B", "aSHF" = "#90D5E4", "CMs_V" = "#EDAE49")

u41 <- plotEmbedding(
  proj_Mef2c_v13_E85_subset, 
  name = "Cluster_Labels_Subset_celltypegroup", 
  embedding = "UMAP_Harmony_Combined_sub", 
  pal = pal,
  size = 1,
  labelAsFactors=F, 
  labelMeans=T,
  labelSize = 4,
  baseSize = 10,
  colorTitle = ""
) + theme(legend.text = element_text(size=10)) + guides(colour = guide_legend(override.aes = list(size=3)))

plotPDF(u41, 
        name = "UMAP_HarmonyComb_E85_Subset_res1_Seurat-labeled_celltypegroup", 
        ArchRProj = proj_Mef2c_v13_E85_subset, 
        addDOC = F
)

#Next, label cells by genotype - samples 22 and 23 are WT, samples 21 and 24 are KO
geno_WT <- gsub("sample22|sample23", "WT", proj_Mef2c_v13_E85_subset$Sample)
geno <- gsub("sample21|sample24", "KO", geno_WT)
proj_Mef2c_v13_E85_subset$Genotype <- geno

#Next, create CellTypeByGenotype column for grouping
proj_Mef2c_v13_E85_subset <- addCellColData(ArchRProj = proj_Mef2c_v13_E85_subset, 
                                              data = paste0(proj_Mef2c_v13_E85_subset$Cluster_Labels_Subset_celltypegroup,"_x_",
                                                            proj_Mef2c_v13_E85_subset$Genotype), 
                                              name = "CellTypeByGenotype", 
                                              cells = getCellNames(proj_Mef2c_v13_E85_subset), 
                                              force = TRUE
                                              )

#Create CellTypeBySample column to allow easy viewing of number of cells per cluster
proj_Mef2c_v13_E85_subset <- addCellColData(ArchRProj = proj_Mef2c_v13_E85_subset, 
                                              data = paste0(proj_Mef2c_v13_E85_subset$Cluster_Labels_Subset_celltypegroup,"_x_",
                                                            proj_Mef2c_v13_E85_subset$Sample), 
                                              name = "CellTypeBySample", 
                                              cells = getCellNames(proj_Mef2c_v13_E85_subset), 
                                              force = TRUE
                                              )

###------------Add Group Coverages aka Create Pseudobulk Reps------------------

#Check cell count tables to inform pseudobulking params
table(proj_Mef2c_v13_E85_subset$CellTypeByGenotype)
table(proj_Mef2c_v13_E85_subset$CellTypeBySample)

#Now addGroupCoverages
proj_Mef2c_v13_E85_subset <- addGroupCoverages(
  ArchRProj = proj_Mef2c_v13_E85_subset,
  groupBy = "CellTypeByGenotype",
  useLabels = TRUE,
  sampleLabels = "Sample",
  minCells = 20,
  maxCells = 500,
  maxFragments = 25 * 10^6,
  minReplicates = 2,
  maxReplicates = 2,
  sampleRatio = 0.8,
  kmerLength = 6,
  threads = getArchRThreads(), 
  returnGroups = FALSE,
  parallelParam = NULL,
  force = TRUE,
  verbose = TRUE,
  logFile = createLogFile("addGroupCoverages")
)

###-------------------------Add Peak Set and Peak Matrix------------------

#addReproduciblePeakSet
proj_Mef2c_v13_E85_subset <- addReproduciblePeakSet(
  ArchRProj = proj_Mef2c_v13_E85_subset,
  groupBy = "CellTypeByGenotype",
  #Note: with reproducibility = 2, both my psuedobulk reps for each group will
  #be required to contain a peak call at the locus - this will yield strict peak calls 
  reproducibility = "2",
  peaksPerCell = 500,
  maxPeaks = 150000,
  minCells = 20,
  excludeChr = c("chrM",
                 "chrY"
  ),
  pathToMacs2 =  "/opt/homebrew/bin/macs2",
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

#Group info used for peak calls - from addGroupCoverages

#                    Group nCells nCellsUsed nReplicates nMin nMax maxPeaks
#aSHF_x_KO       aSHF_x_KO    106        106           2   20   86    53000
#aSHF_x_WT       aSHF_x_WT    203        203           2   95  108   101500
#C16_x_KO         C16_x_KO      7          7           2    4    7     3500
#C16_x_WT         C16_x_WT      6          6           2    4    6     3000
#CMs_IFT_x_KO CMs_IFT_x_KO    526        526           2  151  375   150000
#CMs_IFT_x_WT CMs_IFT_x_WT    334        334           2  140  194   150000
#CMs_OFT_x_KO CMs_OFT_x_KO    107        107           2   20   87    53500
#CMs_OFT_x_WT CMs_OFT_x_WT    315        315           2  156  159   150000
#CMs_V_x_KO     CMs_V_x_KO    223        223           2   54  169   111500
#CMs_V_x_WT     CMs_V_x_WT    519        519           2  238  281   150000
#LPM_x_KO         LPM_x_KO    480        480           2  137  343   150000
#LPM_x_WT         LPM_x_WT    738        738           2  306  432   150000
#NA_x_KO           NA_x_KO     79         79           2   20   59    39500
#NA_x_WT           NA_x_WT     62         62           2   21   41    31000
#PhM_x_KO         PhM_x_KO    572        572           2  153  419   150000
#PhM_x_WT         PhM_x_WT    751        751           2  371  380   150000
#PostM_x_KO     PostM_x_KO    101        101           2   20   81    50500
#PostM_x_WT     PostM_x_WT    143        143           2   71   72    71500
#pSHF_x_KO       pSHF_x_KO    309        309           2   61  248   150000
#pSHF_x_WT       pSHF_x_WT    407        407           2  184  223   150000

#addPeakMatrix
proj_Mef2c_v13_E85_subset <- addPeakMatrix(proj_Mef2c_v13_E85_subset)

#Check that PeakMatrix now appears in available matrices
getAvailableMatrices(proj_Mef2c_v13_E85_subset)


###------------Export Group Bigwigs (optional)-----------------------------------

getGroupBW(
  ArchRProj = proj_Mef2c_v13_E85_subset,
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
proj_Mef2c_v13_E85_subset <- addMotifAnnotations(ArchRProj = proj_Mef2c_v13_E85_subset, 
                                                 motifSet = "vierstra", 
                                                 collection = "archetype", 
                                                 annoName = "Motif", 
                                                 force = T
                                                 )


###-------------Add Peak2Gene Linkages----------------------------------------
proj_Mef2c_v13_E85_subset <- addPeak2GeneLinks(
  ArchRProj = proj_Mef2c_v13_E85_subset,
  reducedDims = "Harmony_Combined_sub",
  useMatrix = "GeneExpressionMatrix"
)


###--------------Save and Load ArchR Project------------------------------------
saveArchRProject(ArchRProj = proj_Mef2c_v13_E85, overwrite = TRUE, dropCells = FALSE, load = TRUE)
saveArchRProject(ArchRProj = proj_Mef2c_v13_E85_subset, overwrite = TRUE, dropCells = FALSE, load = TRUE)

#to load
proj_Mef2c_v13_E85 <- loadArchRProject("~/Mef2c_ArchR_working/Mef2c_v13_E85_ATAC_and_GEX")
proj_Mef2c_v13_E85_subset <- loadArchRProject("~/Mef2c_ArchR_working/Mef2c_v13_E85_subset")
setwd("~/Mef2c_ArchR_working")



