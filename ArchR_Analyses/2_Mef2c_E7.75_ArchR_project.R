###----------------Description-------------------------------------

#This script is used to perform integrated analyses of snRNA-seq and snATAC-seq
#data in ArchR for the four E7.75 Mef2c Multiome samples 

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
#These are created in sript "1_Mef2c_ArchR_create_arrows.R"
load("~/Genome_Gene_Annotations/genomeAnnotation_mm10v13.Robj")
load("~/Genome_Gene_Annotations/geneAnnotation_mm10v13.Robj")
setwd("~/Mef2c_ArchR_working")


###--------------Create ArchR Project----------------------
#Arrow Files previously created in sript "1_Mef2c_ArchR_create_arrows.R"
#Move/copy arrow files to "~/Mef2c_ArchR_working/Original_Arrows"

#Create ArchR Project
proj_Mef2c_v13_E775_ATAC_only <- ArchRProject(
  ArrowFiles = c("Original_Arrows/sample01.arrow", 
                 "Original_Arrows/sample02.arrow", 
                 "Original_Arrows/sample03.arrow", 
                 "Original_Arrows/sample04.arrow"
                 ),
  outputDirectory = "Mef2c_v13_E775_ATAC_only",
  copyArrows = FALSE,
  geneAnnotation = geneAnnotation_mm10v13,
  genomeAnnotation = genomeAnnotation_mm10v13,
  showLogo = FALSE,
  threads = 4
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

seRNA_Mef2c_v13_E775 <- import10xFeatureMatrix(
  input = c("Mef2c_inputs/sample01_filtered_feature_bc_matrix.h5", 
            "Mef2c_inputs/sample02_filtered_feature_bc_matrix.h5", 
            "Mef2c_inputs/sample03_filtered_feature_bc_matrix.h5", 
            "Mef2c_inputs/sample04_filtered_feature_bc_matrix.h5"
            ),
  names = c("sample01", "sample02", "sample03", "sample04"),
  strictMatch = FALSE,
  verbose = TRUE,
  featureType = "Gene Expression"
)

#Note: at this point there are 30,376 cells in the ArchR project that passed
#QC based on the ATAC data and  28,304 cells in the seRNA ojbect with gene 
#expression data. We will need to subset ArchR project to contain only cells 
#that have info in both.
overlap <- intersect(proj_Mef2c_v13_E775_ATAC_only$cellNames, seRNA_Mef2c_v13_E775@colData@rownames)
#There are 25,792 cells that exist in both the ArchR project and the seRNA object

#Add the seRNA data to ArchR project
#Note: will get warning that not all cells in input exist in seRNA - that's 
#okay, we will subset the overlapping cells next
proj_Mef2c_v13_E775_ATAC_only <- addGeneExpressionMatrix(
  input = proj_Mef2c_v13_E775_ATAC_only, 
  seRNA = seRNA_Mef2c_v13_E775, 
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

proj_Mef2c_v13_E775 <- subsetArchRProject(
  ArchRProj = proj_Mef2c_v13_E775_ATAC_only,
  cells = overlap,
  outputDirectory = "Mef2c_v13_E775_ATAC_and_GEX",
  dropCells = TRUE,
  force = TRUE
)

#Inspect project - should be 25,792 cells
proj_Mef2c_v13_E775

#IMPORTANT: Re-install "release_1.0.2" branch before proceeding
detach("package:ArchR", unload = TRUE)
devtools::install_github("GreenleafLab/ArchR", ref="release_1.0.2", repos = BiocManager::repositories())
library(ArchR)


###---------QC and additional TSSEnrichment and nFrag Filtering----------------

#Based on QC metrics, I think I want to further refine cutoff for TSSEnrichment
#Let's build the QC plot TSSEnrichment vs. log10 Unique Frags now that cells have
#been filtered

#Get nFrags and TSSEnrichment scores
df <- getCellColData(proj_Mef2c_v13_E775, select = c("log10(nFrags)", "TSSEnrichment"))
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
ggsave("Mef2c_E775_TSSbyLog10nFrags_initial.pdf", plot = p, device = "pdf", path = "~/Mef2c_ArchR_working/QualityControl", width = 5, height = 5, dpi = 300)

#Also want to filter out cells with high/low GEX UMI/Gene counts - let's look at GEX data
df2 <- getCellColData(proj_Mef2c_v13_E775, select = c("Gex_nUMI", "Gex_nGenes"))
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
ggsave("Mef2c_v13_E775_nUMIbynGenes_initial.pdf", plot = p2, device = "pdf", path = "~/Mef2c_ArchR_working/QualityControl", width = 5, height = 5, dpi = 300)

#Before filtering, visualize subset of cells with new cutoffs
#TSS Enrichiment >= 12, TSS Enrichiment <= 30 nFrags >= 3160 (log10 3.5), Gex_nGenes >= 2000, Gex_nGenes <= 8000
idxPass <- which(proj_Mef2c_v13_E775$TSSEnrichment >= 12 & 
                   proj_Mef2c_v13_E775$TSSEnrichment <= 30 & 
                   proj_Mef2c_v13_E775$nFrags >=3160 & 
                   proj_Mef2c_v13_E775$Gex_nGenes >= 2000 &    #Same cutoffs as Seurat analysis
                   proj_Mef2c_v13_E775$Gex_nGenes <= 8000)     #Same cutoffs as Seurat analysis
cellsPass <- proj_Mef2c_v13_E775$cellNames[idxPass]
df3 <- getCellColData(proj_Mef2c_v13_E775[cellsPass, ], select = c("log10(nFrags)", "TSSEnrichment"))
p3 <- ggPoint(
  x = df3[,1], 
  y = df3[,2], 
  colorDensity = TRUE,
  continuousSet = "sambaNight",
  xlabel = "Log10 Unique Fragments",
  ylabel = "TSS Enrichment",
  xlim = c(log10(500), quantile(df3[,1], probs = 1.0)+0.5),
  ylim = c(0, quantile(df3[,2], probs = 1.0)+5)
) + geom_hline(yintercept = 12, lty = "dashed") + geom_hline(yintercept = 30, lty = "dashed") + geom_vline(xintercept = 3.5, lty = "dashed")
#Save plot
ggsave("Mef2c_v13_E775_TSSbyLog10nFrags_filtered.pdf", plot = p3, device = "pdf", path = "~/Mef2c_ArchR_working/QualityControl", width = 5, height = 5, dpi = 300)

#And the GEX plot
df4 <- getCellColData(proj_Mef2c_v13_E775[cellsPass, ], select = c("Gex_nUMI", "Gex_nGenes"))
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
ggsave("Mef2c_v13_E775_nUMIbynGenes_filtered.pdf", plot = p4, device = "pdf", path = "~/Mef2c_ArchR_working/QualityControl", width = 5, height = 5, dpi = 300)

#Looks good, now filter these cells from project 
proj_Mef2c_v13_E775 <- proj_Mef2c_v13_E775[cellsPass, ]
proj_Mef2c_v13_E775
#Now at 19,071 cells


###----------------Doublet Filtering-----------------------------------------

#First add doublet scores
dir.create("~/Mef2c_ArchR_working/Mef2c_v13_E775_ATAC_and_GEX/DoubletSummaries")
proj_Mef2c_v13_E775 <- addDoubletScores(proj_Mef2c_v13_E775, outDir = paste(getOutputDirectory(proj_Mef2c_v13_E775), "DoubletSummaries", sep = "/"))

#Filter doublets
proj_Mef2c_v13_E775 <- filterDoublets(proj_Mef2c_v13_E775, filterRatio = 2)

#Filtering 1973 cells from ArchRProject!
  #sample01 : 593 of 5448 (10.9%)
  #sample02 : 905 of 6730 (13.4%)
  #sample03 : 212 of 3260 (6.5%)
  #sample04 : 263 of 3633 (7.2%)

#Now at  17,098 cells
saveArchRProject(ArchRProj = proj_Mef2c_v13_E775, overwrite = TRUE, dropCells = TRUE, load = TRUE)


###--------------Dim Reduction - LSI ATAC -----------------------------------
proj_Mef2c_v13_E775 <- addIterativeLSI(
  ArchRProj = proj_Mef2c_v13_E775, 
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
  force = TRUE
)

###--------------Dim Reduction - LSI RNA -----------------------------------
proj_Mef2c_v13_E775 <- addIterativeLSI(
  ArchRProj = proj_Mef2c_v13_E775, 
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
  force = TRUE
)

###------------Combined Dims-------------------------------------------------
proj_Mef2c_v13_E775 <- addCombinedDims(
  ArchRProj = proj_Mef2c_v13_E775, 
  reducedDims = c("LSI_ATAC", "LSI_RNA"), 
  name =  "LSI_Combined"
)

###------------Harmony Correction for Combined Dims-------------------------
proj_Mef2c_v13_E775 <- addHarmony(
  ArchRProj = proj_Mef2c_v13_E775, 
  reducedDims = "LSI_Combined", 
  name = "Harmony_Combined", 
  groupBy = "Sample", 
  force = T
)

###-------------UMAPs--------------------------------------------------------
#ATAC UMAP
proj_Mef2c_v13_E775 <- addUMAP(
  ArchRProj = proj_Mef2c_v13_E775, 
  reducedDims = "LSI_ATAC", 
  name = "UMAP_ATAC", 
  minDist = 0.5, 
  force = TRUE
)

#RNA UMAP
proj_Mef2c_v13_E775 <- addUMAP(
  ArchRProj = proj_Mef2c_v13_E775, 
  reducedDims = "LSI_RNA", 
  name = "UMAP_RNA", 
  minDist = 0.5, 
  force = TRUE
)

#Combined UMAP
proj_Mef2c_v13_E775 <- addUMAP(
  ArchRProj = proj_Mef2c_v13_E775, 
  reducedDims = "LSI_Combined", 
  name = "UMAP_Combined", 
  minDist = 0.5, 
  force = TRUE
)

#Harmony Combined UMAP
proj_Mef2c_v13_E775 <- addUMAP(
  ArchRProj = proj_Mef2c_v13_E775, 
  reducedDims = "Harmony_Combined", 
  name = "UMAP_Harmony_Combined", 
  minDist = 0.5, 
  force = TRUE
)

#Add clusters from Harmony Combined dims
proj_Mef2c_v13_E775 <- addClusters(
  input = proj_Mef2c_v13_E775, 
  reducedDims = "Harmony_Combined", 
  maxClusters = 30,
  name = "Clus_Harm_Comb_res0.5", 
  resolution = 0.5, 
  force = TRUE
)

#Plot embeddings
u1 <- plotEmbedding(
  ArchRProj = proj_Mef2c_v13_E775, 
  name = "Clus_Harm_Comb_res0.5", 
  embedding = "UMAP_ATAC", 
  size = 0.5, 
  labelAsFactors=F, 
  labelMeans=F,
  colorTitle = "",
  baseSize = 16,
) + theme(legend.text = element_text(size=16)) + guides(colour = guide_legend(override.aes = list(size=3)))

u2 <- plotEmbedding(
  ArchRProj = proj_Mef2c_v13_E775, 
  name = "Clus_Harm_Comb_res0.5", 
  embedding = "UMAP_RNA", 
  size = 0.5, 
  labelAsFactors=F, 
  labelMeans=F,
  colorTitle = "",
  baseSize = 16,
) + theme(legend.text = element_text(size=16)) + guides(colour = guide_legend(override.aes = list(size=3)))

u3 <- plotEmbedding(
  ArchRProj = proj_Mef2c_v13_E775, 
  name = "Clus_Harm_Comb_res0.5", 
  embedding = "UMAP_Combined", 
  size = 0.5, 
  labelAsFactors=F, 
  labelMeans=F,
  colorTitle = "",
  baseSize = 16,
) + theme(legend.text = element_text(size=16)) + guides(colour = guide_legend(override.aes = list(size=3)))

u4 <- plotEmbedding(
  ArchRProj = proj_Mef2c_v13_E775, 
  name = "Clus_Harm_Comb_res0.5", 
  embedding = "UMAP_Harmony_Combined", 
  size = 0.5, 
  labelAsFactors=F, 
  labelMeans=F,
  colorTitle = "",
  baseSize = 16,
) + theme(legend.text = element_text(size=16)) + guides(colour = guide_legend(override.aes = list(size=3)))

plotPDF(u1, u2, u3, u4, name = "UMAPs_ATAC_RNA_Combined_HarmonyComb_Mef2c_E775_res0.5", width = 10, height = 10, ArchRProj = proj_Mef2c_v13_E775, addDOC = F)

#UMAPs colored by sample
#First, label cells by sample name - samples 01 and 04 are WT, samples 02 and 03 are KO
sn1 <- gsub("sample01", "E7.75_WT1", proj_Mef2c_v13_E775$Sample)
sn2 <- gsub("sample02", "E7.75_KO1", sn1)
sn3 <- gsub("sample03", "E7.75_KO2", sn2)
sn4 <- gsub("sample04", "E7.75_WT2", sn3)
proj_Mef2c_v13_E775$SampleNames <- sn4

jon_cols = c("E7.75_WT1" = "lightskyblue",
             "E7.75_KO1" = "blue",
             "E7.75_KO2" = "darkblue",
             "E7.75_WT2" = "turquoise1"
             )

u5 <- plotEmbedding(
  proj_Mef2c_v13_E775, 
  name = "SampleNames", 
  embedding = "UMAP_ATAC", 
  size = 0.5, 
  labelAsFactors=F, 
  labelMeans=F, 
  pal = jon_cols,
  colorTitle = ""
) + theme(legend.text = element_text(size=10)) + guides(colour = guide_legend(override.aes = list(size=3)))

u6 <- plotEmbedding(
  proj_Mef2c_v13_E775, 
  name = "SampleNames", 
  embedding = "UMAP_RNA", 
  size = 0.5, 
  labelAsFactors=F, 
  labelMeans=F, 
  pal = jon_cols,
  colorTitle = ""
) + theme(legend.text = element_text(size=10)) + guides(colour = guide_legend(override.aes = list(size=3)))

u7 <- plotEmbedding(
  proj_Mef2c_v13_E775, 
  name = "SampleNames", 
  embedding = "UMAP_Combined", 
  size = 0.5, 
  labelAsFactors=F, 
  labelMeans=F, 
  pal = jon_cols,
  colorTitle = ""
) + theme(legend.text = element_text(size=10)) + guides(colour = guide_legend(override.aes = list(size=3)))

u8 <- plotEmbedding(
  proj_Mef2c_v13_E775, 
  name = "SampleNames", 
  embedding = "UMAP_Harmony_Combined", 
  size = 0.5, 
  labelAsFactors=F, 
  labelMeans=F, 
  pal = jon_cols,
  colorTitle = ""
) + theme(legend.text = element_text(size=10)) + guides(colour = guide_legend(override.aes = list(size=3)))

plotPDF(u5, u6, u7, u8, name = "UMAPs_ATAC_RNA_Combined_HarmonyComb_Mef2c_E775_bySample", ArchRProj = proj_Mef2c_v13_E775, addDOC = F)

#Check doublet scores and other QC metrics

#QC Plots
u9 <- plotEmbedding(proj_Mef2c_v13_E775, name = "DoubletScore", embedding = "UMAP_Harmony_Combined", size = 0.5, labelAsFactors=F, labelMeans=F)
u10 <- plotEmbedding(proj_Mef2c_v13_E775, name = "DoubletEnrichment", embedding = "UMAP_Harmony_Combined", size = 0.5, labelAsFactors=F, labelMeans=F)
u11 <- plotEmbedding(proj_Mef2c_v13_E775, name = "nFrags", embedding = "UMAP_Harmony_Combined", size = 0.5, labelAsFactors=F, labelMeans=F)
u12 <- plotEmbedding(proj_Mef2c_v13_E775, name = "TSSEnrichment", embedding = "UMAP_Harmony_Combined", size = 0.5, labelAsFactors=F, labelMeans=F)
u13 <- plotEmbedding(proj_Mef2c_v13_E775, name = "NucleosomeRatio", embedding = "UMAP_Harmony_Combined", size = 0.5, labelAsFactors=F, labelMeans=F)
u14 <- plotEmbedding(proj_Mef2c_v13_E775, name = "Gex_nUMI", embedding = "UMAP_Harmony_Combined", size = 0.5, labelAsFactors=F, labelMeans=F)
u15 <- plotEmbedding(proj_Mef2c_v13_E775, name = "Gex_nGenes", embedding = "UMAP_Harmony_Combined", size = 0.5, labelAsFactors=F, labelMeans=F)
plotPDF(u9, u10, u11, u12, u13, u14, u15, name = "additional_QC_Combined_Mef2c_E775", ArchRProj = proj_Mef2c_v13_E775, addDOC = F)


###--------------Identifying Clusters----------------------------------------

#First check Seurat Labels on ArchR UMAPs------------------------------------

library(Seurat)
#Load in Seurat object, in this case "mef2c_v13_E775.Robj"

#Get barcodes from ArchR project
archr.barcodes.list <- as.data.frame(proj_Mef2c_v13_E775$cellNames)[,1]
#Get barcodes and cell type labels from Seurat project
seurat.df <- as.data.frame(mef2c_v13_E775$cell_type) 
#Make Seurat barcodes match ArchR formatting
corrected_barcodes <- gsub("_","#",rownames(seurat.df))
rownames(seurat.df) <- corrected_barcodes
seurat.barcodes.list <- rownames(seurat.df)
#Find matching barcodes
matches <- intersect(archr.barcodes.list, seurat.barcodes.list)
#For matched barcodes, add Seurat cluster labels to ArchR metadata
proj_Mef2c_v13_E775$seurat_labels <- NA
proj_Mef2c_v13_E775@cellColData[matches,"seurat_labels"] <- seurat.df[matches,1]

u50 <- plotEmbedding(
  proj_Mef2c_v13_E775, 
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
        name = "UMAP_HarmonyComb_E775_seurat_labels", 
        ArchRProj = proj_Mef2c_v13_E775, 
        addDOC = F
)

#Then look at ArchR feature plots on ArchR UMAPs-----------------------------
markerGenes <- c(
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
  "Bmper", "Grm8"                       #C16 (from Seurat)
)

FP <- plotEmbedding(
  ArchRProj = proj_Mef2c_v13_E775,
  colorBy = "GeneExpressionMatrix",
  name = markerGenes,
  embedding = "UMAP_Harmony_Combined"
)

Fplot1 <- lapply(FP[1:2], function(x){
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
FP1 <- do.call(cowplot::plot_grid, c(list(ncol = 2),Fplot1))

Fplot2 <- lapply(FP[3:4], function(x){
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

Fplot3 <- lapply(FP[5:6], function(x){
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

Fplot4 <- lapply(FP[7:8], function(x){
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

Fplot5 <- lapply(FP[9:10], function(x){
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

Fplot6 <- lapply(FP[11:13], function(x){
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

Fplot7 <- lapply(FP[14], function(x){
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
FP7 <- do.call(cowplot::plot_grid, c(list(ncol = 1),Fplot7))

Fplot8 <- lapply(FP[15:17], function(x){
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

Fplot9 <- lapply(FP[18:20], function(x){
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

Fplot10 <- lapply(FP[21:23], function(x){
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

Fplot11 <- lapply(FP[24:25], function(x){
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

Fplot12 <- lapply(FP[26:27], function(x){
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

Fplot13 <- lapply(FP[28:30], function(x){
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
FP13 <- do.call(cowplot::plot_grid, c(list(ncol = 2),Fplot13))

Fplot14 <- lapply(FP[31:35], function(x){
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
FP14 <- do.call(cowplot::plot_grid, c(list(ncol = 2),Fplot14))

Fplot15 <- lapply(FP[36:37], function(x){
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
FP15 <- do.call(cowplot::plot_grid, c(list(ncol = 2),Fplot15))

Fplot16 <- lapply(FP[38:39], function(x){
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
FP16 <- do.call(cowplot::plot_grid, c(list(ncol = 2),Fplot16))

Fplot17 <- lapply(FP[40:41], function(x){
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
FP17 <- do.call(cowplot::plot_grid, c(list(ncol = 2),Fplot17))

plotPDF(FP1, FP2, FP3, FP4, FP5, FP6, FP7, FP8, FP9, FP10, FP11, FP12,
        FP13, FP14, FP15, FP16, FP17,
        name = "FeaturePlots_E775_full_HarmonyComb_UMAP", 
        ArchRProj = proj_Mef2c_v13_E775, 
        addDOC = F
)

#Now label E7.75 - Full Dataset Cluster Labels ------------------------------

#C1 = Visceral Endoderm
#C2 = Blood
#C3 = Trophoblast/Placental
#C4 = Endothelium/Endocardium
#C5 = Extraembryonic Mesoderm
#C6 = Extraembryonic Ectoderm 1
#C7 = Definitive Endoderm
#C8 = Axial Mesoderm
#C9 = CPs/CMs
#C10 = Anterior Mesoderm
#C11 = Posterior Mesoderm
#C12 = Allantois
#C13 = Extraembryonic Ectoderm 2
#C14 = C14  
#C15 = Ectoderm 1
#C16 = Extraembryonic Ectoderm 3
#C17 = Ectoderm 2
#C18 = Ectoderm 3 / PS

oldlabels <- c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10",
               "C11", "C12", "C13", "C14", "C15", "C16", "C17", "C18"
               )
newlabels <- c("VE", "Blood", "T/P", "Endo", "ExM", "ExEct1", "DE",
               "AxM", "CPs/CMs", "AntM", "PostM", "Allantois", "ExEct2",
               "C14", "Ect1", "ExEct3", "Ect2", "Ect3/PS"
               )

proj_Mef2c_v13_E775$Cluster_Labels <- mapLabels(proj_Mef2c_v13_E775$Clus_Harm_Comb_res0.5, newLabels = newlabels, oldLabels = oldlabels)

u21 <- plotEmbedding(
  proj_Mef2c_v13_E775, 
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
        ArchRProj = proj_Mef2c_v13_E775, 
        addDOC = F
        )


###--------------Subset E7.75--------------------------------------------------
#Subset clusters of interest, namely CPs/CMs, Anterior Mesoderm, Extraembryonic Mesoderm

subset_ids <- which(proj_Mef2c_v13_E775$Cluster_Labels == "CPs/CMs" | 
                    proj_Mef2c_v13_E775$Cluster_Labels == "AntM" |
                    proj_Mef2c_v13_E775$Cluster_Labels == "ExM"  
                    )
subset_cells <- proj_Mef2c_v13_E775$cellNames[subset_ids]
proj_Mef2c_v13_E775_subset <- subsetArchRProject(
                            ArchRProj = proj_Mef2c_v13_E775,
                            cells = subset_cells,
                            outputDirectory = "~/Mef2c_ArchR_working/Mef2c_v13_E775_subset",
                            force = TRUE,
                            dropCells = TRUE
                            )

proj_Mef2c_v13_E775_subset
#Subset project contains 2644 cells


###--------------Subset Dim Reduction - LSI ATAC -----------------------------------
proj_Mef2c_v13_E775_subset <- addIterativeLSI(
  ArchRProj = proj_Mef2c_v13_E775_subset, 
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
  force = TRUE
)

###--------------Subset Dim Reduction - LSI RNA -----------------------------------
proj_Mef2c_v13_E775_subset <- addIterativeLSI(
  ArchRProj = proj_Mef2c_v13_E775_subset, 
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
  force = TRUE
)

###------------Subset Combined Dims-------------------------------------------------
proj_Mef2c_v13_E775_subset <- addCombinedDims(
  ArchRProj = proj_Mef2c_v13_E775_subset, 
  reducedDims = c("LSI_ATAC_sub", "LSI_RNA_sub"), 
  name =  "LSI_Combined_sub"
)

###------------Harmony Correction for Combined Dims-------------------------
proj_Mef2c_v13_E775_subset <- addHarmony(
  ArchRProj = proj_Mef2c_v13_E775_subset, 
  reducedDims = "LSI_Combined_sub", 
  name = "Harmony_Combined_sub", 
  groupBy = "Sample", 
  force = T
)

###-------------Subset UMAPs--------------------------------------------------------

#ATAC UMAP
proj_Mef2c_v13_E775_subset <- addUMAP(
  ArchRProj = proj_Mef2c_v13_E775_subset, 
  reducedDims = "LSI_ATAC_sub", 
  name = "UMAP_ATAC_sub", 
  minDist = 0.5, 
  force = TRUE
)

#RNA UMAP
proj_Mef2c_v13_E775_subset <- addUMAP(
  ArchRProj = proj_Mef2c_v13_E775_subset, 
  reducedDims = "LSI_RNA_sub", 
  name = "UMAP_RNA_sub", 
  minDist = 0.5, 
  force = TRUE
)

#Combined UMAP
proj_Mef2c_v13_E775_subset <- addUMAP(
  ArchRProj = proj_Mef2c_v13_E775_subset, 
  reducedDims = "LSI_Combined_sub", 
  name = "UMAP_Combined_sub", 
  minDist = 0.5, 
  force = TRUE
)

#Harmony Combined UMAP
proj_Mef2c_v13_E775_subset <- addUMAP(
  ArchRProj = proj_Mef2c_v13_E775_subset, 
  reducedDims = "Harmony_Combined_sub", 
  name = "UMAP_Harmony_Combined_sub", 
  minDist = 0.5, 
  force = TRUE
)

#Add clusters from Harmony Combined dims 
proj_Mef2c_v13_E775_subset <- addClusters(
  input = proj_Mef2c_v13_E775_subset, 
  reducedDims = "Harmony_Combined_sub", 
  maxClusters = 30,
  name = "Subset_Clus_Harm_Comb_res2", 
  resolution = 2, 
  force = TRUE
)

#Plot embeddings
u31 <- plotEmbedding(
  ArchRProj = proj_Mef2c_v13_E775_subset, 
  name = "Subset_Clus_Harm_Comb_res2", 
  embedding = "UMAP_ATAC_sub", 
  size = 1, 
  labelAsFactors=F, 
  labelMeans=F,
  colorTitle = "",
  baseSize = 16,
) + theme(legend.text = element_text(size=16)) + guides(colour = guide_legend(override.aes = list(size=3)))

u32 <- plotEmbedding(
  ArchRProj = proj_Mef2c_v13_E775_subset, 
  name = "Subset_Clus_Harm_Comb_res2", 
  embedding = "UMAP_RNA_sub", 
  size = 1, 
  labelAsFactors=F, 
  labelMeans=F,
  colorTitle = "",
  baseSize = 16,
) + theme(legend.text = element_text(size=16)) + guides(colour = guide_legend(override.aes = list(size=3)))

u33 <- plotEmbedding(
  ArchRProj = proj_Mef2c_v13_E775_subset, 
  name = "Subset_Clus_Harm_Comb_res2", 
  embedding = "UMAP_Combined_sub", 
  size = 1, 
  labelAsFactors=F, 
  labelMeans=F,
  colorTitle = "",
  baseSize = 16,
) + theme(legend.text = element_text(size=16)) + guides(colour = guide_legend(override.aes = list(size=3)))

u34 <- plotEmbedding(
  ArchRProj = proj_Mef2c_v13_E775_subset, 
  name = "Subset_Clus_Harm_Comb_res2", 
  embedding = "UMAP_Harmony_Combined_sub", 
  size = 1, 
  labelAsFactors=F, 
  labelMeans=F,
  colorTitle = "",
  baseSize = 16,
) + theme(legend.text = element_text(size=16)) + guides(colour = guide_legend(override.aes = list(size=3)))

plotPDF(u31, u32, u33, u34, name = "UMAPs_ATAC_RNA_Combined_HarmonyComb_E775_Subset_res2", width = 10, height = 10, ArchRProj = proj_Mef2c_v13_E775_subset, addDOC = F)

jon_cols = c("E7.75_WT1" = "lightskyblue",
             "E7.75_KO1" = "blue",
             "E7.75_KO2" = "darkblue",
             "E7.75_WT2" = "turquoise1"
             )

u35 <- plotEmbedding(
  proj_Mef2c_v13_E775_subset, 
  name = "SampleNames", 
  embedding = "UMAP_ATAC_sub", 
  size = 1, 
  labelAsFactors=F, 
  labelMeans=F, 
  pal = jon_cols,
  colorTitle = ""
) + theme(legend.text = element_text(size=10)) + guides(colour = guide_legend(override.aes = list(size=3)))

u36 <- plotEmbedding(
  proj_Mef2c_v13_E775_subset, 
  name = "SampleNames", 
  embedding = "UMAP_RNA_sub", 
  size = 1, 
  labelAsFactors=F, 
  labelMeans=F, 
  pal = jon_cols,
  colorTitle = ""
) + theme(legend.text = element_text(size=10)) + guides(colour = guide_legend(override.aes = list(size=3)))

u37 <- plotEmbedding(
  proj_Mef2c_v13_E775_subset, 
  name = "SampleNames", 
  embedding = "UMAP_Combined_sub", 
  size = 1, 
  labelAsFactors=F, 
  labelMeans=F, 
  pal = jon_cols,
  colorTitle = ""
) + theme(legend.text = element_text(size=10)) + guides(colour = guide_legend(override.aes = list(size=3)))

u38 <- plotEmbedding(
  proj_Mef2c_v13_E775_subset, 
  name = "SampleNames", 
  embedding = "UMAP_Harmony_Combined_sub", 
  size = 1, 
  labelAsFactors=F, 
  labelMeans=F, 
  pal = jon_cols,
  colorTitle = ""
) + theme(legend.text = element_text(size=10)) + guides(colour = guide_legend(override.aes = list(size=3)))

plotPDF(u35, u36, u37, u38, name = "UMAPs_ATAC_RNA_Combined_HarmonyComb_E775_Subset_bySample", ArchRProj = proj_Mef2c_v13_E775_subset, addDOC = F)


###---------------Identifying Subset Clusters-------------------------------

#First check Seurat Labels on ArchR UMAPs------------------------------------

library(Seurat)
#Load in Seurat object, in this case "mef2c_v13_E775_subset"

#Get barcodes from ArchR project
archr.barcodes.list <- as.data.frame(proj_Mef2c_v13_E775_subset$cellNames)[,1]
#Get barcodes and cell type labels from Seurat project
seurat.df <- as.data.frame(mef2c_v13_E775_subset$cell_type_subset) 
#Make Seurat barcodes match ArchR formatting
corrected_barcodes <- gsub("_","#",rownames(seurat.df))
rownames(seurat.df) <- corrected_barcodes
seurat.barcodes.list <- rownames(seurat.df)
#Find matching barcodes
matches <- intersect(archr.barcodes.list, seurat.barcodes.list)
#For matched barcodes, add Seurat cluster labels to ArchR metadata
proj_Mef2c_v13_E775_subset$seurat_labels_subset <- NA
proj_Mef2c_v13_E775_subset@cellColData[matches,"seurat_labels_subset"] <- seurat.df[matches,1]

u51 <- plotEmbedding(
  proj_Mef2c_v13_E775_subset, 
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
        name = "UMAP_HarmonyComb_E775_subset_seurat_labels", 
        ArchRProj = proj_Mef2c_v13_E775_subset, 
        addDOC = F
)

#Then look at ArchR feature plots on ArchR UMAPs-----------------------------
markerGenes_sub <- c(
  "Hand1", "Bmp4",                      #ExMeso
  "Sobp", "Otx2", "Cyp26c1",            #Cranial Meso
  "Isl1", "Fgf8", "Fgf10", "Mef2c",     #SHF
  "Gata4", "Gata6",                     #Lateral Mesoderm
  "Fst", "Meox1",                       #Somitic Meso
  "Nkx1-2", "Sox2", "T",                #NMPs
  "Tbx1",                               #Paraxial Meso
  "Tbx5", "Nkx2-5", "Mab21l2",          #FHF/JCF
  "Runx1", "Itga4", "Cd44",             #HSCs
  "Tnnt2", "Ttn",                       #CMs
  "Slc7a8", "Slc2a2", "Cubn", "Amn"     #KPs                             
)

FPs <- plotEmbedding(
  ArchRProj = proj_Mef2c_v13_E775_subset,
  colorBy = "GeneExpressionMatrix",
  name = markerGenes_sub,
  embedding = "UMAP_Harmony_Combined_sub",
  plotAs = "points",
  size = 1
)

Fplot1s <- lapply(FPs[1:2], function(x){
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

Fplot2s <- lapply(FPs[3:5], function(x){
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

Fplot3s <- lapply(FPs[6:9], function(x){
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

Fplot4s <- lapply(FPs[10:11], function(x){
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

Fplot5s <- lapply(FPs[12:13], function(x){
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

Fplot6s <- lapply(FPs[14:16], function(x){
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

Fplot7s <- lapply(FPs[17], function(x){
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
FP7s <- do.call(cowplot::plot_grid, c(list(ncol = 1),Fplot7s))

Fplot8s <- lapply(FPs[18:20], function(x){
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

Fplot9s <- lapply(FPs[21:23], function(x){
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

Fplot10s <- lapply(FPs[24:25], function(x){
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

Fplot11s <- lapply(FPs[26:29], function(x){
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
FP11s <- do.call(cowplot::plot_grid, c(list(ncol = 2),Fplot11s))

plotPDF(FP1s, FP2s, FP3s, FP4s, FP5s, FP6s, FP7s, FP8s, FP9s, FP10s, FP11s,
        name = "FeaturePlots_E775_subset_HarmonyComb_UMAP", 
        ArchRProj = proj_Mef2c_v13_E775_subset, 
        addDOC = F
)

#Now label E7.75 - Subset Cluster Labels ------------------------------

#C1 = SHF1
#C2 = LPM
#C3 = SoM1
#C4 = SoM2
#C5 = PrxM
#C6 = CrM1
#C7 = CrM2
#C8 = HSCs
#C9 = ExM1
#C10 = ExM2
#C11 = EXM3
#C12 = ExM4
#C13 = ExM5
#C14 = ExM6
#C15 = CMs
#C16 = SHF2/JCF
#C17 = FHF

oldlabels2 <- c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10",
               "C11", "C12", "C13", "C14", "C15", "C16", "C17"
               )
newlabels2 <- c("SHF1", "LPM", "SoM1", "SoM2", 
               "PrxM", "CrM1", "CrM2", "HSCs",
               "ExM1", "ExM2", "ExM3", "ExM4",  
               "ExM5", "ExM6", "CMs", "SHF2/JCF", "FHF"
               )

proj_Mef2c_v13_E775_subset$Cluster_Labels_Subset <- mapLabels(proj_Mef2c_v13_E775_subset$Subset_Clus_Harm_Comb_res2, newLabels = newlabels2, oldLabels = oldlabels2)

u40 <- plotEmbedding(
  proj_Mef2c_v13_E775_subset, 
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
        name = "UMAP_HarmonyComb_E775_Subset_res2_labeled", 
        ArchRProj = proj_Mef2c_v13_E775_subset, 
        addDOC = F
)


###----------Prepare to Create Pseudobulk Replicates--------------------------------

#First, create identity labels with CPs grouped together, and all other cell types with multiple clusters
#i.e. combine FHF, SHF1, SHF2/JCF clusters - excluding CMs b/c only contain WT cells
oldlabels3 <- c("SHF1", "LPM", "SoM1", "SoM2", 
                "PrxM", "CrM1", "CrM2", "HSCs",
                "ExM1", "ExM2", "ExM3", "ExM4",  
                "ExM5", "ExM6", "CMs", "SHF2/JCF", "FHF"
)
newlabels3 <- c("CPs", "LPM", "SoM", "SoM", 
                "PrxM", "CrM", "CrM", "HSCs",
                "ExM", "ExM", "ExM", "ExM",  
                "ExM", "ExM", "CMs", "CPs", "CPs"
)

proj_Mef2c_v13_E775_subset$Cluster_Labels_Subset_celltypesgrouped <- mapLabels(proj_Mef2c_v13_E775_subset$Cluster_Labels_Subset, newLabels = newlabels3, oldLabels = oldlabels3)

u41 <- plotEmbedding(
  proj_Mef2c_v13_E775_subset, 
  name = "Cluster_Labels_Subset_celltypesgrouped", 
  embedding = "UMAP_Harmony_Combined_sub", 
  size = 1,
  labelAsFactors=F, 
  labelMeans=T,
  labelSize = 4,
  baseSize = 10,
  colorTitle = ""
) + theme(legend.text = element_text(size=10)) + guides(colour = guide_legend(override.aes = list(size=3)))

plotPDF(u41, 
        name = "UMAP_HarmonyComb_E775_Subset_res2_labeled_celltypesgrouped", 
        ArchRProj = proj_Mef2c_v13_E775_subset, 
        addDOC = F
)

#Next, label cells by genotype - samples 01 and 04 are WT, samples 02 and 03 are KO
geno_WT <- gsub("sample01|sample04", "WT", proj_Mef2c_v13_E775_subset$Sample)
geno <- gsub("sample02|sample03", "KO", geno_WT)
proj_Mef2c_v13_E775_subset$Genotype <- geno

#Next, create CellTypeByGenotype column for grouping
proj_Mef2c_v13_E775_subset <- addCellColData(ArchRProj = proj_Mef2c_v13_E775_subset, 
                                              data = paste0(proj_Mef2c_v13_E775_subset$Cluster_Labels_Subset_celltypesgrouped,"_x_",
                                                            proj_Mef2c_v13_E775_subset$Genotype), 
                                              name = "CellTypeByGenotype", 
                                              cells = getCellNames(proj_Mef2c_v13_E775_subset), 
                                              force = TRUE
                                              )

#Create CellTypeBySample column to allow easy viewing of number of cells per cluster
proj_Mef2c_v13_E775_subset <- addCellColData(ArchRProj = proj_Mef2c_v13_E775_subset, 
                                              data = paste0(proj_Mef2c_v13_E775_subset$Cluster_Labels_Subset_celltypesgrouped,"_x_",
                                                            proj_Mef2c_v13_E775_subset$Sample), 
                                              name = "CellTypeBySample", 
                                              cells = getCellNames(proj_Mef2c_v13_E775_subset), 
                                              force = TRUE
                                              )


###------------Add Group Coverages aka Create Pseudobulk Reps------------------

#Check cell count tables to inform pseudobulking params
table(proj_Mef2c_v13_E775_subset$CellTypeByGenotype)
table(proj_Mef2c_v13_E775_subset$CellTypeBySample)

#Now addGroupCoverages
proj_Mef2c_v13_E775_subset <- addGroupCoverages(
  ArchRProj = proj_Mef2c_v13_E775_subset,
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
proj_Mef2c_v13_E775_subset <- addReproduciblePeakSet(
  ArchRProj = proj_Mef2c_v13_E775_subset,
  groupBy = "CellTypeByGenotype",
  #Note: with reproducibility = 2, both my pseodobulk reps for each group will
  #be required to contain a peak call at the locus - this will yield strict peak calls
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

#Group info used for peak calls - from addGroupCoverages
#              Group nCells nCellsUsed nReplicates nMin nMax maxPeaks
#CMs_x_WT   CMs_x_WT     43         43           2   20   23    21500
#CPs_x_KO   CPs_x_KO    270        270           2   77  193   135000
#CPs_x_WT   CPs_x_WT    121        121           2   26   95    60500
#CrM_x_KO   CrM_x_KO    243        243           2   77  166   121500
#CrM_x_WT   CrM_x_WT    193        193           2   65  128    96500
#ExM_x_KO   ExM_x_KO    565        565           2  181  384   150000
#ExM_x_WT   ExM_x_WT    424        424           2   35  389   150000
#HSCs_x_KO HSCs_x_KO     36         33           2   20   20    16500
#HSCs_x_WT HSCs_x_WT     52         52           2   20   32    26000
#LPM_x_KO   LPM_x_KO     47         47           2   20   27    23500
#LPM_x_WT   LPM_x_WT     69         69           2   30   39    34500
#PrxM_x_KO PrxM_x_KO     88         88           2   24   64    44000
#PrxM_x_WT PrxM_x_WT     66         66           2   29   37    33000
#SoM_x_KO   SoM_x_KO    162        162           2   74   88    81000
#SoM_x_WT   SoM_x_WT    265        265           2  104  161   132500

#addPeakMatrix
proj_Mef2c_v13_E775_subset <- addPeakMatrix(proj_Mef2c_v13_E775_subset)

#Check that PeakMatrix now appears in available matrices
getAvailableMatrices(proj_Mef2c_v13_E775_subset)


###------------Export Group Bigwigs (optional)----------------------------------

getGroupBW(
  ArchRProj = proj_Mef2c_v13_E775_subset,
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
proj_Mef2c_v13_E775_subset <- addMotifAnnotations(ArchRProj = proj_Mef2c_v13_E775_subset, motifSet = "vierstra", collection = "archetype", annoName = "Motif", force = T)


###-------------Add Peak2Gene Linkages----------------------------------------
proj_Mef2c_v13_E775_subset <- addPeak2GeneLinks(
  ArchRProj = proj_Mef2c_v13_E775_subset,
  reducedDims = "Harmony_Combined_sub",
  useMatrix = "GeneExpressionMatrix"
)


###--------------Save and Load ArchR Project------------------------------------
saveArchRProject(ArchRProj = proj_Mef2c_v13_E775, overwrite = TRUE, dropCells = FALSE, load = TRUE)
saveArchRProject(ArchRProj = proj_Mef2c_v13_E775_subset, overwrite = TRUE, dropCells = FALSE, load = TRUE)

#to load
proj_Mef2c_v13_E775 <- loadArchRProject("~/Mef2c_ArchR_working/Mef2c_v13_E775_ATAC_and_GEX")
proj_Mef2c_v13_E775_subset <- loadArchRProject("~/Mef2c_ArchR_working/Mef2c_v13_E775_subset")
setwd("~/Mef2c_v13_ArchR_working")



