###----------------Description-------------------------------------

#This script is used to perform integrated analyses of snRNA-seq and snATAC-seq
#data in ArchR for the four E9 Mef2c Multiome samples 

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

proj_Mef2c_v13_E9_ATAC_only <- ArchRProject(
  ArrowFiles = c("Original_Arrows/sample05.arrow", 
                 "Original_Arrows/sample06.arrow", 
                 "Original_Arrows/sample07.arrow", 
                 "Original_Arrows/sample08.arrow"
                 ),
  outputDirectory = "Mef2c_E9_ATAC_only",
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

seRNA_Mef2c_v13_E9 <- import10xFeatureMatrix(
  input = c("Mef2c_inputs/sample05_filtered_feature_bc_matrix.h5", 
            "Mef2c_inputs/sample06_filtered_feature_bc_matrix.h5", 
            "Mef2c_inputs/sample07_filtered_feature_bc_matrix.h5", 
            "Mef2c_inputs/sample08_filtered_feature_bc_matrix.h5"
            ),
  names = c("sample05", "sample06", "sample07", "sample08"),
  strictMatch = FALSE,
  verbose = TRUE,
  featureType = "Gene Expression"
)

#Note: at this point there are 50505 cells in the ArchR project that passed
#QC based on the ATAC data and  52768 cells in the seRNA ojbect with gene 
#expression data. We will need to subset ArchR project to contain only cells 
#that have info in both.
overlap <- intersect(proj_Mef2c_v13_E9_ATAC_only$cellNames, seRNA_Mef2c_v13_E9@colData@rownames)
#There are 47,845 cells that exist in both the ArchR project and the seRNA object

#Add the seRNA data to ArchR project
#Note: will get warning that not all cells in input exist in seRNA - that's 
#okay, we will subset the overlapping cells next
proj_Mef2c_v13_E9_ATAC_only <- addGeneExpressionMatrix(
  input = proj_Mef2c_v13_E9_ATAC_only, 
  seRNA = seRNA_Mef2c_v13_E9, 
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

proj_Mef2c_v13_E9 <- subsetArchRProject(
  ArchRProj = proj_Mef2c_v13_E9_ATAC_only,
  cells = overlap,
  outputDirectory = "Mef2c_v13_E9_ATAC_and_GEX",
  dropCells = TRUE,
  force = TRUE
)

#Inspect project - should be 47,845 cells
proj_Mef2c_v13_E9

#IMPORTANT: Re-install "release_1.0.2" branch before proceeding
detach("package:ArchR", unload = TRUE)
devtools::install_github("GreenleafLab/ArchR", ref="release_1.0.2", repos = BiocManager::repositories())
library(ArchR)


###---------QC and additional TSSEnrichment and nFrag Filtering----------------

#Based on QC metrics, I think I want to further refine cutoff for TSSEnrichment
#Let's build the QC plot TSSEnrichment vs. log10 Unique Frags now that cells have
#been filtered

#Get nFrags and TSSEnrichment scores
df <- getCellColData(proj_Mef2c_v13_E9, select = c("log10(nFrags)", "TSSEnrichment"))
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
ggsave("Mef2c_v13_E9_TSSbyLog10nFrags_initial.pdf", plot = p, device = "pdf", path = "~/Mef2c_ArchR_working/QualityControl", width = 5, height = 5, dpi = 300)

#Also want to filter out cells with high/low GEX UMI/Gene counts - let's look at GEX data
df2 <- getCellColData(proj_Mef2c_v13_E9, select = c("Gex_nUMI", "Gex_nGenes"))
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
ggsave("Mef2c_v13_E9_nUMIbynGenes_initial.pdf", plot = p2, device = "pdf", path = "~/Mef2c_ArchR_working/QualityControl", width = 5, height = 5, dpi = 300)

#Before filtering, visualize subset of cells with new cutoffs
#TSS Enrichment >= 13, TSS Enrichment <= 30, nFrags >= 4000 (log10 3.6), Gex_nGenes >= 500, Gex_nGenes <= 9000
idxPass <- which(proj_Mef2c_v13_E9$TSSEnrichment >= 13 & 
                   proj_Mef2c_v13_E9$TSSEnrichment <= 30 & 
                   proj_Mef2c_v13_E9$nFrags >= 4000 & 
                   proj_Mef2c_v13_E9$Gex_nGenes >= 500 &   #Same cutoffs as Seurat analysis
                   proj_Mef2c_v13_E9$Gex_nGenes <= 9000    #Same cutoffs as Seurat analysis
                 )     
cellsPass <- proj_Mef2c_v13_E9$cellNames[idxPass]
df3 <- getCellColData(proj_Mef2c_v13_E9[cellsPass, ], select = c("log10(nFrags)", "TSSEnrichment"))
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
ggsave("Mef2c_v13_E9_TSSbyLog10nFrags_filtered.pdf", plot = p3, device = "pdf", path = "~/Mef2c_ArchR_working/QualityControl", width = 5, height = 5, dpi = 300)

#And the GEX plot
df4 <- getCellColData(proj_Mef2c_v13_E9[cellsPass, ], select = c("Gex_nUMI", "Gex_nGenes"))
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
ggsave("Mef2c_v13_E9_nUMIbynGenes_filtered.pdf", plot = p4, device = "pdf", path = "~/Mef2c_ArchR_working/QualityControl", width = 5, height = 5, dpi = 300)

#Looks good, now filter these cells from project 
proj_Mef2c_v13_E9 <- proj_Mef2c_v13_E9[cellsPass, ]
proj_Mef2c_v13_E9
#Now at 38985 cells


###----------------Doublet Filtering-----------------------------------------

#First add doublet scores
dir.create("~/Mef2c_ArchR_working/Mef2c_v13_E9_ATAC_and_GEX/DoubletSummaries")
proj_Mef2c_v13_E9 <- addDoubletScores(proj_Mef2c_v13_E9, outDir = paste(getOutputDirectory(proj_Mef2c_v13_E9), "DoubletSummaries", sep = "/"))

#Filter doublets
proj_Mef2c_v13_E9 <- filterDoublets(proj_Mef2c_v13_E9, filterRatio = 1)

#Filtering 3961 cells from ArchRProject!
  #sample05 : 991 of 9955 (10%)
  #sample06 : 408 of 6393 (6.4%)
  #sample07 : 1264 of 11243 (11.2%)
  #sample08 : 1298 of 11394 (11.4%)

#Now at 35,024 cells 

saveArchRProject(ArchRProj = proj_Mef2c_v13_E9, overwrite = TRUE, dropCells = TRUE, load = TRUE)


###--------------Dim Reduction - LSI ATAC -----------------------------------
proj_Mef2c_v13_E9 <- addIterativeLSI(
  ArchRProj = proj_Mef2c_v13_E9, 
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
proj_Mef2c_v13_E9 <- addIterativeLSI(
  ArchRProj = proj_Mef2c_v13_E9, 
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
proj_Mef2c_v13_E9 <- addCombinedDims(
  ArchRProj = proj_Mef2c_v13_E9, 
  reducedDims = c("LSI_ATAC", "LSI_RNA"), 
  name =  "LSI_Combined"
)

###------------Harmony Correction for Combined Dims-------------------------
proj_Mef2c_v13_E9 <- addHarmony(
  ArchRProj = proj_Mef2c_v13_E9, 
  reducedDims = "LSI_Combined", 
  name = "Harmony_Combined", 
  groupBy = "Sample", 
  force = T
)

###-------------UMAPs--------------------------------------------------------
#ATAC UMAP
proj_Mef2c_v13_E9 <- addUMAP(
  ArchRProj = proj_Mef2c_v13_E9, 
  reducedDims = "LSI_ATAC", 
  name = "UMAP_ATAC", 
  minDist = 0.5, 
  force = TRUE,
  seed = 1
)

#RNA UMAP
proj_Mef2c_v13_E9 <- addUMAP(
  ArchRProj = proj_Mef2c_v13_E9, 
  reducedDims = "LSI_RNA", 
  name = "UMAP_RNA", 
  minDist = 0.5, 
  force = TRUE,
  seed = 1
)

#Combined UMAP
proj_Mef2c_v13_E9 <- addUMAP(
  ArchRProj = proj_Mef2c_v13_E9, 
  reducedDims = "LSI_Combined", 
  name = "UMAP_Combined", 
  minDist = 0.5, 
  force = TRUE,
  seed = 1
)

#Harmony Combined UMAP
proj_Mef2c_v13_E9 <- addUMAP(
  ArchRProj = proj_Mef2c_v13_E9, 
  reducedDims = "Harmony_Combined", 
  name = "UMAP_Harmony_Combined", 
  minDist = 0.5, 
  force = TRUE,
  seed = 1
)

#Add clusters from Harmony Combined dims
proj_Mef2c_v13_E9 <- addClusters(
  input = proj_Mef2c_v13_E9, 
  reducedDims = "Harmony_Combined", 
  maxClusters = 30,
  name = "Clus_Harm_Comb_res0.5", 
  resolution = 0.5, 
  force = TRUE,
  seed = 1
)

#Plot embeddings
u1 <- plotEmbedding(
  ArchRProj = proj_Mef2c_v13_E9, 
  name = "Clus_Harm_Comb_res0.5", 
  embedding = "UMAP_ATAC", 
  size = 0.5, 
  labelAsFactors=F, 
  labelMeans=F,
  colorTitle = "",
  baseSize = 16,
) + theme(legend.text = element_text(size=16)) + guides(colour = guide_legend(override.aes = list(size=3)))

u2 <- plotEmbedding(
  ArchRProj = proj_Mef2c_v13_E9, 
  name = "Clus_Harm_Comb_res0.5", 
  embedding = "UMAP_RNA", 
  size = 0.5, 
  labelAsFactors=F, 
  labelMeans=F,
  colorTitle = "",
  baseSize = 16,
) + theme(legend.text = element_text(size=16)) + guides(colour = guide_legend(override.aes = list(size=3)))

u3 <- plotEmbedding(
  ArchRProj = proj_Mef2c_v13_E9, 
  name = "Clus_Harm_Comb_res0.5", 
  embedding = "UMAP_Combined", 
  size = 0.5, 
  labelAsFactors=F, 
  labelMeans=F,
  colorTitle = "",
  baseSize = 16,
) + theme(legend.text = element_text(size=16)) + guides(colour = guide_legend(override.aes = list(size=3)))

u4 <- plotEmbedding(
  ArchRProj = proj_Mef2c_v13_E9, 
  name = "Clus_Harm_Comb_res0.5", 
  embedding = "UMAP_Harmony_Combined", 
  size = 0.5, 
  labelAsFactors=F, 
  labelMeans=F,
  colorTitle = "",
  baseSize = 16,
) + theme(legend.text = element_text(size=16)) + guides(colour = guide_legend(override.aes = list(size=3)))

plotPDF(u1, u2, u3, u4, name = "UMAPs_ATAC_RNA_Combined_HarmonyComb_res0.5", width = 10, height = 10, ArchRProj = proj_Mef2c_v13_E9, addDOC = F)

#UMAPs colored by sample
#First, label cells by sample name - samples 05 and 08 are WT, samples 06 and 07 are KO
sn1 <- gsub("sample05", "E9.0_WT1", proj_Mef2c_v13_E9$Sample)
sn2 <- gsub("sample06", "E9.0_KO1", sn1)
sn3 <- gsub("sample07", "E9.0_KO2", sn2)
sn4 <- gsub("sample08", "E9.0_WT2", sn3)
proj_Mef2c_v13_E9$SampleNames <- sn4

jon_cols = c("E9.0_WT1" = "pink","E9.0_KO1" = "red", "E9.0_KO2" = "red4", "E9.0_WT2" = "plum1")

u5 <- plotEmbedding(
  proj_Mef2c_v13_E9, 
  name = "SampleNames", 
  embedding = "UMAP_ATAC", 
  size = 0.5, 
  labelAsFactors=F, 
  labelMeans=F, 
  pal = jon_cols,
  colorTitle = ""
) + theme(legend.text = element_text(size=10)) + guides(colour = guide_legend(override.aes = list(size=3)))

u6 <- plotEmbedding(
  proj_Mef2c_v13_E9, 
  name = "SampleNames", 
  embedding = "UMAP_RNA", 
  size = 0.5, 
  labelAsFactors=F, 
  labelMeans=F, 
  pal = jon_cols,
  colorTitle = ""
) + theme(legend.text = element_text(size=10)) + guides(colour = guide_legend(override.aes = list(size=3)))

u7 <- plotEmbedding(
  proj_Mef2c_v13_E9, 
  name = "SampleNames", 
  embedding = "UMAP_Combined", 
  size = 0.5, 
  labelAsFactors=F, 
  labelMeans=F, 
  pal = jon_cols,
  colorTitle = ""
) + theme(legend.text = element_text(size=10)) + guides(colour = guide_legend(override.aes = list(size=3)))

u8 <- plotEmbedding(
  proj_Mef2c_v13_E9, 
  name = "SampleNames", 
  embedding = "UMAP_Harmony_Combined", 
  size = 0.5, 
  labelAsFactors=F, 
  labelMeans=F, 
  pal = jon_cols,
  colorTitle = ""
) + theme(legend.text = element_text(size=10)) + guides(colour = guide_legend(override.aes = list(size=3)))

plotPDF(u5, u6, u7, u8, name = "UMAPs_ATAC_RNA_Combined_HarmonyComb_bySample", ArchRProj = proj_Mef2c_v13_E9, addDOC = F)

#Check doublet scores and other QC metrics

#QC Plots
u9 <- plotEmbedding(proj_Mef2c_v13_E9, name = "DoubletScore", embedding = "UMAP_Harmony_Combined", size = 0.5, labelAsFactors=F, labelMeans=F)
u10 <- plotEmbedding(proj_Mef2c_v13_E9, name = "DoubletEnrichment", embedding = "UMAP_Harmony_Combined", size = 0.5, labelAsFactors=F, labelMeans=F)
u11 <- plotEmbedding(proj_Mef2c_v13_E9, name = "nFrags", embedding = "UMAP_Harmony_Combined", size = 0.5, labelAsFactors=F, labelMeans=F)
u12 <- plotEmbedding(proj_Mef2c_v13_E9, name = "TSSEnrichment", embedding = "UMAP_Harmony_Combined", size = 0.5, labelAsFactors=F, labelMeans=F)
u13 <- plotEmbedding(proj_Mef2c_v13_E9, name = "NucleosomeRatio", embedding = "UMAP_Harmony_Combined", size = 0.5, labelAsFactors=F, labelMeans=F)
u14 <- plotEmbedding(proj_Mef2c_v13_E9, name = "Gex_nUMI", embedding = "UMAP_Harmony_Combined", size = 0.5, labelAsFactors=F, labelMeans=F)
u15 <- plotEmbedding(proj_Mef2c_v13_E9, name = "Gex_nGenes", embedding = "UMAP_Harmony_Combined", size = 0.5, labelAsFactors=F, labelMeans=F)
plotPDF(u9, u10, u11, u12, u13, u14, u15, name = "additional_QC_Combined", ArchRProj = proj_Mef2c_v13_E9, addDOC = F)


###--------------Identifying Clusters----------------------------------------

#First check Seurat Labels on ArchR UMAPs------------------------------------

library(Seurat)
#Load in Seurat object, in this case "mef2c_v13_E9"  

#Get barcodes from ArchR project
archr.barcodes.list <- as.data.frame(proj_Mef2c_v13_E9$cellNames)[,1]
#Get barcodes and cell type labels from Seurat project
seurat.df <- as.data.frame(mef2c_v13_E9$cell_type) 
#Make Seurat barcodes match ArchR formatting
corrected_barcodes <- gsub("_","#",rownames(seurat.df))
rownames(seurat.df) <- corrected_barcodes
seurat.barcodes.list <- rownames(seurat.df)
#Find matching barcodes
matches <- intersect(archr.barcodes.list, seurat.barcodes.list)
#For matched barcodes, add Seurat cluster labels to ArchR metadata
proj_Mef2c_v13_E9$seurat_labels <- NA
proj_Mef2c_v13_E9@cellColData[matches,"seurat_labels"] <- seurat.df[matches,1]

u50 <- plotEmbedding(
  proj_Mef2c_v13_E9, 
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
        name = "UMAP_Comb_E9_seurat_labels", 
        ArchRProj = proj_Mef2c_v13_E9, 
        addDOC = F
)


#Then look at ArchR feature plots on ArchR UMAPs-----------------------------

markerGenes <- c(
  "Pbx1", "Foxp2", "Dach1", "Meis1",    #Mixed Mesoderm
  "Sox2", "Fgf14", "Nrcam",             #Neural Tube
  "Meox1", "Fst",                       #Somitic Mesoderm
  "Sox17", "Epcam", "Cdh1", "Cldn6",    #Definitive Endoderm
  "Isl1", "Tbx1", "Fgf10",              #SHF
  "Tnnt2", "Ttn", "Myl7", "Actc1",      #CMs
  "Wt1", "Tbx18",                       #Proepicardium
  "Tfap2b", "Sox10", "Ets1",            #Neural Crest
  "Cdh5", "Kdr", "Flt1",                #Endothelium
  "Hba-a1", "Hba-a2", "Hbb-y",          #Blood
  "Ttr", "Afp", "Apoa1",                #Yolk Sac
  "Wnt6", "Flt1",                       #Extraembryonic Endoderm 
  "Nog", "Shh"                          #Notochord
)

FP <- plotEmbedding(
  ArchRProj = proj_Mef2c_v13_E9,
  colorBy = "GeneExpressionMatrix",
  name = markerGenes,
  embedding = "UMAP_Harmony_Combined"
)

Fplot1 <- lapply(FP[1:4], function(x){
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

Fplot2 <- lapply(FP[5:7], function(x){
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

Fplot3 <- lapply(FP[8:9], function(x){
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

Fplot4 <- lapply(FP[10:13], function(x){
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
FP4 <- do.call(cowplot::plot_grid, c(list(ncol = 3),Fplot4))

Fplot5 <- lapply(FP[14:16], function(x){
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
FP5 <- do.call(cowplot::plot_grid, c(list(ncol = 3),Fplot5))

Fplot6 <- lapply(FP[17:20], function(x){
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

Fplot7 <- lapply(FP[21:22], function(x){
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

Fplot8 <- lapply(FP[23:25], function(x){
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
FP8 <- do.call(cowplot::plot_grid, c(list(ncol = 1),Fplot8))

Fplot9 <- lapply(FP[26:28], function(x){
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
FP9 <- do.call(cowplot::plot_grid, c(list(ncol = 3),Fplot9))

Fplot10 <- lapply(FP[29:31], function(x){
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
FP10 <- do.call(cowplot::plot_grid, c(list(ncol = 3),Fplot10))

Fplot11 <- lapply(FP[32:34], function(x){
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
FP11 <- do.call(cowplot::plot_grid, c(list(ncol = 3),Fplot11))

Fplot12 <- lapply(FP[35:36], function(x){
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

Fplot13 <- lapply(FP[37:38], function(x){
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

plotPDF(FP1, FP2, FP3, FP4, FP5, FP6, FP7, FP8, FP9, FP10, FP11, FP12, FP13, 
        name = "FeaturePlots_E9_full_Comb_UMAP", 
        ArchRProj = proj_Mef2c_v13_E9, 
        addDOC = F
)

#Now label E9 - Full Dataset Cluster Labels ------------------------------

#C1 = Blood
#C2 = Yolk Sac
#C3 = Pharyngeal Arch Cells (Mixed Germ Layer contributions)
#C4 = Neural Tube 1
#C5 = Neural Tube 2
#C6 = Neural Tube 3
#C7 = Neural Tube 4
#C8 = Neural Tube 5
#C9 = Neural Crest
#C10 = Endoderm 1 
#C11 = Endothelium 
#C12 = CMs 
#C13 = Somitic Mesoderm
#C14 = Pharyngeal Mesoderm 
#C15 = Pericardium/Mesenchyme 
#C16 = SHF
#C17 = Mixed Mesoderm
#C18 = Splanchnic Mesoderm 
#C19 = Endoderm 2 
#C20 = Endoderm 3
#C21 = Endoderm 4
#C22 = Endoderm 5
#C23 = Notochord

oldlabels <- c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10",
               "C11", "C12", "C13", "C14", "C15", "C16", "C17", "C18", "C19", "C20", "C21", "C22", "C23"
               )

newlabels <- c("Blood", "YS", "PhA", "NT1", "NT2", "NT3", "NT4", "NT5",
               "NC", "En1", "Endo", "CMs", "SoM", "PhM",
               "Pe/Me", "SHF", "MixM", "SpM", "En2","En3","En4","En5","Noto"
               )

proj_Mef2c_v13_E9$Cluster_Labels <- mapLabels(proj_Mef2c_v13_E9$Clus_Harm_Comb_res0.5, newLabels = newlabels, oldLabels = oldlabels)

u21 <- plotEmbedding(
  proj_Mef2c_v13_E9, 
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
        ArchRProj = proj_Mef2c_v13_E9, 
        addDOC = F
        )


###--------------Subset E9--------------------------------------------------
#Subset clusters of interest, namely CMs, SHF, and Pe/Me

subset_ids <- which(proj_Mef2c_v13_E9$Cluster_Labels == "CMs" | 
                      proj_Mef2c_v13_E9$Cluster_Labels == "SHF" |
                      proj_Mef2c_v13_E9$Cluster_Labels == "Pe/Me" 
                    )
subset_cells <- proj_Mef2c_v13_E9$cellNames[subset_ids]
proj_Mef2c_v13_E9_subset <- subsetArchRProject(
                            ArchRProj = proj_Mef2c_v13_E9,
                            cells = subset_cells,
                            outputDirectory = "~/Mef2c_ArchR_working/Mef2c_v13_E9_subset",
                            force = TRUE,
                            dropCells = TRUE
                            )

proj_Mef2c_v13_E9_subset
#Subset project contains 7299 cells


###--------------Subset Dim Reduction - LSI ATAC -----------------------------------
proj_Mef2c_v13_E9_subset <- addIterativeLSI(
  ArchRProj = proj_Mef2c_v13_E9_subset, 
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
proj_Mef2c_v13_E9_subset <- addIterativeLSI(
  ArchRProj = proj_Mef2c_v13_E9_subset, 
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
proj_Mef2c_v13_E9_subset <- addCombinedDims(
  ArchRProj = proj_Mef2c_v13_E9_subset, 
  reducedDims = c("LSI_ATAC_sub", "LSI_RNA_sub"), 
  name =  "LSI_Combined_sub"
)

###------------Harmony Correction for Combined Dims-------------------------
proj_Mef2c_v13_E9_subset <- addHarmony(
  ArchRProj = proj_Mef2c_v13_E9_subset, 
  reducedDims = "LSI_Combined_sub", 
  name = "Harmony_Combined_sub", 
  groupBy = "Sample", 
  force = T
)

###-------------Subset UMAPs--------------------------------------------------------

#ATAC UMAP
proj_Mef2c_v13_E9_subset <- addUMAP(
  ArchRProj = proj_Mef2c_v13_E9_subset, 
  reducedDims = "LSI_ATAC_sub", 
  name = "UMAP_ATAC_sub", 
  minDist = 0.5,
  seed = 1,
  force = TRUE
)

#RNA UMAP
proj_Mef2c_v13_E9_subset <- addUMAP(
  ArchRProj = proj_Mef2c_v13_E9_subset, 
  reducedDims = "LSI_RNA_sub", 
  name = "UMAP_RNA_sub", 
  minDist = 0.5, 
  seed = 1,
  force = TRUE
)

#Combined UMAP
proj_Mef2c_v13_E9_subset <- addUMAP(
  ArchRProj = proj_Mef2c_v13_E9_subset, 
  reducedDims = "LSI_Combined_sub", 
  name = "UMAP_Combined_sub", 
  minDist = 0.5,
  seed = 1,
  force = TRUE
)

#Harmony Combined UMAP
proj_Mef2c_v13_E9_subset <- addUMAP(
  ArchRProj = proj_Mef2c_v13_E9_subset, 
  reducedDims = "Harmony_Combined_sub", 
  name = "UMAP_Harmony_Combined_sub", 
  minDist = 0.5, 
  seed = 1,
  force = TRUE
)

#Add clusters from Harmony Combined dims 
proj_Mef2c_v13_E9_subset <- addClusters(
  input = proj_Mef2c_v13_E9_subset, 
  reducedDims = "Harmony_Combined_sub", 
  maxClusters = 30,
  name = "Subset_Clus_Harm_Comb_res1", 
  resolution = 1, 
  seed = 1,
  force = TRUE
)

#Plot embeddings
u31 <- plotEmbedding(
  ArchRProj = proj_Mef2c_v13_E9_subset, 
  name = "Subset_Clus_Harm_Comb_res1", 
  embedding = "UMAP_ATAC_sub", 
  size = 1, 
  labelAsFactors=F, 
  labelMeans=F,
  colorTitle = "",
  baseSize = 16,
) + theme(legend.text = element_text(size=16)) + guides(colour = guide_legend(override.aes = list(size=3)))

u32 <- plotEmbedding(
  ArchRProj = proj_Mef2c_v13_E9_subset, 
  name = "Subset_Clus_Harm_Comb_res1", 
  embedding = "UMAP_RNA_sub", 
  size = 1, 
  labelAsFactors=F, 
  labelMeans=F,
  colorTitle = "",
  baseSize = 16,
) + theme(legend.text = element_text(size=16)) + guides(colour = guide_legend(override.aes = list(size=3)))

u33 <- plotEmbedding(
  ArchRProj = proj_Mef2c_v13_E9_subset, 
  name = "Subset_Clus_Harm_Comb_res1", 
  embedding = "UMAP_Combined_sub", 
  size = 1, 
  labelAsFactors=F, 
  labelMeans=F,
  colorTitle = "",
  baseSize = 16,
) + theme(legend.text = element_text(size=16)) + guides(colour = guide_legend(override.aes = list(size=3)))

u34 <- plotEmbedding(
  ArchRProj = proj_Mef2c_v13_E9_subset, 
  name = "Subset_Clus_Harm_Comb_res1", 
  embedding = "UMAP_Harmony_Combined_sub", 
  size = 1, 
  labelAsFactors=F, 
  labelMeans=F,
  colorTitle = "",
  baseSize = 16,
) + theme(legend.text = element_text(size=16)) + guides(colour = guide_legend(override.aes = list(size=3)))

plotPDF(u31, u32, u33, u34, name = "UMAPs_ATAC_RNA_Combined_HarmonyComb_E9_Subset_res1", width = 10, height = 10, ArchRProj = proj_Mef2c_v13_E9_subset, addDOC = F)

jon_cols = c("E9.0_WT1" = "pink","E9.0_KO1" = "red", "E9.0_KO2" = "red4", "E9.0_WT2" = "plum1")

u35 <- plotEmbedding(
  proj_Mef2c_v13_E9_subset, 
  name = "SampleNames", 
  embedding = "UMAP_ATAC_sub", 
  size = 1, 
  labelAsFactors=F, 
  labelMeans=F, 
  pal = jon_cols,
  colorTitle = ""
) + theme(legend.text = element_text(size=10)) + guides(colour = guide_legend(override.aes = list(size=3)))

u36 <- plotEmbedding(
  proj_Mef2c_v13_E9_subset, 
  name = "SampleNames", 
  embedding = "UMAP_RNA_sub", 
  size = 1, 
  labelAsFactors=F, 
  labelMeans=F, 
  pal = jon_cols,
  colorTitle = ""
) + theme(legend.text = element_text(size=10)) + guides(colour = guide_legend(override.aes = list(size=3)))

u37 <- plotEmbedding(
  proj_Mef2c_v13_E9_subset, 
  name = "SampleNames", 
  embedding = "UMAP_Combined_sub", 
  size = 1, 
  labelAsFactors=F, 
  labelMeans=F, 
  pal = jon_cols,
  colorTitle = ""
) + theme(legend.text = element_text(size=10)) + guides(colour = guide_legend(override.aes = list(size=3)))

u38 <- plotEmbedding(
  proj_Mef2c_v13_E9_subset, 
  name = "SampleNames", 
  embedding = "UMAP_Harmony_Combined_sub", 
  size = 1, 
  labelAsFactors=F, 
  labelMeans=F, 
  pal = jon_cols,
  colorTitle = ""
) + theme(legend.text = element_text(size=10)) + guides(colour = guide_legend(override.aes = list(size=3)))

plotPDF(u35, u36, u37, u38, name = "UMAPs_ATAC_RNA_Combined_HarmonyComb_E9_Subset_bySample", ArchRProj = proj_Mef2c_v13_E9_subset, addDOC = F)


###---------------Identifying Subset Clusters-------------------------------

#First check Seurat Labels on ArchR UMAPs------------------------------------

library(Seurat)
#Load in Seurat object, in this case "mef2c_v13_E9_subset"

#Get barcodes from ArchR project
archr.barcodes.list <- as.data.frame(proj_Mef2c_v13_E9_subset$cellNames)[,1]
#Get barcodes and cell type labels from Seurat project
seurat.df <- as.data.frame(mef2c_v13_E9_subset$cell_type_subset) 
#Make Seurat barcodes match ArchR formatting
corrected_barcodes <- gsub("_","#",rownames(seurat.df))
rownames(seurat.df) <- corrected_barcodes
seurat.barcodes.list <- rownames(seurat.df)
#Find matching barcodes
matches <- intersect(archr.barcodes.list, seurat.barcodes.list)
#For matched barcodes, add Seurat cluster labels to ArchR metadata
proj_Mef2c_v13_E9_subset$seurat_labels_subset <- NA
proj_Mef2c_v13_E9_subset@cellColData[matches,"seurat_labels_subset"] <- seurat.df[matches,1]

u51 <- plotEmbedding(
  proj_Mef2c_v13_E9_subset, 
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
        name = "UMAP_HarmonyComb_E9_subset_seurat_labels", 
        ArchRProj = proj_Mef2c_v13_E9_subset, 
        addDOC = F
)

#Then look at ArchR feature plots on ArchR UMAPs-----------------------------
markerGenes_sub <- c(
  "Tnnt2", "Ttn", "Myl7",                       #CMs
  "Nppa", "Myl2", "Irx4",                       #Ventricles
  "Hand1", "Cited1",                            #LV
  "Cck", "Irx2", "Pln",                         #RV
  "Tbx5", "Angpt1", "Vsnl1", "Nr2f2",           #Atria
  "Rspo3", "Bmp2", "Tbx2",                      #AVC
  "Wnt11", "Rgs5", "Sema3c", "Nrp2", "Tdgf1",   #OFT
  "Wt1", "Tbx18",                               #Proepicardium
  "Isl1", "Tbx1", "Fgf8", "Fgf10",              #aSHF
  "Hoxa1", "Hoxb1", "Foxf1", "Aldh1a2"          #pSHF
)

FPs <- plotEmbedding(
  ArchRProj = proj_Mef2c_v13_E9_subset,
  colorBy = "GeneExpressionMatrix",
  name = markerGenes_sub,
  embedding = "UMAP_Harmony_Combined_sub",
  plotAs = "points",
  size = 1
)

Fplot1s <- lapply(FPs[1:3], function(x){
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

Fplot2s <- lapply(FPs[4:6], function(x){
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
FP2s <- do.call(cowplot::plot_grid, c(list(ncol = 3),Fplot2s))

Fplot3s <- lapply(FPs[7:8], function(x){
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
FP3s <- do.call(cowplot::plot_grid, c(list(ncol = 3),Fplot3s))

Fplot4s <- lapply(FPs[9:11], function(x){
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

Fplot5s <- lapply(FPs[12:15], function(x){
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
FP5s <- do.call(cowplot::plot_grid, c(list(ncol = 3),Fplot5s))

Fplot6s <- lapply(FPs[16:18], function(x){
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

Fplot7s <- lapply(FPs[19:23], function(x){
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

Fplot8s <- lapply(FPs[24:25], function(x){
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
FP8s <- do.call(cowplot::plot_grid, c(list(ncol = 3),Fplot8s))

Fplot9s <- lapply(FPs[26:29], function(x){
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
FP9s <- do.call(cowplot::plot_grid, c(list(ncol = 3),Fplot9s))

Fplot10s <- lapply(FPs[30:33], function(x){
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
FP10s <- do.call(cowplot::plot_grid, c(list(ncol = 3),Fplot10s))

plotPDF(FP1s, FP2s, FP3s, FP4s, FP5s, FP6s, FP7s, FP8s, FP9s, FP10s,
        name = "FeaturePlots_E9_subset_HarmonyComb_UMAP", 
        ArchRProj = proj_Mef2c_v13_E9_subset, 
        addDOC = F
)

#Now label E9 - Subset Cluster Labels ------------------------------

#C1 = CMs_A1
#C2 = C2
#C3 = CMs_AVC
#C4 = CMs_OFT
#C5 = CMs_A2
#C6 = CMs_V1
#C7 = CMs_V2
#C8 = Me1
#C9 = Me2
#C10 = VP1
#C11 = Pe
#C12 = VP2
#C13 = SHF1
#C14 = SHF2
#C15 =  SHF3
#C16 = SHF4

oldlabels2 <- c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10",
               "C11", "C12", "C13", "C14", "C15", "C16"
               )

newlabels2 <- c("CMs_A1", "C2", "CMs_AVC", "CMs_OFT", 
               "CMs_A2", "CMs_V1", "CMs_V2", "Me1",
               "Me2", "VP1", "Pe", "VP2", "SHF1", "SHF2", "SHF3", "SHF4"  
               )

proj_Mef2c_v13_E9_subset$Cluster_Labels_Subset <- mapLabels(proj_Mef2c_v13_E9_subset$Subset_Clus_Harm_Comb_res1, newLabels = newlabels2, oldLabels = oldlabels2)

u40 <- plotEmbedding(
  proj_Mef2c_v13_E9_subset, 
  name = "Cluster_Labels_Subset", 
  embedding = "UMAP_Harmony_Combined_sub", 
  size = 0.5,
  labelAsFactors=F, 
  labelMeans=T,
  labelSize = 3,
  baseSize = 10,
  colorTitle = ""
) + theme(legend.text = element_text(size=10)) + guides(colour = guide_legend(override.aes = list(size=3)))

plotPDF(u40, 
        name = "UMAP_HarmonyComb_E9_Subset_res1_labeled", 
        ArchRProj = proj_Mef2c_v13_E9_subset, 
        addDOC = F
)


###----------Prepare to Create Pseudobulk Replicates--------------------------------

#First, create identity labels with all cell type clusters grouped together
#For these analyses, we will group CMs_A and CMs_AVC into a single CMs_IFT Group 

oldlabels3 <- c("CMs_A1", "C2", "CMs_AVC", "CMs_OFT", 
"CMs_A2", "CMs_V1", "CMs_V2", "Me1",
"Me2", "VP1", "Pe", "VP2", "SHF1", "SHF2", "SHF3", "SHF4"  )

newlabels3 <- c("CMs_IFT", "C2", "CMs_IFT", "CMs_OFT", 
                "CMs_IFT", "CMs_V", "CMs_V", "Me",
                "Me", "VP", "Pe", "VP", "SHF", "SHF", "SHF", "SHF"
)

proj_Mef2c_v13_E9_subset$Cluster_Labels_Subset_celltypegroup <- mapLabels(proj_Mef2c_v13_E9_subset$Cluster_Labels_Subset, newLabels = newlabels3, oldLabels = oldlabels3)

u41 <- plotEmbedding(
  proj_Mef2c_v13_E9_subset, 
  name = "Cluster_Labels_Subset_celltypegroup", 
  embedding = "UMAP_Harmony_Combined_sub", 
  size = 0.5,
  labelAsFactors=F, 
  labelMeans=T,
  labelSize = 3,
  baseSize = 10,
  colorTitle = ""
) + theme(legend.text = element_text(size=10)) + guides(colour = guide_legend(override.aes = list(size=3)))

plotPDF(u41, 
        name = "UMAP_HarmonyComb_E9_Subset_res1_labeled_celltypegroup", 
        ArchRProj = proj_Mef2c_v13_E9_subset, 
        addDOC = F
)

#Next, label cells by genotype - samples 05 and 08 are WT, samples 06 and 07 are KO
geno_WT <- gsub("sample05|sample08", "WT", proj_Mef2c_v13_E9_subset$Sample)
geno <- gsub("sample06|sample07", "KO", geno_WT)
proj_Mef2c_v13_E9_subset$Genotype <- geno

#Next, create CellTypeByGenotype column for grouping
proj_Mef2c_v13_E9_subset <- addCellColData(ArchRProj = proj_Mef2c_v13_E9_subset, 
                                              data = paste0(proj_Mef2c_v13_E9_subset$Cluster_Labels_Subset_celltypegroup,"_x_",
                                                            proj_Mef2c_v13_E9_subset$Genotype), 
                                              name = "CellTypeByGenotype", 
                                              cells = getCellNames(proj_Mef2c_v13_E9_subset), 
                                              force = TRUE
                                              )

#Create CellTypeBySample column to allow easy viewing of number of cells per cluster
proj_Mef2c_v13_E9_subset <- addCellColData(ArchRProj = proj_Mef2c_v13_E9_subset, 
                                              data = paste0(proj_Mef2c_v13_E9_subset$Cluster_Labels_Subset_celltypegroup,"_x_",
                                                            proj_Mef2c_v13_E9_subset$Sample), 
                                              name = "CellTypeBySample", 
                                              cells = getCellNames(proj_Mef2c_v13_E9_subset), 
                                              force = TRUE
                                              )

###------------Add Group Coverages aka Create Pseudobulk Reps------------------

#Check cell count tables to inform pseudobulking params
table(proj_Mef2c_v13_E9_subset$CellTypeByGenotype)
table(proj_Mef2c_v13_E9_subset$CellTypeBySample)

#Now addGroupCoverages
proj_Mef2c_v13_E9_subset <- addGroupCoverages(
  ArchRProj = proj_Mef2c_v13_E9_subset,
  groupBy = "CellTypeByGenotype",
  useLabels = TRUE,
  sampleLabels = "Sample",
  minCells = 20,
  maxCells = 700,
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
proj_Mef2c_v13_E9_subset <- addReproduciblePeakSet(
  ArchRProj = proj_Mef2c_v13_E9_subset,
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

#                                 Group nCells nCellsUsed nReplicates nMin nMax maxPeaks
#C2_x_KO           C2_x_KO    149        149           2   52   97    74500
#C2_x_WT           C2_x_WT    150        150           2   73   77    75000
#CMs_IFT_x_KO CMs_IFT_x_KO    910        910           2  299  611   150000
#CMs_IFT_x_WT CMs_IFT_x_WT    477        477           2  164  313   150000
#CMs_OFT_x_KO CMs_OFT_x_KO    148        148           2   40  108    74000
#CMs_OFT_x_WT CMs_OFT_x_WT    365        365           2  160  205   150000
#CMs_V_x_KO     CMs_V_x_KO    415        415           2  162  253   150000
#CMs_V_x_WT     CMs_V_x_WT    614        614           2  235  379   150000
#Me_x_KO           Me_x_KO    324        324           2   77  247   150000
#Me_x_WT           Me_x_WT    408        408           2  149  259   150000
#Pe_x_KO           Pe_x_KO    363        363           2  157  206   150000
#Pe_x_WT           Pe_x_WT    173        173           2   86   87    86500
#SHF_x_KO         SHF_x_KO    867        867           2  327  540   150000
#SHF_x_WT         SHF_x_WT    907        907           2  316  591   150000
#VP_x_KO           VP_x_KO    561        561           2  215  346   150000
#VP_x_WT           VP_x_WT    468        468           2  183  285   150000

#addPeakMatrix
proj_Mef2c_v13_E9_subset <- addPeakMatrix(proj_Mef2c_v13_E9_subset)

#Check that PeakMatrix now appears in available matrices
getAvailableMatrices(proj_Mef2c_v13_E9_subset)


###------------Export Group Bigwigs (optional)-----------------------------------

getGroupBW(
  ArchRProj = proj_Mef2c_v13_E9_subset,
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
proj_Mef2c_v13_E9_subset <- addMotifAnnotations(ArchRProj = proj_Mef2c_v13_E9_subset, 
                                                motifSet = "vierstra", 
                                                collection = "archetype", 
                                                annoName = "Motif", 
                                                force = T
                                                )


###-------------Add Peak2Gene Linkages----------------------------------------
proj_Mef2c_v13_E9_subset <- addPeak2GeneLinks(
  ArchRProj = proj_Mef2c_v13_E9_subset,
  reducedDims = "Harmony_Combined_sub",
  useMatrix = "GeneExpressionMatrix"
)


###--------------Save and Load ArchR Project------------------------------------
saveArchRProject(ArchRProj = proj_Mef2c_v13_E9, overwrite = TRUE, dropCells = FALSE, load = TRUE)
saveArchRProject(ArchRProj = proj_Mef2c_v13_E9_subset, overwrite = TRUE, dropCells = FALSE, load = TRUE)

#to load
proj_Mef2c_v13_E9 <- loadArchRProject("~/Mef2c_ArchR_working/Mef2c_v13_E9_ATAC_and_GEX")
proj_Mef2c_v13_E9_subset <- loadArchRProject("~/Mef2c_ArchR_working/Mef2c_v13_E9_subset")
setwd("~/Mef2c_ArchR_working")


