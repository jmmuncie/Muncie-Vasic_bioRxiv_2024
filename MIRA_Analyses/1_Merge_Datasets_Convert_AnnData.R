#This script is used to generate a single Seurat object for Mef2c Multiome 
#data, containing cells of interest from E7.75, E8.5, and E9 timepoints.
#The input data was processed via a separate set of scripts - already QC'd, 
#subset, etc. So here we are mostly just further subsetting and merging objects

#Before doing this the first time, you will need to create the "SCEtoAD" conda
#environment that will be used by sceasy. To do that, install miniconda
#(https://docs.conda.io/en/latest/miniconda.html)
#and then run the following line of code in the terminal:

#conda create -n SCEtoAD -c bioconda anndata


#Load libraries
library(Seurat)
library(SingleCellExperiment)
library(LoomExperiment)
library(anndata)
library(sceasy)
library(reticulate)
use_condaenv("SCEtoAD")
library(ArchR); library(Seurat); library(dplyr)
library(GenomeInfoDb); library(ensembldb)
library(ggplot2); library(patchwork)
set.seed(1234)
library(BSgenome) 
library(hexbin)

#Load in starting Seurat objects
setwd("~/Desktop/MIRA_analysis/Mef2c")
load("Seurat_Objects_Inputs/mef2c_v13_E775_subset_scTrans.Robj")
load("Seurat_Objects_Inputs/mef2c_v13_E85_scTrans_v3_harmony_subset.Robj")
load("Seurat_Objects_Inputs/mef2c_v13_E9_subset_scTrans.Robj")

#Examine UMAPs to check cell types/clustering
#E7.75
Idents(mef2c_v13_E775_subset) <- "cell_type_subset"
DimPlot(mef2c_v13_E775_subset, 
        reduction = "umap", 
        label = T, 
        repel = T, 
        label.size = 4, 
        pt.size = 0.5
        ) + NoLegend()
#E8.5
Idents(mef2c_v13_E85_v3_harmony_subset) <- "harmony_cell_type_subset"
DimPlot(mef2c_v13_E85_v3_harmony_subset, 
        reduction = "umap", 
        label = T, 
        repel = T, 
        label.size = 4, 
        pt.size = 0.5
        ) + NoLegend()
#E9
Idents(mef2c_v13_E9_subset) <- "cell_type_subset"
DimPlot(mef2c_v13_E9_subset, 
        reduction = "umap", 
        label = T, 
        repel = T, 
        label.size = 4, 
        pt.size = 0.5
        ) + NoLegend()


#Subset cell types of interest for MIRA
#E7.75 - brings 516 cells
Idents(mef2c_v13_E775_subset) <- "cell_type_subset"
mef2c_E775_MIRA_subset <- subset(mef2c_v13_E775_subset, 
                                 idents = c("CMs/FHF", "JCF", "SHF1", "SHF2"
                                            )
                                 )
#E8.5 - brings 4369 cells
Idents(mef2c_v13_E85_v3_harmony_subset) <- "harmony_cell_type_subset"
mef2c_E85_MIRA_subset <- subset(mef2c_v13_E85_v3_harmony_subset, 
                                 idents = c("CMs-V", "CMs-OFT", "CMs-IFT1",
                                            "CMs-IFT2", "CMs-IFT3", "aSHF",
                                            "pSHF"
                                            )
                                )
#E9 - brings 4364 cells
#Note: in this iteration, excluding VP, Pe, and SHF from E9 - not so interested
#in these cell types at E9, at this stage, higher priority to simplify for MIRA
Idents(mef2c_v13_E9_subset) <- "cell_type_subset"
mef2c_E9_MIRA_subset <- subset(mef2c_v13_E9_subset, 
                                idents = c("CMs-V", "CMs-OFT", "CMs-AVC",
                                           "CMs-A1", "CMs-A2"
                                           )
                                )
  
#Now merge the three objects - 9,249 cells
#Note that by default, merge will combine Seurat objects based on the raw 
#count matrices, erasing any previously normalized and scaled data matrices.
mef2c_MIRA_initial <- merge(mef2c_E775_MIRA_subset, 
                            y = c(mef2c_E85_MIRA_subset, mef2c_E9_MIRA_subset), 
                            project = "mef2c_MIRA"
                            )
mef2c_MIRA_initial
#But it looks like this object does still contain the SCT assay, so when we 
#convert to AnnData we will need to make sure we bring over the raw counts. 
#Let's set the default assay to RNA to help avoid any issues. 
DefaultAssay(mef2c_MIRA_initial) <- "RNA"

#Note that the orig.idents brought over from the merged objects list both the 
#timepoint and the sample number from the original objects. However, when I 
#initially created the E9 object (samples05-08), I mislabeled them as "E775" -
#let's correct those here
table(mef2c_MIRA_initial$orig.ident)
class.vec1 <- gsub("E775_sample05","E9_sample05", (mef2c_MIRA_initial$orig.ident))
class.vec2 <- gsub("E775_sample06","E9_sample06", class.vec1)
class.vec3 <- gsub("E775_sample07","E9_sample07", class.vec2)
class.vec4 <- gsub("E775_sample08","E9_sample08", class.vec3)
mef2c_MIRA_initial$orig.ident <- class.vec4
table(mef2c_MIRA_initial$orig.ident)
#Now the orig.ident labels are correct. We should make sure to bring these over
#to our AnnData object so we can use them to label timepoint/sample/genotype
#as needed in MIRA

#Save R object 
save(mef2c_MIRA_initial, file = "~/Desktop/MIRA_analysis/Mef2c/Mef2c_Seurat_Merged_initial.Robj")


###--------------Create ArchR Project----------------------

#Arrow Files previously created in "Mef2c_v13_ArchR_create_arrows.R"
#Copied the original arrow files over into path indicated below to generate
#the inital project - deleted these after creating the overlapped project 
#(intersected with the Seurat object) to save space

#Load genome and gene annotations
load("~/Desktop/Mef2c_v13_ArchR_working/Genome_Gene_Annotations/genomeAnnotation_mm10v13.Robj")
load("~/Desktop/Mef2c_v13_ArchR_working/Genome_Gene_Annotations/geneAnnotation_mm10v13.Robj")

#Create ArchR Project
setwd("~/Desktop/MIRA_analysis/Mef2c")
proj_Mef2c_MIRA_initial <- ArchRProject(
  ArrowFiles = c("Original_Arrows/sample01.arrow", 
                 "Original_Arrows/sample02.arrow", 
                 "Original_Arrows/sample03.arrow", 
                 "Original_Arrows/sample04.arrow",
                 "Original_Arrows/sample05.arrow", 
                 "Original_Arrows/sample06.arrow", 
                 "Original_Arrows/sample07.arrow", 
                 "Original_Arrows/sample08.arrow",
                 "Original_Arrows/sample21.arrow", 
                 "Original_Arrows/sample22.arrow", 
                 "Original_Arrows/sample23.arrow", 
                 "Original_Arrows/sample24.arrow"
  ),
  outputDirectory = "Mef2c_MIRA_ArchR_proj_initial",
  copyArrows = FALSE,
  geneAnnotation = geneAnnotation_mm10v13,
  genomeAnnotation = genomeAnnotation_mm10v13,
  showLogo = FALSE,
  threads = 4
)

#Inspect ArchR Project
proj_Mef2c_MIRA_initial
#View which data matrices are available in the project
getAvailableMatrices(proj_Mef2c_MIRA_initial)


###----------Find overlaps with cells in Seurat object------------------

seurat_obj_names <- gsub("_", "#", mef2c_MIRA_initial@assays$RNA@data@Dimnames[[2]])
overlap <- intersect(proj_Mef2c_MIRA_initial$cellNames, seurat_obj_names)
#There are 8,640 cells that exist in both the ArchR project and the Seurat object

#IMPORTANT: Need to install and use "dev_emptyChr" branch 
devtools::install_github("GreenleafLab/ArchR", ref="dev_emptyChr", repos = BiocManager::repositories(), force =TRUE)
detach("package:ArchR", unload = TRUE)
library(ArchR)

proj_Mef2c_MIRA_final <- subsetArchRProject(
  ArchRProj = proj_Mef2c_MIRA_initial,
  cells = overlap,
  outputDirectory = "Mef2c_MIRA_ArchR_proj_final",
  dropCells = TRUE,
  force = TRUE
)

#IMPORTANT: Re-install "release_1.0.2" branch before proceeding
devtools::install_github("GreenleafLab/ArchR", ref="release_1.0.2", repos = BiocManager::repositories())
detach("package:ArchR", unload = TRUE)
library(ArchR)

proj_Mef2c_MIRA_final


###---------QC and additional TSSEnrichment and nFrag Filtering----------------

#Based on QC metrics, I think I want to further refine cutoff for TSSEnrichment
#Let's build the QC plot TSSEnrichment vs. log10 Unique Frags now that cells have
#been filtered
#Get nFrags and TSSEnrichment scores
df <- getCellColData(proj_Mef2c_MIRA_final, select = c("log10(nFrags)", "TSSEnrichment"))
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
#View plot
p
#Save plot
ggsave("Mef2c_MIRA_TSSbyLog10nFrags_initial.pdf", plot = p, device = "pdf", path = "~/Desktop/MIRA_analysis/Mef2c/ArchR_QC_Plots", width = 5, height = 5, dpi = 300)

#Before filtering, visualize subset of cells with new cutoffs
#TSS Enrichiment >= 12, TSS Enrichiment <= 30 nFrags >= 5011 (log10 3.7)
idxPass <- which(proj_Mef2c_MIRA_final$TSSEnrichment >= 12 & proj_Mef2c_MIRA_final$TSSEnrichment <= 30 & proj_Mef2c_MIRA_final$nFrags >= 5011)
cellsPass <- proj_Mef2c_MIRA_final$cellNames[idxPass]
df2 <- getCellColData(proj_Mef2c_MIRA_final[cellsPass, ], select = c("log10(nFrags)", "TSSEnrichment"))
p2 <- ggPoint(
  x = df2[,1], 
  y = df2[,2], 
  colorDensity = TRUE,
  continuousSet = "sambaNight",
  xlabel = "Log10 Unique Fragments",
  ylabel = "TSS Enrichment",
  xlim = c(log10(500), quantile(df2[,1], probs = 1.0)+0.5),
  ylim = c(0, quantile(df2[,2], probs = 1.0)+5)
) + geom_hline(yintercept = 12, lty = "dashed") + geom_hline(yintercept = 30, lty = "dashed") + geom_vline(xintercept = 3.7, lty = "dashed")
#View plot
p2
#Save plot
ggsave("Mef2c_MIRA_TSSbyLog10nFrags_filtered.pdf", plot = p2, device = "pdf", path = "~/Desktop/MIRA_analysis/Mef2c/ArchR_QC_Plots", width = 5, height = 5, dpi = 300)

#Looks good, now filter these cells from project 
proj_Mef2c_MIRA_final <- proj_Mef2c_MIRA_final[cellsPass, ]
proj_Mef2c_MIRA_final
#Now at 7,267 cells


###----------------Doublet Filtering-----------------------------------------
#NOTE: WAS NOT ABLE TO FILTER DOUBLETS DUE TO TOO FEW CELLS IN SAMPLE 04

#First add doublet scores
#dir.create("~/Desktop/MIRA_analysis/Mef2c/Mef2c_MIRA_ArchR_proj_final/DoubletSummaries")
#proj_Mef2c_MIRA_final <- addDoubletScores(proj_Mef2c_MIRA_final, outDir = paste(getOutputDirectory(proj_Mef2c_MIRA_final), "DoubletSummaries", sep = "/"))
#Note: gave error at sample 04, I think because too few cells, not sure if will be able to filter doublets

#Save separate instance of project prior to filtering doublets
#Will delete this to save disk space after confirming downstream analyses look good
#proj_Mef2c_MIRA_ArchR_proj_final_before_doub_filter <- proj_Mef2c_MIRA_final
#saveArchRProject(
#  ArchRProj = proj_Mef2c_MIRA_ArchR_proj_final_before_doub_filter, 
#  outputDirectory = "Mef2c_MIRA_ArchR_proj_final_before_doub_filter",
#  overwrite = TRUE, 
#  dropCells = F, 
#  load = TRUE
#)

#If you choose to filter, run the following
#proj_Mef2c_MIRA_final <- filterDoublets(proj_Mef2c_MIRA_final, filterRatio = 1)

#Filtering ___ cells from ArchRProject!

#Now at ____ cells
#saveArchRProject(ArchRProj = proj_Mef2c_MIRA_final, overwrite = TRUE, dropCells = TRUE, load = TRUE)


###--------------Dim Reduction - LSI ATAC -----------------------------------
proj_Mef2c_MIRA_final <- addIterativeLSI(
  ArchRProj = proj_Mef2c_MIRA_final, 
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


###-------------UMAP--------------------------------------------------------
#ATAC UMAP
proj_Mef2c_MIRA_final <- addUMAP(
  ArchRProj = proj_Mef2c_MIRA_final, 
  reducedDims = "LSI_ATAC", 
  name = "UMAP_ATAC", 
  minDist = 0.5, 
  force = TRUE
)
#Clusters
proj_Mef2c_MIRA_final <- addClusters(
  input = proj_Mef2c_MIRA_final, 
  reducedDims = "LSI_ATAC", 
  maxClusters = 30,
  name = "Clusters_res0.5", 
  resolution = 0.5, 
  force = TRUE
)
#Plot Embedding
u1 <- plotEmbedding(
  ArchRProj = proj_Mef2c_MIRA_final, 
  name = "Clusters_res0.5", 
  embedding = "UMAP_ATAC", 
  size = 0.5, 
  labelAsFactors=F, 
  labelMeans=F,
  colorTitle = "",
  baseSize = 16,
) + theme(legend.text = element_text(size=16)) + guides(colour = guide_legend(override.aes = list(size=3)))
#Save
plotPDF(u1, name = "UMAPs_ATAC_res0.5", width = 10, height = 10, ArchRProj = proj_Mef2c_MIRA_final, addDOC = F)

#Label by timepoint/genotype and check UMAP
sn1 <- gsub("sample01", "E7.75_WT", proj_Mef2c_MIRA_final$Sample)
sn2 <- gsub("sample02", "E7.75_KO", sn1)
sn3 <- gsub("sample03", "E7.75_KO", sn2)
sn4 <- gsub("sample04", "E7.75_WT", sn3)
sn5 <- gsub("sample05", "E9_WT", sn4)
sn6 <- gsub("sample06", "E9_KO", sn5)
sn7 <- gsub("sample07", "E9_KO", sn6)
sn8 <- gsub("sample08", "E9_WT", sn7)
sn9 <- gsub("sample21", "E8.5_KO", sn8)
sn10 <- gsub("sample22", "E8.5_WT", sn9)
sn11 <- gsub("sample23", "E8.5_WT", sn10)
sn12 <- gsub("sample24", "E8.5_KO", sn11)
proj_Mef2c_MIRA_final$SampleNames <- sn12
#Plot Embedding
u2 <- plotEmbedding(
  ArchRProj = proj_Mef2c_MIRA_final, 
  name = "SampleNames", 
  embedding = "UMAP_ATAC", 
  size = 0.5, 
  labelAsFactors=F, 
  labelMeans=F,
  colorTitle = "",
  baseSize = 16,
) + theme(legend.text = element_text(size=16)) + guides(colour = guide_legend(override.aes = list(size=3)))
#Save
plotPDF(u2, name = "UMAPs_ATAC_bySample", width = 10, height = 10, ArchRProj = proj_Mef2c_MIRA_final, addDOC = F)

#Next, create ClusterBySampleNames (timepoint/genotype) column for grouping
proj_Mef2c_MIRA_final <- addCellColData(ArchRProj = proj_Mef2c_MIRA_final, 
                                        data = paste0(proj_Mef2c_MIRA_final$Clusters_res0.5,"_x_",
                                                      proj_Mef2c_MIRA_final$SampleNames), 
                                        name = "ClusterBySampleNames", 
                                        cells = getCellNames(proj_Mef2c_MIRA_final), 
                                        force = TRUE
)
#Create ClusterBySample column to allow easy viewing of number of cells per cluster
proj_Mef2c_MIRA_final <- addCellColData(ArchRProj = proj_Mef2c_MIRA_final, 
                                        data = paste0(proj_Mef2c_MIRA_final$Clusters_res0.5,"_x_",
                                                      proj_Mef2c_MIRA_final$Sample), 
                                        name = "ClusterBySample", 
                                        cells = getCellNames(proj_Mef2c_MIRA_final), 
                                        force = TRUE
)
table(proj_Mef2c_MIRA_final$ClusterBySampleNames)
table(proj_Mef2c_MIRA_final$ClusterBySample)


#IMPORTANT: Need to install and use "dev_emptyChr" branch 
devtools::install_github("GreenleafLab/ArchR", ref="dev_emptyChr", repos = BiocManager::repositories())
detach("package:ArchR", unload = TRUE)
library(ArchR)


#Now addGroupCoverages
#Took about 40 min for subset proj of 10,500 cells and fairly large number of cluster/genotype groups
proj_Mef2c_MIRA_final <- addGroupCoverages(
  ArchRProj = proj_Mef2c_MIRA_final,
  groupBy = "ClusterBySampleNames",
  useLabels = TRUE,
  sampleLabels = "Sample",
  minCells = 20,
  maxCells = 1000,
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
  excludeChr = c("chrM", "chrY"),
  logFile = createLogFile("addGroupCoverages")
)


#IMPORTANT: Re-install "release_1.0.2" branch before proceeding
devtools::install_github("GreenleafLab/ArchR", ref="release_1.0.2", repos = BiocManager::repositories())
detach("package:ArchR", unload = TRUE)
library(ArchR)


#addReproduciblePeakSet
#Took about 30 minutes for subset of 10,500 cells

proj_Mef2c_MIRA_final <- addReproduciblePeakSet(
  ArchRProj = proj_Mef2c_MIRA_final,
  groupBy = "ClusterBySampleNames",
  #Note: with reproducibility = 2, both my pseodobulk reps for each group will
  #be required to contain a peak call at the locus - this will give a strict
  #first pass at the analysis 
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
#                        Group nCells nCellsUsed nReplicates nMin nMax maxPeaks
#C1_x_E7.75_KO C1_x_E7.75_KO    244        244           2   60  184   122000
#C1_x_E7.75_WT C1_x_E7.75_WT    158        158           2   43  115    79000
#C1_x_E8.5_KO   C1_x_E8.5_KO      1          1           2    1    1      500
#C1_x_E8.5_WT   C1_x_E8.5_WT      2          2           2    2    2     1000
#C2_x_E8.5_KO   C2_x_E8.5_KO    369        369           2   73  296   150000
#C2_x_E8.5_WT   C2_x_E8.5_WT    431        431           2  196  235   150000
#C2_x_E9_KO       C2_x_E9_KO      4          4           2    3    4     2000
#C2_x_E9_WT       C2_x_E9_WT      7          7           2    5    7     3500
#C3_x_E7.75_WT C3_x_E7.75_WT      1          1           2    1    1      500
#C3_x_E8.5_KO   C3_x_E8.5_KO     37         33           2   20   20    16500
#C3_x_E8.5_WT   C3_x_E8.5_WT     52         52           2   20   32    26000
#C3_x_E9_KO       C3_x_E9_KO    206        206           2   45  161   103000
#C3_x_E9_WT       C3_x_E9_WT    175        175           2   70  105    87500
#C4_x_E8.5_KO   C4_x_E8.5_KO    154        154           2   36  118    77000
#C4_x_E8.5_WT   C4_x_E8.5_WT    335        335           2  148  187   150000
#C4_x_E9_KO       C4_x_E9_KO      3          3           2    2    3     1500
#C4_x_E9_WT       C4_x_E9_WT     16         14           2    9    9     7000
#C5_x_E7.75_WT C5_x_E7.75_WT      2          2           2    2    2     1000
#C5_x_E8.5_KO   C5_x_E8.5_KO    180        180           2   77  103    90000
#C5_x_E8.5_WT   C5_x_E8.5_WT    180        180           2   73  107    90000
#C5_x_E9_KO       C5_x_E9_KO      3          3           2    2    3     1500
#C5_x_E9_WT       C5_x_E9_WT     16         14           2    9    9     7000
#C6_x_E8.5_KO   C6_x_E8.5_KO    207        207           2   51  156   103500
#C6_x_E8.5_WT   C6_x_E8.5_WT    259        259           2  120  139   129500
#C6_x_E9_WT       C6_x_E9_WT     29         29           2   20   20    14500
#C7_x_E8.5_KO   C7_x_E8.5_KO     10         10           2    7    8     5000
#C7_x_E8.5_WT   C7_x_E8.5_WT      7          7           2    5    7     3500
#C7_x_E9_KO       C7_x_E9_KO    650        650           2  192  458   150000
#C7_x_E9_WT       C7_x_E9_WT    205        205           2   91  114   102500
#C8_x_E8.5_KO   C8_x_E8.5_KO    554        554           2  161  393   150000
#C8_x_E8.5_WT   C8_x_E8.5_WT      8          8           2    6    6     4000
#C8_x_E9_KO       C8_x_E9_KO      3          3           2    2    3     1500
#C8_x_E9_WT       C8_x_E9_WT      1          1           2    1    1      500
#C9_x_E7.75_WT C9_x_E7.75_WT      2          2           2    2    2     1000
#C9_x_E8.5_WT   C9_x_E8.5_WT    279        279           2  135  144   139500
#C9_x_E9_KO       C9_x_E9_KO      1          1           2    1    1      500
#C9_x_E9_WT       C9_x_E9_WT    295        295           2  126  169   147500
#C10_x_E8.5_KO C10_x_E8.5_KO     89         89           2   20   69    44500
#C10_x_E8.5_WT C10_x_E8.5_WT      1          1           2    1    1      500
#C10_x_E9_KO     C10_x_E9_KO     86         86           2   27   59    43000
#C10_x_E9_WT     C10_x_E9_WT      2          2           2    2    2     1000
#C11_x_E8.5_KO C11_x_E8.5_KO      2          2           2    2    2     1000
#C11_x_E9_KO     C11_x_E9_KO    667        667           2  236  431   150000
#C11_x_E9_WT     C11_x_E9_WT     31         30           2   20   20    15000
#C12_x_E8.5_WT C12_x_E8.5_WT    511        511           2  233  278   150000
#C12_x_E9_KO     C12_x_E9_KO      4          4           2    3    4     2000
#C12_x_E9_WT     C12_x_E9_WT    514        514           2  204  310   150000
#C13_x_E8.5_WT C13_x_E8.5_WT      3          3           2    2    3     1500
#C13_x_E9_KO     C13_x_E9_KO     13         13           2    9    9     6500
#C13_x_E9_WT     C13_x_E9_WT    258        258           2   93  165   129000


#addPeakMatrix
proj_Mef2c_MIRA_final <- addPeakMatrix(proj_Mef2c_MIRA_final)
getAvailableMatrices(proj_Mef2c_MIRA_final)

#getPeakMatrix
peaks <- getMatrixFromProject(proj_Mef2c_MIRA_final, useMatrix = "PeakMatrix")

#Convert to SingleCellExperiment format
mef2c_MIRA_peaks <- as(peaks, "SingleCellExperiment")

#Replace Peak idxs (1, 2, 3, ...) with actual ranges
seq_names <- as.vector(peaks@rowRanges@seqnames)
peak_ranges <- mapply(paste0, seq_names, ":", start(peaks@rowRanges@ranges), "-", end(peaks@rowRanges@ranges))
names(peak_ranges) <- peak_ranges
mef2c_MIRA_peaks@rowRanges$idx <- peak_ranges

#Use sceasy to convert from SCE to AnnData
setwd("~/Desktop/MIRA_analysis/Mef2c/")
sceasy::convertFormat(mef2c_MIRA_peaks, from="sce", to="anndata",
                      outFile="Mef2c_Peaks_AnnData.h5ad")

#Save ArchRProj
saveArchRProject(ArchRProj = proj_Mef2c_MIRA_final, overwrite = TRUE, dropCells = FALSE, load = TRUE)


###-------------------FINAL OBJECT-------------------------------------------

#After creating and doing QC/DoubFilt on the full ArchR project, find overlaps between
#the cells in the ArchR proj and the initial Seurat obj to filter out cells
#that did not pass the ATAC QC/DoubFilt - after this, ArchR proj and Seurat obj should
#contain exactly the same set of cell barcodes

#Note that by examining number of cells per sample in the Seurat object compared
#to those remaining in the ArchR project after QC, we can confirm that
#no sample(s) were overwhelmingly biased by QC filters
table(mef2c_MIRA_initial$orig.ident)
table(proj_Mef2c_MIRA_final$Sample)
#              inArchR     inSeurat     %filtered
# sample01        119           135          11.8
# sample02        184           258          28.6
# sample03         60            69          13.0
# sample04         44            54          18.5
# sample05        657           894          26.5
# sample06        514           741          30.6
# sample07       1128          1662          32.1
# sample08        892          1067          16.4
# sample21       1178          1414          16.6
# sample22        922          1069          13.7
# sample23       1146          1395          17.8
# sample24        424           491          13.6
# Total          7267          9249          21.4

#Find overlapping cells and subset
archr_names <- gsub("#", "_", proj_Mef2c_MIRA_final$cellNames)
overlap <- intersect(archr_names, mef2c_MIRA_initial@assays$RNA@data@Dimnames[[2]])
mef2c_MIRA_final <- subset(mef2c_MIRA_initial, cells = overlap)
#Seurat object now contains 7267 cells

#Clean up and create new metadata info
unique(colnames(mef2c_MIRA_final@meta.data))
mef2c_MIRA_final$SCT_snn_res.0.4 <- NULL
mef2c_MIRA_final$SCT_snn_res.0.5 <- NULL
mef2c_MIRA_final$SCT_snn_res.0.6 <- NULL
mef2c_MIRA_final$SCT_snn_res.0.7 <- NULL
mef2c_MIRA_final$SCT_snn_res.1 <- NULL
mef2c_MIRA_final$SCT_snn_res.1.2 <- NULL
mef2c_MIRA_final$SCT_snn_res.1.5 <- NULL
mef2c_MIRA_final$SCT_snn_res.2 <- NULL
mef2c_MIRA_final$seurat_clusters <- NULL
mef2c_MIRA_final$cell_type <- NULL
mef2c_MIRA_final$cell_type_subset <- NULL
mef2c_MIRA_final$harmony_cell_type <-NULL
mef2c_MIRA_final$harmony_cell_type_subset <- NULL
mef2c_MIRA_final$cell_type_pool_x_genotype <- NULL

#combine pooled cell type labels
#Get vector of labels from "harmony_cell_type_subset_pool" 
#These come from E8.5 dataset
harm_labels <- mef2c_MIRA_final$harmony_cell_type_subset_pool
#Get vector of labels from "cell_type_subset_pool"
#These come from E7.75 and E9 datasets
cell_labels <- mef2c_MIRA_final$cell_type_subset_pooled
#ID the na values in cell_labels
na_idx <- is.na(cell_labels)
#Replace those with the harm_labels, i.e. merge in the E8.5 cell type labels
cell_labels[na_idx] <- harm_labels[na_idx]
mef2c_MIRA_final$Cell_Type <- cell_labels
#Remove now redundant "harmony_cell_type_subset_pool" and "cell_type_subset_pooled" metadata
mef2c_MIRA_final$harmony_cell_type_subset_pool <- NULL
mef2c_MIRA_final$cell_type_subset_pooled <- NULL

#Create Timepoint label
time <- mef2c_MIRA_final$Sample_Name
time[time=="E7.75_WT1" | time=="E7.75_WT2" | time=="E7.75_KO1" | time=="E7.75_KO2"] <- "E7.75"
time[time=="E8.5_WT1" | time=="E8.5_WT2" | time=="E8.5_KO1" | time=="E8.5_KO2"] <- "E8.5"
time[time=="E9.0_WT1" | time=="E9.0_WT2" | time=="E9.0_KO1" | time=="E9.0_KO2"] <- "E9"
mef2c_MIRA_final$Timepoint <- time

#Create combined Timepoint_Cell_Type label
newlabel <- paste0(mef2c_MIRA_final$Timepoint,"_", mef2c_MIRA_final$Cell_Type)
mef2c_MIRA_final$Timepoint_x_Cell_Type <- newlabel

#Save Seurat object
save(mef2c_MIRA_final, file = "~/Desktop/MIRA_analysis/Mef2c/Mef2c_Seurat_Merged_final.Robj")

#Use sceasy to convert from Seurat to AnnData
setwd("~/Desktop/MIRA_analysis/Mef2c/")
sceasy::convertFormat(mef2c_MIRA_final, from="seurat", to="anndata",
                      outFile="Mef2c_RNA_AnnData.h5ad")
