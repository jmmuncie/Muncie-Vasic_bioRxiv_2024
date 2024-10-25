#This script is used to create ArchR Arrow Files for all 12 samples from the
#Mef2c multiome experiment (2 WT and 2 KO embryos at each timepoint E7.75, 
#E8.5, and E9)

#Script written by Jonathon M. Muncie-Vasic
#edited 10/25/2024


###---------------Load Libraries and Set Threads------------------------------

#To install ArchR release_1.0.2 branch:
devtools::install_github("GreenleafLab/ArchR", ref="release_1.0.2", repos = BiocManager::repositories())

#See https://www.archrproject.com/index.html for full ArchR installation instructions

library(ArchR); library(Seurat); library(dplyr)
library(GenomeInfoDb); library(ensembldb)
library(ggplot2); library(patchwork)
set.seed(1234)
library(BSgenome) 
library(hexbin)

#Set ArchR threads to 4
addArchRThreads(threads = 4, force = FALSE)


###---------Build Genome and Gene Annotations for Alignment----------------

#Builing GRanges object for mm10 blacklist
#mm10 blacklist1 regions source: https://github.com/Boyle-Lab/Blacklist/blob/master/lists/mm10-blacklist.v2.bed.gz
#mm10 blacklist2 regions source (mito genome peaks that map to nuclear genome): https://github.com/caleblareau/mitoblacklist/blob/master/peaks/mm10_peaks.narrowPeak
#Make sure blacklist files are in working directory
setwd("~/blacklist regions")
#Merging these two blacklists, as were used to generate the precompiled mm10 genome in ArchR
blacklist1 <- read.delim("mm10-blacklist.v2.bed", sep="\t", header = FALSE)
colnames(blacklist1) <- c("chr", "start", "end", "region")
black_gr1 <- makeGRangesFromDataFrame(blacklist1)
blacklist2 <- read.delim("mm10-blacklist.mito", sep="\t", header = FALSE)
colnames(blacklist2) <- c("chr", "start", "end", "region")
black_gr2 <- makeGRangesFromDataFrame(blacklist2)
black_gr3 <- c(black_gr1,black_gr2)

#Build genomeAnnotation from custom forged BSgenome and blacklist
library(BSgenome.Mmusculus.GENCODE.mm10v13) 
genomeAnnotation_mm10v13 = createGenomeAnnotation(genome = "BSgenome.Mmusculus.GENCODE.mm10v13", blacklist = black_gr3, filter = FALSE, filterChr = NULL)

#Build geneAnnotation from genes_v13.gtf file - which is the mm10 gtf file 
#I initially downloaded from 10x with the unlocalized scaffolds removed and the
#F6-eGFP transgene annotations added
#Make sure "genes_v13.gtf" is in working directory 
setwd("~/gtf_for_geneAnnotation")
anno_path <- ""
db <- ensDbFromGtf(gtf=paste0(anno_path,"genes_v13.gtf"), 
                   organism = "Mus_musculus", genomeVersion = "mm10", version = 1.0)
edb <- EnsDb(db); rm(db)

gene.ranges <- genes(edb); length(gene.ranges$gene_id)
tss.ranges <- resize(gene.ranges, 1, "start")
exon.ranges = exons(edb)

geneAnnotation_mm10v13 = createGeneAnnotation(
  TSS = tss.ranges, 
  exons = exon.ranges, 
  genes = gene.ranges
)
dbDisconnect()

save(genomeAnnotation_mm10v13, file = "~/Genome_Gene_Annotations/genomeAnnotation_mm10v13.Robj")
save(geneAnnotation_mm10v13, file = "~/Genome_Gene_Annotations/geneAnnotation_mm10v13.Robj")


###--------------Create Arrow Files and ArchR Project----------------------

#Point to input files
#These are the fragments.tsv.gz and fragments.tsv.gz.tbi files for each sample
#that are output from the Cellranger-arc counts pipeline

#Set working directory one level up from input files 
#e.g. if input files are stored in "~/Mef2c_ArchR_working/Mef2c_inputs", then:
setwd("~/Mef2c_ArchR_working")
#Note: getInputFiles looks for file names with pattern "fragments.tsv.gz"
inputFiles <- getInputFiles("Mef2c_inputs")

#Create Arrow Files
#Dont set filterTSS too high because you can always increase later, default is 4
#Set subThreading = FALSE after having error creating HDF5 file - recommended
#by similar error on ArchR github
ArrowFiles_Mef2c_v13 <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  minTSS = 4, 
  minFrags = 1000, 
  excludeChr = c("chrM", "chrY"),
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  genomeAnnotation = genomeAnnotation_mm10v13,
  geneAnnotation = geneAnnotation_mm10v13,
  subThreading = FALSE
)

#Inspect ArrowFiles - make sure they are saved to a folder and you can
#access them to create downstream ArchR projects
ArrowFiles_Mef2c_v13

#Check quality control plots generated during creation of Arrow Files before 
#proceeding to creating ArchR projects!



