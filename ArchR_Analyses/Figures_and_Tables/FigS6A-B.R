#Load libraries
library(ArchR); library(Seurat); library(dplyr)
library(GenomeInfoDb); library(ensembldb)
library(ggplot2); library(patchwork)
set.seed(1234)
library(BSgenome) 
library(hexbin)

#Load ArchR project
proj_Mef2c_v13_E85_subset <- loadArchRProject("~/Mef2c_ArchR_working/Mef2c_v13_E85_subset")

#Set working directory
setwd("~/Mef2c_ArchR_working")


#Get DARs
markerTest_CMs_IFT <- getMarkerFeatures(
  ArchRProj = proj_Mef2c_v13_E85_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "CellTypeByGenotype",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "CMs_IFT_x_KO",
  bgdGroups = "CMs_IFT_x_WT"
)

#DARs up in KO
markerTest_CMs_IFT_list_upinKO <- getMarkers(markerTest_CMs_IFT, cutOff = "FDR <= 0.15 & Log2FC >= 0.5", returnGR = TRUE)
#Let's convert the upinKO list to a GRanges object, sorted by Log2FC
markerpeaksgr_IFT_KO <- markerTest_CMs_IFT_list_upinKO$`CMs_IFT_x_KO`
order_IFT_KO <- order(markerpeaksgr_IFT_KO$Log2FC, decreasing = T, start(markerpeaksgr_IFT_KO))
markerpeaks_sorted_IFT_KO <- markerpeaksgr_IFT_KO[order_IFT_KO]


#KO DARs with NR Motifs - for Fig S6A

#Look for AC0422, AC0414, AC0416, and AC0434 motifs within peaks
#Rationale: these four nuclear receptor (likely Nr2f2) motifs show up in the Top 10 motifs enriched in KO IFT CM DARs
#First, use getMatches to find which peaks contain the specific motif annotation you are looking for
mp_AC0422 <- getMatches(ArchRProj = proj_Mef2c_v13_E85_subset,
                        name = "Motif",
                        annoName = "AC0422" #AC0422-ESRRA/NR5A/Nuclear_receptor
)
#Returns a sparse matrix with True/False logical for every peak in the peak set
#Next we want to subset only the matching peaks 
sub_mp_AC0422 <- mp_AC0422[which(mp_AC0422@assays@data@listData$matches), ]
#Then, find matches between the subset and our differential peaks 
match_KO_DARs_with_AC0422 <- findOverlaps(markerpeaks_sorted_IFT_KO, sub_mp_AC0422)
KO_DARs_with_AC0422 <- markerpeaks_sorted_IFT_KO[match_KO_DARs_with_AC0422@from]
#We see that 170 of the 890 differential peaks contain the AC0422 motif

#Repeat for the other motifs
mp_AC0414 <- getMatches(ArchRProj = proj_Mef2c_v13_E85_subset,
                        name = "Motif",
                        annoName = "AC0414" #AC0414-NR4A/NR2C/Nuclear_receptor
)
sub_mp_AC0414 <- mp_AC0414[which(mp_AC0414@assays@data@listData$matches), ]
match_KO_DARs_with_AC0414 <- findOverlaps(markerpeaks_sorted_IFT_KO, sub_mp_AC0414)
KO_DARs_with_AC0414 <- markerpeaks_sorted_IFT_KO[match_KO_DARs_with_AC0414@from]
#We see that 214 of the 890 differential peaks contain the AC0414 motif

mp_AC0416 <- getMatches(ArchRProj = proj_Mef2c_v13_E85_subset,
                        name = "Motif",
                        annoName = "AC0416" #AC0416-RXRA/NR2C/Nuclear_receptor
)
sub_mp_AC0416 <- mp_AC0416[which(mp_AC0416@assays@data@listData$matches), ]
match_KO_DARs_with_AC0416 <- findOverlaps(markerpeaks_sorted_IFT_KO, sub_mp_AC0416)
KO_DARs_with_AC0416 <- markerpeaks_sorted_IFT_KO[match_KO_DARs_with_AC0416@from]
#We see that 173 of the 890 differential peaks contain the AC0416 motif

mp_AC0434 <- getMatches(ArchRProj = proj_Mef2c_v13_E85_subset,
                        name = "Motif",
                        annoName = "AC0434" #AC0434-NR2F/RXRG/Nuclear_receptor
)
sub_mp_AC0434 <- mp_AC0434[which(mp_AC0434@assays@data@listData$matches), ]
match_KO_DARs_with_AC0434 <- findOverlaps(markerpeaks_sorted_IFT_KO, sub_mp_AC0434)
KO_DARs_with_AC0434 <- markerpeaks_sorted_IFT_KO[match_KO_DARs_with_AC0434@from]
#We see that 134 of the 890 differential peaks contain the AC0434 motif

#Combine the four Granges objects and remove duplicates to get peaks that contain at least one of the four NR motifs
KO_DARs_with_NRmotifs <- c(KO_DARs_with_AC0422, KO_DARs_with_AC0414, KO_DARs_with_AC0416, KO_DARs_with_AC0434)
order_KO_DARs_NR <- order(KO_DARs_with_NRmotifs$Log2FC, decreasing = T, start(KO_DARs_with_NRmotifs))
KO_DARs_with_NRmotifs <- KO_DARs_with_NRmotifs[order_KO_DARs_NR]
#Get indices for duplicated ranges
dup_indices <- which(duplicated(KO_DARs_with_NRmotifs))
#Get indices of unique entries (removes duplicated indices from the full sequence of indices)
unique_indices <- setdiff(seq_along(KO_DARs_with_NRmotifs), dup_indices)
#Subset Granges object using unique indices
KO_DARs_with_NRmotifs_unique <- KO_DARs_with_NRmotifs[unique_indices,]
#Note: this is 353 peaks 

#Intersect these peaks with the full peak set to get nearestGene and peakType info
pgr <- getPeakSet(proj_Mef2c_v13_E85_subset)
#Subset matching peaks
sub_pgr_KO_DARs_NR_indexes <- findOverlaps(KO_DARs_with_NRmotifs_unique, pgr)
sub_pgr_KO_DARs_NR <- pgr[sub_pgr_KO_DARs_NR_indexes@to]

#Now create data frame with DAR info that will become bed file
KO_DARs_with_NRmotifs_bed_df <- data.frame(
  chr = KO_DARs_with_NRmotifs_unique@seqnames,
  start = KO_DARs_with_NRmotifs_unique@ranges@start,
  end = KO_DARs_with_NRmotifs_unique@ranges@start + 500,
  mm10PeakRange = paste0(KO_DARs_with_NRmotifs_unique@seqnames,':',KO_DARs_with_NRmotifs_unique@ranges@start,'-',KO_DARs_with_NRmotifs_unique@ranges@start + 500),
  nearestGene = sub_pgr_KO_DARs_NR@elementMetadata$nearestGene, 
  peakType = sub_pgr_KO_DARs_NR@elementMetadata$peakType, 
  peakLog2FC = KO_DARs_with_NRmotifs_unique$Log2FC,
  peakFDR = KO_DARs_with_NRmotifs_unique$FDR
)
KO_DARs_with_NRmotifs_bed_df[, 2:3] <- lapply(KO_DARs_with_NRmotifs_bed_df[, 2:3], as.integer)

write.table(KO_DARs_with_NRmotifs_bed_df, "~/Mef2c_ArchR_working/FigS6A_E85_IFT-CMs_DARs_upinKO_with_NRmotifs_mm10.bed", sep="\t", row.names=FALSE, col.names=FALSE, quote = FALSE)


#KO DARs with GATA Motifs - for Fig S6B

#Look for AC0252 and AC0254 Gata motifs within peaks
#Rationale: these two motifs show up in the Top 10 motifs enriched in KO IFT CM DARs
#First, use getMatches to find which peaks contain the specific motif annotation you are looking for
mp_AC0252 <- getMatches(ArchRProj = proj_Mef2c_v13_E85_subset,
                        name = "Motif",
                        annoName = "AC0252" #AC0252-GATA/GATA1/TALIGATA
)
#Returns a sparse matrix with True/False logical for every peak in the peak set
#Next we want to subset only the matching peaks 
sub_mp_AC0252 <- mp_AC0252[which(mp_AC0252@assays@data@listData$matches), ]
#Then, find matches between the subset and our differential peaks 
match_KO_DARs_with_AC0252 <- findOverlaps(markerpeaks_sorted_IFT_KO, sub_mp_AC0252)
KO_DARs_with_AC0252 <- markerpeaks_sorted_IFT_KO[match_KO_DARs_with_AC0252@from]
#We see that 415 of the 890 differential peaks contain the AC0252 motif

#Repeat for the AC0254 motif
mp_AC0254 <- getMatches(ArchRProj = proj_Mef2c_v13_E85_subset,
                        name = "Motif",
                        annoName = "AC0254" #AC0254-GATA/TRPSGATA
)
sub_mp_AC0254 <- mp_AC0254[which(mp_AC0254@assays@data@listData$matches), ]
match_KO_DARs_with_AC0254 <- findOverlaps(markerpeaks_sorted_IFT_KO, sub_mp_AC0254)
KO_DARs_with_AC0254 <- markerpeaks_sorted_IFT_KO[match_KO_DARs_with_AC0254@from]
#We see that 173 of the 890 differential peaks contain the AC0254 motif

#Combine the two Granges objects and remove duplicates to get peaks that contain either motif
KO_DARs_with_GATAmotifs <- c(KO_DARs_with_AC0252, KO_DARs_with_AC0254)
order_KO_DARs_GATA <- order(KO_DARs_with_GATAmotifs$Log2FC, decreasing = T, start(KO_DARs_with_GATAmotifs))
KO_DARs_with_GATAmotifs <- KO_DARs_with_GATAmotifs[order_KO_DARs_GATA]
#Get indices for duplicated ranges
dup_indices <- which(duplicated(KO_DARs_with_GATAmotifs))
#Get indices of unique entries (removes duplicated indices from the full sequence of indices)
unique_indices <- setdiff(seq_along(KO_DARs_with_GATAmotifs), dup_indices)
#Subset Granges object using unique indices
KO_DARs_with_GATAmotifs_unique <- KO_DARs_with_GATAmotifs[unique_indices,]
#Note: this is 415 peaks - so every peak that had the AC0254 motif annotated, was contained in the list of peaks with AC0252 annotated

#Intersect these peaks with the full peak set to get nearestGene and peakType info
pgr <- getPeakSet(proj_Mef2c_v13_E85_subset)
#Subset matching peaks
sub_pgr_KO_DARs_GATA_indexes <- findOverlaps(KO_DARs_with_GATAmotifs_unique, pgr)
sub_pgr_KO_DARs_GATA <- pgr[sub_pgr_KO_DARs_GATA_indexes@to]

#Now create data frame with DAR info that will become bed file
KO_DARs_with_GATAmotifs_bed_df <- data.frame(
  chr = KO_DARs_with_GATAmotifs_unique@seqnames,
  start = KO_DARs_with_GATAmotifs_unique@ranges@start,
  end = KO_DARs_with_GATAmotifs_unique@ranges@start + 500,
  mm10PeakRange = paste0(KO_DARs_with_GATAmotifs_unique@seqnames,':',KO_DARs_with_GATAmotifs_unique@ranges@start,'-',KO_DARs_with_GATAmotifs_unique@ranges@start + 500),
  nearestGene = sub_pgr_KO_DARs_GATA@elementMetadata$nearestGene, 
  peakType = sub_pgr_KO_DARs_GATA@elementMetadata$peakType, 
  peakLog2FC = KO_DARs_with_GATAmotifs_unique$Log2FC,
  peakFDR = KO_DARs_with_GATAmotifs_unique$FDR
)
KO_DARs_with_GATAmotifs_bed_df[, 2:3] <- lapply(KO_DARs_with_GATAmotifs_bed_df[, 2:3], as.integer)

write.table(KO_DARs_with_GATAmotifs_bed_df, "~/Mef2c_ArchR_working/FigS6B_E85_IFT-CMs_DARs_upinKO_with_GATAmotifs_mm10.bed", sep="\t", row.names=FALSE, col.names=FALSE, quote = FALSE)


