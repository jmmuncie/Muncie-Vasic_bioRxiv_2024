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


#Get Marker Peaks
markerTest_CMs_IFT <- getMarkerFeatures(
  ArchRProj = proj_Mef2c_v13_E85_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "CellTypeByGenotype",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "CMs_IFT_x_KO",
  bgdGroups = "CMs_IFT_x_WT"
)

markerTest_CMs_IFT_list_upinKO <- getMarkers(markerTest_CMs_IFT, cutOff = "FDR <= 0.15 & Log2FC >= 0.5", returnGR = TRUE)
markerTest_CMs_IFT_list_upinWT <- getMarkers(markerTest_CMs_IFT, cutOff = "FDR <= 0.15 & Log2FC <= -0.5", returnGR = TRUE)

#Let's convert the upinWT list to a GRanges object, sorted by Log2FC
markerpeaksgr_IFT_WT <- markerTest_CMs_IFT_list_upinWT$`CMs_IFT_x_KO`
order_IFT_WT <- order(markerpeaksgr_IFT_WT$Log2FC, start(markerpeaksgr_IFT_WT))
markerpeaks_sorted_IFT_WT <- markerpeaksgr_IFT_WT[order_IFT_WT]

#Let's convert the upinKO list to a GRanges object, sorted by Log2FC
markerpeaksgr_IFT_KO <- markerTest_CMs_IFT_list_upinKO$`CMs_IFT_x_KO`
order_IFT_KO <- order(markerpeaksgr_IFT_KO$Log2FC, decreasing = T, start(markerpeaksgr_IFT_KO))
markerpeaks_sorted_IFT_KO <- markerpeaksgr_IFT_KO[order_IFT_KO]


#Retrieve P2G links
p2g_df <- getPeak2GeneLinks(
  ArchRProj = proj_Mef2c_v13_E85_subset,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)
#Get peak ids from p2g_df
peakids <- p2g_df$idxATAC
#Get gene ids from p2g_df
geneids <- p2g_df$idxRNA

#Get the GR object that consists of all peaks 
gr_peaks <- metadata(p2g_df)[[1]]
#Subset only the peaks that are linked to genes
linked_peaks <- gr_peaks[peakids]
#Add idxATAC and idxRNA to linked_peaks GRanges object to keep track 
linked_peaks$idxATAC <- peakids
linked_peaks$idxRNA <- geneids


#Now we want to find matches between our differential peaks up in KO
#and the linked peaks:
#'findOverlaps()' will identify matching Granges 
match_KO <- findOverlaps(markerpeaks_sorted_IFT_KO, linked_peaks)
#Suset both lists for matches
markerpeaks_sorted_IFT_KO_matched <- markerpeaks_sorted_IFT_KO[match_KO@from]
linkedpeaks_matched_IFT_KO <- linked_peaks[match_KO@to]
#Note that these two granges above match, only that the markerpeaks object 
#contains the metadata including log2fc from the differential peak analysis 

#Finally, we need to go back and pull the linked differential peaks from the 
#linked peaks data frame, which contains metadata for the p2g links that we 
#will want to use
p2g_diff_df_KO <- p2g_df[match_KO@to,]

#Now let's add all the relevant metadata to the linkedpeaks_matched GRanges object

#First let's add the gene names
linkedpeaks_matched_IFT_KO$GeneName <- metadata(p2g_df)$geneSet$name[linkedpeaks_matched_IFT_KO$idxRNA]
#Next, add Log2FC and FDR values from the diff peak testing
peaktest_KO <- markerpeaks_sorted_IFT_KO[match_KO@from, ]
linkedpeaks_matched_IFT_KO$Log2FC_DiffPeakTest <- peaktest_KO$Log2FC
linkedpeaks_matched_IFT_KO$FDR_DiffPeakTest <- peaktest_KO$FDR
#Add Correlation and FDR values from p2g linkage
linkedpeaks_matched_IFT_KO$Correlation_p2gLink <- p2g_diff_df_KO$Correlation
linkedpeaks_matched_IFT_KO$FDR_p2gLink <- p2g_diff_df_KO$FDR


#Manually construct loops object for only loops that are linked to 
#differential peaks

resolution <- 1
peakSummits <- resize(linkedpeaks_matched_IFT_KO, 1, "center") 
geneStarts <- resize(metadata(p2g_df)$geneSet[linkedpeaks_matched_IFT_KO$idxRNA], 1, "start") 

if(!is.null(resolution)){
  summitTiles <- floor(start(peakSummits) / resolution) * resolution + floor(resolution / 2)
  geneTiles <- floor(start(geneStarts) / resolution) * resolution + floor(resolution / 2)
}else{
  summitTiles <- start(peakSummits)
  geneTiles <- start(geneStarts)
}

.constructGR <- function(
    seqnames = NULL, 
    start = NULL, 
    end = NULL, 
    ignoreStrand = TRUE
){
  df <- data.frame(seqnames, start, end)
  idx <- which(df[,2] > df[,3])
  df[idx,2:3] <-  df[idx,3:2]
  if(!ignoreStrand){
    strand <- rep("+",nrow(df))
    strand[idx] <- "-" 
  }else{
    strand <- rep("*",nrow(df))
  }
  gr <- GRanges(df[,1], IRanges(df[,2],df[,3]), strand = strand)
  return(gr)
}

loops <- .constructGR( 
  seqnames = seqnames(peakSummits), 
  start = summitTiles, 
  end = geneTiles
) 
mcols(loops)$value <-linkedpeaks_matched_IFT_KO$Correlation_p2gLink 
mcols(loops)$FDR <- linkedpeaks_matched_IFT_KO$FDR_p2gLink 

loops <- loops[order(mcols(loops)$value, decreasing=TRUE)] 
loops <- unique(loops) 
loops <- loops[width(loops) > 0] 
loops <- sort(sortSeqlevels(loops)) 

loops <- SimpleList(Peak2GeneLinks = loops) 

#Custom palette for figure panel
pal = c("CMs_OFT_x_WT" = "#D1495B", "CMs_OFT_x_KO" = "#95555B", "CMs_V_x_WT" = "#EDAE49", "CMs_V_x_KO" = "#A37D41", "CMs_IFT_x_WT" = "#20A39E", "CMs_IFT_x_KO" = "#306281")

#Plot browser track
tracks <- plotBrowserTrack(
  ArchRProj = proj_Mef2c_v13_E85_subset, 
  groupBy = "CellTypeByGenotype", 
  useGroups = c("CMs_OFT_x_WT", "CMs_OFT_x_KO", "CMs_V_x_WT", "CMs_V_x_KO", "CMs_IFT_x_WT", "CMs_IFT_x_KO"),
  pal = pal,
  geneSymbol = c( "Wnt2"),
  upstream = 45000,
  downstream = 25000,
  features =  getMarkers(markerTest_CMs_IFT, cutOff = "FDR <= 0.15 & Log2FC >= 0.5", returnGR = TRUE)["CMs_IFT_x_KO"],
  ylim = c(0,0.98),
  loops = loops
)
plotPDF(tracks, name = "Fig3F_Wnt2_locus", width = 8, height = 7, ArchRProj = proj_Mef2c_v13_E85_subset, addDOC = F)

