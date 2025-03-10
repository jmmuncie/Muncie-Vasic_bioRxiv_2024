###----------------Description-------------------------------------

#This script is written to create P2G-DEG correlation plots 

###---------------Load Libraries------------------------------
library(ArchR); library(Seurat); library(dplyr)
library(GenomeInfoDb); library(ensembldb)
library(ggplot2); library(patchwork); library(ggrepel); library(tune)
set.seed(1234)
library(BSgenome) 
library(hexbin)

###-------------Load ArchR Project------------------------------------
proj_Mef2c_v13_E85_subset <- loadArchRProject("~/Mef2c_ArchR_working/Mef2c_v13_E85_subset")
setwd("~/Mef2c_ArchR_working")


###----------------Pairwise Marker Features CMs IFT KO v WT----------------------

markerTest_CMs_IFT <- getMarkerFeatures(
  ArchRProj = proj_Mef2c_v13_E85_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "CellTypeByGenotype",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "CMs_IFT_x_KO",
  bgdGroups = "CMs_IFT_x_WT"
)


###--------------------WT DARs--------------------------------------

markerTest_CMs_IFT_list_upinWT <- getMarkers(markerTest_CMs_IFT, cutOff = "FDR <= 0.15 & Log2FC <= -0.5", returnGR = TRUE)

#Let's convert the upinWT list to a GRanges object, sorted by Log2FC
markerpeaksgr_IFT_WT <- markerTest_CMs_IFT_list_upinWT$`CMs_IFT_x_KO`
order_IFT_WT <- order(markerpeaksgr_IFT_WT$Log2FC, start(markerpeaksgr_IFT_WT))
markerpeaks_sorted_IFT_WT <- markerpeaksgr_IFT_WT[order_IFT_WT]


###----------------WT DARs linked to GEX (P2G)---------------------

#Retrieve P2G links from dataset
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

#Now we want to find matches between our WT IFT DARs and the linked peaks:
#'findOverlaps()' will identify matching Granges 
match_WT_DARs_P2G <- findOverlaps(markerpeaks_sorted_IFT_WT, linked_peaks)
#Suset both lists for matches
markerpeaks_sorted_IFT_WT_matched <- markerpeaks_sorted_IFT_WT[match_WT_DARs_P2G@from]
linkedpeaks_matched_IFT_WT <- linked_peaks[match_WT_DARs_P2G@to]
#Note that these two granges above match, only that the markerpeaks object 
#contains the metadata inlcuding log2fc from the differential peak analysis

#Intersect these peaks with the full peak set to get nearestGene and peakType info
#Subset matching peaks
pgr <- getPeakSet(proj_Mef2c_v13_E85_subset)
sub_pgr_WT_IFT_DARs_linked_indexes <- findOverlaps(linkedpeaks_matched_IFT_WT, pgr)
sub_pgr_WT_IFT_DARs_linked <- pgr[sub_pgr_WT_IFT_DARs_linked_indexes@to]

#Finally, we need to go back and pull the linked differential peaks from the 
#linked peaks data frame, which contains metadata for the p2g links that we 
#will want to use
P2G_WT_DARs_df <- p2g_df[match_WT_DARs_P2G@to,]

#Now let's add all the relevant metadata to the linkedpeaks_matched_IFT_WT GRanges object
#First let's add the P2G gene names
linkedpeaks_matched_IFT_WT$GeneName <- metadata(p2g_df)$geneSet$name[linkedpeaks_matched_IFT_WT$idxRNA]
#Next, add Log2FC and FDR values from the diff peak testing
linkedpeaks_matched_IFT_WT$Log2FC_DiffPeakTest <- markerpeaks_sorted_IFT_WT_matched$Log2FC
linkedpeaks_matched_IFT_WT$FDR_DiffPeakTest <- markerpeaks_sorted_IFT_WT_matched$FDR
#Add Correlation and FDR values from p2g linkage
linkedpeaks_matched_IFT_WT$Correlation_p2gLink <- P2G_WT_DARs_df$Correlation
linkedpeaks_matched_IFT_WT$FDR_p2gLink <- P2G_WT_DARs_df$FDR
#Finally, add peakType and nearestGene/TSS info from sub_pgr
linkedpeaks_matched_IFT_WT$peakType <- sub_pgr_WT_IFT_DARs_linked$peakType
linkedpeaks_matched_IFT_WT$nearestGene <- sub_pgr_WT_IFT_DARs_linked$nearestGene
linkedpeaks_matched_IFT_WT$distToGeneStart <- sub_pgr_WT_IFT_DARs_linked$distToGeneStart
linkedpeaks_matched_IFT_WT$nearestTSS <- sub_pgr_WT_IFT_DARs_linked$nearestTSS
linkedpeaks_matched_IFT_WT$distToTSS <- sub_pgr_WT_IFT_DARs_linked$distToTSS

#Create data frame for all WT DARs with P2G link
WT_DARs_P2G_bed_df <- data.frame(
  chr = linkedpeaks_matched_IFT_WT@seqnames,
  start = linkedpeaks_matched_IFT_WT@ranges@start,
  end = linkedpeaks_matched_IFT_WT@ranges@start + 500,
  mm10PeakRange = paste0(linkedpeaks_matched_IFT_WT@seqnames,':',linkedpeaks_matched_IFT_WT@ranges@start,'-',linkedpeaks_matched_IFT_WT@ranges@start + 500),
  peakType = linkedpeaks_matched_IFT_WT$peakType,
  peakLog2FC = linkedpeaks_matched_IFT_WT$Log2FC_DiffPeakTest,
  peakFDR = linkedpeaks_matched_IFT_WT$FDR_DiffPeakTest,
  LinkedGene = linkedpeaks_matched_IFT_WT$GeneName, 
  P2G_Correlation = linkedpeaks_matched_IFT_WT$Correlation_p2gLink,
  P2G_FDR = linkedpeaks_matched_IFT_WT$FDR_p2gLink,
  nearestGene = linkedpeaks_matched_IFT_WT$nearestGene,
  distToNearestGeneStart = linkedpeaks_matched_IFT_WT$distToGeneStart,
  nearestTSS = linkedpeaks_matched_IFT_WT$nearestTSS,
  distToNearestTSS = linkedpeaks_matched_IFT_WT$distToTSS
)


###------------------KO DARs------------------------------------------

markerTest_CMs_IFT_list_upinKO <- getMarkers(markerTest_CMs_IFT, cutOff = "FDR <= 0.15 & Log2FC >= 0.5", returnGR = TRUE)

#Let's convert the upinKO list to a GRanges object, sorted by Log2FC
markerpeaksgr_IFT_KO <- markerTest_CMs_IFT_list_upinKO$`CMs_IFT_x_KO`
order_IFT_KO <- order(markerpeaksgr_IFT_KO$Log2FC, decreasing = T, start(markerpeaksgr_IFT_KO))
markerpeaks_sorted_IFT_KO <- markerpeaksgr_IFT_KO[order_IFT_KO]


###----------------KO DARs linked to GEX (P2G)---------------------

#Now we want to find matches between our KO IFT DARs and the linked peaks
#(which we already retrieved above):

#'findOverlaps()' will identify matching Granges 
match_KO_DARs_P2G <- findOverlaps(markerpeaks_sorted_IFT_KO, linked_peaks)
#Suset both lists for matches
markerpeaks_sorted_IFT_KO_matched <- markerpeaks_sorted_IFT_KO[match_KO_DARs_P2G@from]
linkedpeaks_matched_IFT_KO <- linked_peaks[match_KO_DARs_P2G@to]
#Note that these two granges above match, only that the markerpeaks object 
#contains the metadata inlcuding log2fc from the differential peak analysis

#Intersect these peaks with the full peak set to get nearestGene and peakType info
#Subset matching peaks
sub_pgr_KO_IFT_DARs_linked_indexes <- findOverlaps(linkedpeaks_matched_IFT_KO, pgr)
sub_pgr_KO_IFT_DARs_linked <- pgr[sub_pgr_KO_IFT_DARs_linked_indexes@to]

#Finally, we need to go back and pull the linked differential peaks from the 
#linked peaks data frame, which contains metadata for the p2g links that we 
#will want to use
P2G_KO_DARs_df <- p2g_df[match_KO_DARs_P2G@to,]

#Now let's add all the relevant metadata to the linkedpeaks_matched_IFT_KO GRanges object
#First let's add the P2G gene names
linkedpeaks_matched_IFT_KO$GeneName <- metadata(p2g_df)$geneSet$name[linkedpeaks_matched_IFT_KO$idxRNA]
#Next, add Log2FC and FDR values from the diff peak testing
linkedpeaks_matched_IFT_KO$Log2FC_DiffPeakTest <- markerpeaks_sorted_IFT_KO_matched$Log2FC
linkedpeaks_matched_IFT_KO$FDR_DiffPeakTest <- markerpeaks_sorted_IFT_KO_matched$FDR
#Add Correlation and FDR values from p2g linkage
linkedpeaks_matched_IFT_KO$Correlation_p2gLink <- P2G_KO_DARs_df$Correlation
linkedpeaks_matched_IFT_KO$FDR_p2gLink <- P2G_KO_DARs_df$FDR
#Finally, add peakType and nearestGene/TSS info from sub_pgr
linkedpeaks_matched_IFT_KO$peakType <- sub_pgr_KO_IFT_DARs_linked$peakType
linkedpeaks_matched_IFT_KO$nearestGene <- sub_pgr_KO_IFT_DARs_linked$nearestGene
linkedpeaks_matched_IFT_KO$distToGeneStart <- sub_pgr_KO_IFT_DARs_linked$distToGeneStart
linkedpeaks_matched_IFT_KO$nearestTSS <- sub_pgr_KO_IFT_DARs_linked$nearestTSS
linkedpeaks_matched_IFT_KO$distToTSS <- sub_pgr_KO_IFT_DARs_linked$distToTSS

#Create data frame for all KO DARs with P2G link
KO_DARs_P2G_bed_df <- data.frame(
  chr = linkedpeaks_matched_IFT_KO@seqnames,
  start = linkedpeaks_matched_IFT_KO@ranges@start,
  end = linkedpeaks_matched_IFT_KO@ranges@start + 500,
  mm10PeakRange = paste0(linkedpeaks_matched_IFT_KO@seqnames,':',linkedpeaks_matched_IFT_KO@ranges@start,'-',linkedpeaks_matched_IFT_KO@ranges@start + 500),
  peakType = linkedpeaks_matched_IFT_KO$peakType,
  peakLog2FC = linkedpeaks_matched_IFT_KO$Log2FC_DiffPeakTest,
  peakFDR = linkedpeaks_matched_IFT_KO$FDR_DiffPeakTest,
  LinkedGene = linkedpeaks_matched_IFT_KO$GeneName, 
  P2G_Correlation = linkedpeaks_matched_IFT_KO$Correlation_p2gLink,
  P2G_FDR = linkedpeaks_matched_IFT_KO$FDR_p2gLink,
  nearestGene = linkedpeaks_matched_IFT_KO$nearestGene,
  distToNearestGeneStart = linkedpeaks_matched_IFT_KO$distToGeneStart,
  nearestTSS = linkedpeaks_matched_IFT_KO$nearestTSS,
  distToNearestTSS = linkedpeaks_matched_IFT_KO$distToTSS
)


###----------------Combine Data Frames------------------------------

DARs_P2G <- rbind(WT_DARs_P2G_bed_df, KO_DARs_P2G_bed_df)


###---------------Get DEGs from Seurat------------------------------

#Load Seurat object
load("~/Desktop/Mef2c_v13_Seurat_scTrans_working/Seurat_Objects/mef2c_v13_E85_scTrans_v3_harmony_subset.Robj")

#Get DEGs
#Test for DEG between KO and WT cells in IFT CMs clusters (3 clusters)
Idents(mef2c_v13_E85_v3_harmony_subset) <- "harmony_cell_type_subset"
DEGs <- FindMarkers(mef2c_v13_E85_v3_harmony_subset, asay = "SCT", ident.1 = "KO", ident.2 = "WT", group.by = "Genotype", subset.ident = c("CMs-IFT1", "CMs-IFT2", "CMs-IFT3"))
DEGs <- subset(DEGs, p_val_adj < 0.05)
DEGs <- DEGs[order(DEGs$avg_log2FC, decreasing = TRUE), ]
names(DEGs)[3] <- "pct.KO"
names(DEGs)[4] <- "pct.WT"

#Identify rows in the DARs_P2G df where 8th column (linked gene) matches rownames of DEGs
matching_rows <- DARs_P2G[, 8] %in% rownames(DEGs)

#Extract matching rows from DARs_P2G and corresponding rows from DEGs
DARs_matching <- DARs_P2G[matching_rows, ]
DEGs_matching <- DEGs[match(DARs_P2G[matching_rows, 8], rownames(DEGs)), ]

#Bind together into single dataframe
DARs_P2G_with_DEG_info <- cbind(DARs_matching, DEGs_matching)
#Clear row names (currently row numbers from original DARs_P2G df)
rownames(DARs_P2G_with_DEG_info) <- NULL
#Update column names for DEG info
col_names <- colnames(DARs_P2G_with_DEG_info)
col_names[15:19] <- paste("DEG_", col_names[15:19], sep = "")
colnames(DARs_P2G_with_DEG_info) <- col_names

#Remove all X chromosome data (due to unbalanced sex of WT/KO embryos)
DARs_P2G_with_DEG_info <- DARs_P2G_with_DEG_info[DARs_P2G_with_DEG_info$chr != "chrX", ]

#Save dataframe to a CSV file - Table S4C
write.csv(DARs_P2G_with_DEG_info, "~/Desktop/Working_Directory/TableS4C_E85_IFT-CMs_DARs_P2G_with_DEG_info.csv", row.names = FALSE)


###------------------Create Plot-------------------------------------

#Get the indices for the peaks with the smallest and largest Log2FC Diff Peak Test
smallest_x_indices <- order(DARs_P2G_with_DEG_info[, 6])[1:5]
largest_x_indices <- order(DARs_P2G_with_DEG_info[, 6], decreasing = T)[1:5]

#Get the indices for the peaks with the smallest and largest Avg Log2FC DEG
#Note: for these values we have to deal with repeats (i.e. peaks can link to the 
#same gene, so the top 5 smallest DEG values could all belong to the same gene)
DEG_values <- DARs_P2G_with_DEG_info[, 16]
ordered_indices <- order(DEG_values)
sorted_values <- DEG_values[ordered_indices]
unique_values <- unique(sorted_values)
smallest_unique_values <- unique_values[1:5]
smallest_unique_y_indices <- match(smallest_unique_values, DARs_P2G_with_DEG_info[, 16])
unique_values_decreasing <- rev(unique_values)
largest_unique_values <- unique_values_decreasing[1:5]
largest_unique_y_indices <- match(largest_unique_values, DARs_P2G_with_DEG_info[, 16])

#Add manual labels
Actn2_indices <- which(DARs_P2G_with_DEG_info[, 8] == "Actn2")
first_Actn2 <- Actn2_indices[1]
Myh6_indices <- which(DARs_P2G_with_DEG_info[, 8] == "Myh6")
first_Myh6 <- Myh6_indices[1]
Myh7_indices <- which(DARs_P2G_with_DEG_info[, 8] == "Myh7")
first_Myh7 <- Myh7_indices[1]
Actc1_indices <- which(DARs_P2G_with_DEG_info[, 8] == "Actc1")
first_Actc1 <- Actc1_indices[1]
Gata4_indices <- which(DARs_P2G_with_DEG_info[, 8] == "Gata4")
first_Gata4 <- Gata4_indices[2]
Gata6_indices <- which(DARs_P2G_with_DEG_info[, 8] == "Gata6")
first_Gata6 <- Gata6_indices[2]
Sfrp5_indices <- which(DARs_P2G_with_DEG_info[, 8] == "Sfrp5")
first_Sfrp5 <- Sfrp5_indices[1]
Nkx25_indices <- which(DARs_P2G_with_DEG_info[, 8] == "Nkx2-5")
first_Nkx25 <- Nkx25_indices[1]
Tnnt2_indices <- which(DARs_P2G_with_DEG_info[, 8] == "Tnnt2")
first_Tnnt2 <- Tnnt2_indices[1]
Tnni1_indices <- which(DARs_P2G_with_DEG_info[, 8] == "Tnni1")
first_Tnni1 <- Tnni1_indices[1]
Wnt2_indices <- which(DARs_P2G_with_DEG_info[, 8] == "Wnt2")
first_Wnt2 <- Wnt2_indices[1]

#Compile indices for all points to label
man_indices <- c(smallest_x_indices, smallest_unique_y_indices, largest_x_indices, largest_unique_y_indices,
                 first_Actn2, first_Myh6, first_Myh7, first_Actc1, first_Gata4,
                 first_Gata6, first_Sfrp5, first_Nkx25, first_Tnnt2, first_Tnni1, 
                 first_Wnt2)
unique_indices <- unique(man_indices)

#Add column for label
DARs_P2G_with_DEG_info$label = ""

#Add gene name labels only to points we want to label
DARs_P2G_with_DEG_info$label[unique_indices] = DARs_P2G_with_DEG_info[unique_indices,8]

#Plot
myplot <- 
  ggplot(data = DARs_P2G_with_DEG_info, 
         aes(x = DARs_P2G_with_DEG_info[, 6],
             y = DARs_P2G_with_DEG_info[, 16],
             label = label
         )
  ) +
  geom_point(aes(color = P2G_Correlation), size = 2.5) +
  geom_text_repel(box.padding = 1, max.overlaps = Inf) +
  labs(x = "Log2FC Diff Peak Test", 
       y = "Avg Log2FC DEG", 
       title = "P2G/DEG Correlation Plot") +
  theme_classic() +
  scale_color_gradientn(colors = paletteContinuous(set = "comet")) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  geom_vline(xintercept = 0, linetype = "solid", color = "black") +
  coord_cartesian(xlim = range(c(min(DARs_P2G_with_DEG_info[, 6]), 
                                 max(DARs_P2G_with_DEG_info[, 6]))), 
                  ylim = range(c(min(DARs_P2G_with_DEG_info[, 16]), 
                                 max(DARs_P2G_with_DEG_info[, 16])))
  ) +
  coord_fixed() +
  theme(aspect.ratio = 0.65)

myplot
pdf("~/Desktop/Working_Directory/P2G_DEG_Correlation_Plots/Fig3D_E85_IFT-CMs_P2G_DEG_corr.pdf", width = 10, height = 8)
print(myplot)
dev.off()

