#Load ArchR library and projects  
library(ArchR)
proj_Mef2c_v13_E85_subset <- loadArchRProject("~/Mef2c_ArchR_working/Mef2c_v13_E85_subset")
proj_Mef2c_v13_E85_E9_subset <- loadArchRProject("~/Mef2c_ArchR_working/Mef2c_v13_E85_E9_subset")

#Set working directory
setwd("~/Mef2c_ArchR_working")


###-----------------I) To Determine Number of Overlapping DARS-------------

#IMPORTANT NOTE: Because I had to create the second project with E8.5 and E9 combined
#to call DARs in OFT CMs, the peak sets are not identical. This means when comparing
#between OFT and either or both of the other two segments, I am returning *overlapping*
#peaks, not necessarily identical peaks. For this reason, for comparisons that involve
#OFT, I will write out bed files for the both the OFT peak ranges and the IFT/V peak ranges

#Get list of DARs you want to compare, convert to GRanges objects

#i): DARs accessible in WT/KO IFT CMs
markerTest_CMs_IFT <- getMarkerFeatures(
  ArchRProj = proj_Mef2c_v13_E85_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "CellTypeByGenotype",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "CMs_IFT_x_KO",
  bgdGroups = "CMs_IFT_x_WT"
)
markerTest_CMs_IFT_list_upinWT <- getMarkers(markerTest_CMs_IFT, cutOff = "FDR <= 0.15 & Log2FC <= -0.5", returnGR = TRUE)
markerTest_CMs_IFT_list_upinKO <- getMarkers(markerTest_CMs_IFT, cutOff = "FDR <= 0.15 & Log2FC >= 0.5", returnGR = TRUE)
#Let's convert the lists to GRanges objects, sorted by Log2FC
markerpeaksgr_IFT_WT <- markerTest_CMs_IFT_list_upinWT$`CMs_IFT_x_KO`
markerpeaksgr_IFT_KO <- markerTest_CMs_IFT_list_upinKO$`CMs_IFT_x_KO`
order_IFT_WT <- order(markerpeaksgr_IFT_WT$Log2FC, start(markerpeaksgr_IFT_WT))
order_IFT_KO <- order(markerpeaksgr_IFT_KO$Log2FC, decreasing = T, start(markerpeaksgr_IFT_KO))
markerpeaks_sorted_IFT_WT <- markerpeaksgr_IFT_WT[order_IFT_WT]
markerpeaks_sorted_IFT_KO <- markerpeaksgr_IFT_KO[order_IFT_KO]

#ii): DARs accessible in WT/KO V CMs
markerTest_CMs_V <- getMarkerFeatures(
  ArchRProj = proj_Mef2c_v13_E85_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "CellTypeByGenotype",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "CMs_V_x_KO",
  bgdGroups = "CMs_V_x_WT"
)
markerTest_CMs_V_list_upinWT <- getMarkers(markerTest_CMs_V, cutOff = "FDR <= 0.15 & Log2FC <= -0.5", returnGR = TRUE)
markerTest_CMs_V_list_upinKO <- getMarkers(markerTest_CMs_V, cutOff = "FDR <= 0.15 & Log2FC >= 0.5", returnGR = TRUE)
#Let's convert the lists to GRanges objects, sorted by Log2FC
markerpeaksgr_V_WT <- markerTest_CMs_V_list_upinWT$`CMs_V_x_KO`
markerpeaksgr_V_KO <- markerTest_CMs_V_list_upinKO$`CMs_V_x_KO`
order_V_WT <- order(markerpeaksgr_V_WT$Log2FC, start(markerpeaksgr_V_WT))
order_V_KO <- order(markerpeaksgr_V_KO$Log2FC, decreasing = T, start(markerpeaksgr_V_KO))
markerpeaks_sorted_V_WT <- markerpeaksgr_V_WT[order_V_WT]
markerpeaks_sorted_V_KO <- markerpeaksgr_V_KO[order_V_KO]

#iii): DARs accessible in WT OFT CMs
markerTest_CMs_OFT <- getMarkerFeatures(
  ArchRProj = proj_Mef2c_v13_E85_E9_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "CellTypeByGenotype",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "CMs_OFT_x_KO",
  bgdGroups = "CMs_OFT_x_WT"
)
markerTest_CMs_OFT_list_upinWT <- getMarkers(markerTest_CMs_OFT, cutOff = "FDR <= 0.15 & Log2FC <= -0.5", returnGR = TRUE)
markerTest_CMs_OFT_list_upinKO <- getMarkers(markerTest_CMs_OFT, cutOff = "FDR <= 0.15 & Log2FC >= 0.5", returnGR = TRUE)
#Let's convert the lists to GRanges objects, sorted by Log2FC
markerpeaksgr_OFT_WT <- markerTest_CMs_OFT_list_upinWT$`CMs_OFT_x_KO`
markerpeaksgr_OFT_KO <- markerTest_CMs_OFT_list_upinKO$`CMs_OFT_x_KO`
order_OFT_WT <- order(markerpeaksgr_OFT_WT$Log2FC, start(markerpeaksgr_OFT_WT))
order_OFT_KO <- order(markerpeaksgr_OFT_KO$Log2FC, decreasing = T, start(markerpeaksgr_OFT_KO))
markerpeaks_sorted_OFT_WT <- markerpeaksgr_OFT_WT[order_OFT_WT]
markerpeaks_sorted_OFT_KO <- markerpeaksgr_OFT_KO[order_OFT_KO]

#Now find overlaps between the Granges for each of the pairwise comparisons:
matching_DARs_WT_IFT_V_indices <- findOverlaps(markerpeaks_sorted_IFT_WT, markerpeaks_sorted_V_WT)
matching_DARs_KO_IFT_V_indices <- findOverlaps(markerpeaks_sorted_IFT_KO, markerpeaks_sorted_V_KO)
matching_DARs_WT_IFT_OFT_indices <- findOverlaps(markerpeaks_sorted_IFT_WT, markerpeaks_sorted_OFT_WT)
matching_DARs_KO_IFT_OFT_indices <- findOverlaps(markerpeaks_sorted_IFT_KO, markerpeaks_sorted_OFT_KO)
matching_DARs_WT_OFT_V_indices <- findOverlaps(markerpeaks_sorted_OFT_WT, markerpeaks_sorted_V_WT)
matching_DARs_KO_OFT_V_indices <- findOverlaps(markerpeaks_sorted_OFT_KO, markerpeaks_sorted_V_KO)

#Count number of overlaps for each comparison
length(matching_DARs_WT_IFT_V_indices) #343
length(matching_DARs_KO_IFT_V_indices) #23
length(matching_DARs_WT_IFT_OFT_indices) #376
length(matching_DARs_KO_IFT_OFT_indices) #5
length(matching_DARs_WT_OFT_V_indices) #366
length(matching_DARs_KO_OFT_V_indices) #5

#Subset the matching peaks from the V/IFT project
matching_DARs_WT_IFT_V <- markerpeaks_sorted_IFT_WT[matching_DARs_WT_IFT_V_indices@from]
matching_DARs_KO_IFT_V <- markerpeaks_sorted_IFT_KO[matching_DARs_KO_IFT_V_indices@from]
matching_DARs_WT_IFT_OFT <- markerpeaks_sorted_IFT_WT[matching_DARs_WT_IFT_OFT_indices@from]
matching_DARs_KO_IFT_OFT <- markerpeaks_sorted_IFT_KO[matching_DARs_KO_IFT_OFT_indices@from]
matching_DARs_WT_OFT_V <- markerpeaks_sorted_V_WT[matching_DARs_WT_OFT_V_indices@to]
matching_DARs_KO_OFT_V <- markerpeaks_sorted_V_KO[matching_DARs_KO_OFT_V_indices@to]
#And for OFT comparisons, subset the peak ranges that come from the OFT project
matching_DARs_WT_IFT_OFT_ranges_OFT <- markerpeaks_sorted_OFT_WT[matching_DARs_WT_IFT_OFT_indices@to]
matching_DARs_KO_IFT_OFT_ranges_OFT <- markerpeaks_sorted_OFT_KO[matching_DARs_KO_IFT_OFT_indices@to]
matching_DARs_WT_OFT_V_ranges_OFT <- markerpeaks_sorted_OFT_WT[matching_DARs_WT_OFT_V_indices@from]
matching_DARs_KO_OFT_V_ranges_OFT <- markerpeaks_sorted_OFT_KO[matching_DARs_KO_OFT_V_indices@from]

#Check that lengths are correct
length(matching_DARs_WT_IFT_V) #343
length(matching_DARs_KO_IFT_V) #23
length(matching_DARs_WT_IFT_OFT) #376
length(matching_DARs_WT_IFT_OFT_ranges_OFT) #376
length(matching_DARs_KO_IFT_OFT) #5
length(matching_DARs_KO_IFT_OFT_ranges_OFT) #5
length(matching_DARs_WT_OFT_V) #366
length(matching_DARs_WT_OFT_V_ranges_OFT) #366
length(matching_DARs_KO_OFT_V) #5
length(matching_DARs_KO_OFT_V_ranges_OFT) #5

#Now find overlaps for DARs in all three segments
matching_DARs_WT_IFT_V_OFT_indices <- findOverlaps(matching_DARs_WT_IFT_OFT, matching_DARs_WT_OFT_V)
matching_DARs_WT_IFT_V_OFT_ranges_OFT_indices <- findOverlaps(matching_DARs_WT_IFT_OFT_ranges_OFT, matching_DARs_WT_OFT_V_ranges_OFT)
matching_DARs_KO_IFT_V_OFT_indices <- findOverlaps(matching_DARs_KO_IFT_OFT, matching_DARs_KO_OFT_V)
matching_DARs_KO_IFT_V_OFT_ranges_OFT_indices <- findOverlaps(matching_DARs_KO_IFT_OFT_ranges_OFT, matching_DARs_KO_OFT_V_ranges_OFT)

#Count number of overlaps
#Note: slightly different result for WT depending on if using the ranges from V/IFT project
#or the OFT projec to do the 3-way comparison. That's okay, for now
length(matching_DARs_WT_IFT_V_OFT_indices) #205
length(matching_DARs_WT_IFT_V_OFT_ranges_OFT_indices) #203
length(matching_DARs_KO_IFT_V_OFT_indices) #2
length(matching_DARs_KO_IFT_V_OFT_ranges_OFT_indices) #2

#Subset the matching peaks from the both projects
matching_DARs_WT_IFT_V_OFT <- matching_DARs_WT_IFT_OFT[matching_DARs_WT_IFT_V_OFT_indices@from]
matching_DARs_WT_IFT_V_OFT_ranges_OFT <- matching_DARs_WT_IFT_OFT_ranges_OFT[matching_DARs_WT_IFT_V_OFT_ranges_OFT_indices@from]
matching_DARs_KO_IFT_V_OFT <- matching_DARs_KO_IFT_OFT[matching_DARs_KO_IFT_V_OFT_indices@from]
matching_DARs_KO_IFT_V_OFT_ranges_OFT <- matching_DARs_KO_IFT_OFT_ranges_OFT[matching_DARs_KO_IFT_V_OFT_ranges_OFT_indices@from]

#Check that lengths are corrrect
length(matching_DARs_WT_IFT_V_OFT) #205
length(matching_DARs_WT_IFT_V_OFT_ranges_OFT) #203
length(matching_DARs_KO_IFT_V_OFT) #2
length(matching_DARs_KO_IFT_V_OFT_ranges_OFT) #2

#Now use a for loop to create dataframes and write out beds:
#Set up to loop
DARs <- list(overlapping_DARs_WT_IFT_V = matching_DARs_WT_IFT_V, 
             overlapping_DARs_KO_IFT_V = matching_DARs_KO_IFT_V, 
             overlapping_DARs_WT_IFT_OFT = matching_DARs_WT_IFT_OFT, 
             overlapping_DARs_KO_IFT_OFT = matching_DARs_KO_IFT_OFT, 
             overlapping_DARs_WT_OFT_V = matching_DARs_WT_OFT_V, 
             overlapping_DARs_KO_OFT_V = matching_DARs_KO_OFT_V, 
             overlapping_DARs_WT_IFT_V_OFT = matching_DARs_WT_IFT_V_OFT, 
             overlapping_DARs_KO_IFT_V_OFT = matching_DARs_KO_IFT_V_OFT,
             overlapping_DARs_WT_IFT_OFT_ranges_OFT = matching_DARs_WT_IFT_OFT_ranges_OFT, 
             overlapping_DARs_KO_IFT_OFT_ranges_OFT = matching_DARs_KO_IFT_OFT_ranges_OFT, 
             overlapping_DARs_WT_OFT_V_ranges_OFT = matching_DARs_WT_OFT_V_ranges_OFT, 
             overlapping_DARs_KO_OFT_V_ranges_OFT = matching_DARs_KO_OFT_V_ranges_OFT, 
             overlapping_DARs_WT_IFT_V_OFT_ranges_OFT = matching_DARs_WT_IFT_V_OFT_ranges_OFT, 
             overlapping_DARs_KO_IFT_V_OFT_ranges_OFT = matching_DARs_KO_IFT_V_OFT_ranges_OFT
            )

lengths <- list(overlapping_DARs_WT_IFT_V = "343", 
                overlapping_DARs_KO_IFT_V = "23", 
                overlapping_DARs_WT_IFT_OFT = "376", 
                overlapping_DARs_KO_IFT_OFT = "5", 
                overlapping_DARs_WT_OFT_V = "366", 
                overlapping_DARs_KO_OFT_V = "5", 
                overlapping_DARs_WT_IFT_V_OFT = "205", 
                overlapping_DARs_KO_IFT_V_OFT = "2",
                overlapping_DARs_WT_IFT_OFT_ranges_OFT = "376", 
                overlapping_DARs_KO_IFT_OFT_ranges_OFT = "5", 
                overlapping_DARs_WT_OFT_V_ranges_OFT = "366", 
                overlapping_DARs_KO_OFT_V_ranges_OFT = "5", 
                overlapping_DARs_WT_IFT_V_OFT_ranges_OFT = "203", 
                overlapping_DARs_KO_IFT_V_OFT_ranges_OFT = "2"
               )

#Loop:
for (DAR_names in names(DARs)) {
  bed_df <- data.frame(
    chr = DARs[[DAR_names]]@seqnames,
    start = DARs[[DAR_names]]@ranges@start,
    end = DARs[[DAR_names]]@ranges@start + 500,
    PeakRange = paste0(DARs[[DAR_names]]@seqnames,':',DARs[[DAR_names]]@ranges@start,'-',DARs[[DAR_names]]@ranges@start + 500),
    peakLog2FC = DARs[[DAR_names]]$Log2FC,
    peakFDR = DARs[[DAR_names]]$FDR
  )
  #Write out bed
  write.table(bed_df, paste0('~/Mef2c_ArchR_working/Overlapping_DAR_beds/Mef2c_E85_',lengths[DAR_names],'_',DAR_names,'.bed'), sep="\t", row.names=FALSE, col.names=TRUE, quote = FALSE)
}


###-----------------II) To Generate Venn Diagrams-----------------------


#Note: Ended up using Eulerr shiny app online to create final plots because you can add the counts
#https://eulerr.co/

library(eulerr)

#Mef2c E8.5 WT DARs and overlap

fit_WT <- euler(c('IFT' = 2086, 'V' = 990, 'OFT' = 1028, 'IFT&V' = 343, 'IFT&OFT' = 376, 'OFT&V' = 366, 'IFT&V&OFT' = 205))

pdf(file = "~/Mef2c_ArchR_working/Overlapping_Peaks/VD_Mef2c_E85_Lost_DARs.pdf",   
    width = 3, 
    height = 3)

plot(fit_WT, fill=c('#20A39E', '#EDAE49', '#D1495B'))

dev.off()

#Mef2c E8.5 KO DARs and overlap

fit_KO <- euler(c('IFT' = 890, 'V' = 189, 'OFT' = 108, 'IFT&V' = 23, 'IFT&OFT' = 5, 'OFT&V' = 5, 'IFT&V&OFT' = 2))

pdf(file = "~/Mef2c_ArchR_working/Overlapping_Peaks/VD_Mef2c_E85_Gained_DARs.pdf",   
    width = 3, 
    height = 3)

plot(fit_KO, fill=c('#20A39E', '#EDAE49', '#D1495B'))

dev.off()

