#This R scripts loads in an ArchR project, identifies the most differentially
#accessible regions between WT and KO, and then saves the peaks with
#FDR <=.01

library(cicero)
library(monocle3)
library(ArchR)


timepoint = "775"
load_folder <- paste0("./data/Mef2c_v13_E", timepoint, "_subset/")
proj_mef2c <- loadArchRProject(path = load_folder, force = FALSE, showLogo = TRUE)

markersPeaks <- getMarkerFeatures(
  ArchRProj = proj_mef2c,
  useMatrix = "PeakMatrix",
  groupBy = "Genotype",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01")

write.table(markerList$WT, file = paste0("./data/base_grn_diff_access/E", timepoint , "wt.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(markerList$KO, file = paste0("./data/base_grn_diff_access/E", timepoint , "ko.txt"), sep = "\t", row.names = FALSE, quote = FALSE)