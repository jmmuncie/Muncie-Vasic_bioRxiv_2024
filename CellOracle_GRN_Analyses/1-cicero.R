#This R script loads in an ArchR project and subsets a
#peak-by-cell matrix for a given genotype (e.g., knockout).
#Using Cicero and Monocle3, calculates peak co-accessibility,
#and exports results, including peak information and Cicero connections,
#for downstream gene regulatory network analysis.
#Paths to input data and output folders are specified, and relevant
#parameters (e.g., FDR cutoff, cell type) are set to tailor the
#analysis for genotype-specific insights.

library(cicero)
library(monocle3)
library(ArchR)

load_folder <- "./data/Mef2c_v13_E85_subset/"
proj_mef2c <- loadArchRProject(path = load_folder, force = FALSE, showLogo = TRUE)
save_folder <- "./data/base_grn_outputs/E85/"

# First, get matrix formatted correctly
peak_mat = getMatrixFromProject(proj_mef2c, useMatrix = "PeakMatrix", binarize = TRUE)
peak_assays <- assays(peak_mat)

peakinfo <- data.frame(chr = peak_mat@rowRanges@seqnames, 
                       bp1 = peak_mat@rowRanges@ranges@start, 
                       bp2 = peak_mat@rowRanges@ranges@start+500, 
                       site_name = paste(peak_mat@rowRanges@seqnames,
		       peak_mat@rowRanges@ranges@start,
		       peak_mat@rowRanges@ranges@start+500,  sep="_"))

row.names(peakinfo) <- peakinfo$site_name

cellinfo <- data.frame(cells = peak_mat@colData@rownames)
row.names(cellinfo) <- peak_mat@colData@rownames

matrix_in <- peak_assays$PeakMatrix
row.names(matrix_in) <- row.names(peakinfo)

for (cell_type in c("WT", "KO")) {
  matrix_in <- peak_assays$PeakMatrix
  row.names(matrix_in) <- row.names(peakinfo)
  peak_mat = getMatrixFromProject(proj_mef2c, useMatrix = "PeakMatrix", binarize = TRUE)
  peak_assays <- assays(peak_mat)
  
  mask = proj_mef2c$Genotype == cell_type
  
  temp_matrix = matrix_in[,mask]
  temp_cellinfo = data.frame(cells = cellinfo[mask,])
  row.names(temp_cellinfo) = cellinfo[mask,]
  
  input_cds <-  suppressWarnings(new_cell_data_set(temp_matrix,
                                                   cell_metadata = temp_cellinfo,
                                                   gene_metadata = peakinfo))
  
  input_cds <- detect_genes(input_cds)
  input_cds <- estimate_size_factors(input_cds)
  input_cds <- preprocess_cds(input_cds, method = "LSI")
  input_cds <- reduce_dimension(input_cds, reduction_method = 'UMAP', 
                                preprocess_method = "LSI")
  
  umap_coords <- reducedDims(input_cds)$UMAP
  
  rm(matrix_in, peak_assays, peak_mat)
  
  cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = umap_coords)
  
  chromosome_length <- read.table("./data/mm10_chromosome_length.txt")
  
  conns <- run_cicero(cicero_cds, chromosome_length) 
  
  head(conns)
  
  all_peaks <- row.names(exprs(input_cds))
  write.csv(x = all_peaks, file = paste0(save_folder, cell_type, "_peaks.csv"))
  write.csv(x = conns, file = paste0(save_folder, cell_type, "_cicero_connections.csv"))
}
