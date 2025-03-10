library(cicero)
library(monocle3)
library(ArchR)


## Uncomment if running Cicero on E8.5
# E85: CMs_V_x_KO, CMs_IFT_x_KO, CMs_IFT_x_WT, CMs_OFT_x_WT, CMs_V_x_WT, CMs_OFT_x_KO
load_folder <- "./data/Mef2c_v13_E85_subset/"
save_folder <- "./data/base_grn_outputs/E85/"
proj_mef2c <- loadArchRProject(path = load_folder, force = FALSE, showLogo = TRUE)
celltypes_to_keep <- c("CMs_V_x_KO", "CMs_IFT_x_KO", "CMs_IFT_x_WT", 
                       "CMs_OFT_x_WT", "CMs_V_x_WT", "CMs_OFT_x_KO")


## Uncomment if running Cicero on E9
# E9: CMs_V_x_KO, CMs_IFT_x_WT, CMs_IFT_x_KO, CMs_V_x_WT, CMs_OFT_x_KO, CMs_OFT_x_WT
# load_folder <- "./data/Mef2c_v13_E9_subset/"
# save_folder <- "./data/base_grn_outputs/E9/"
# proj_mef2c <- loadArchRProject(path = load_folder, force = FALSE, showLogo = TRUE)
# celltypes_to_keep <- c("CMs_V_x_KO", "CMs_IFT_x_WT", "CMs_IFT_x_KO",
#                        "CMs_V_x_WT", "CMs_OFT_x_KO", "CMs_OFT_x_WT")


# Extract cell metadata
cell_metadata <- getCellColData(proj_mef2c)
# Get the cell names that match your target cell types
cells_to_keep <- rownames(cell_metadata)[cell_metadata$CellTypeByGenotype %in% celltypes_to_keep]

# First, get matrix formatted correctly
peak_mat = getMatrixFromProject(proj_mef2c, useMatrix = "PeakMatrix", binarize = TRUE)
peak_mat <- peak_mat[, cells_to_keep]

peak_assays <- assays(peak_mat)

peakinfo <- data.frame(chr = peak_mat@rowRanges@seqnames, 
                       bp1 = peak_mat@rowRanges@ranges@start, 
                       bp2 = peak_mat@rowRanges@ranges@start+500, 
                       site_name = paste(peak_mat@rowRanges@seqnames, peak_mat@rowRanges@ranges@start, peak_mat@rowRanges@ranges@start+500,  sep="_"))
row.names(peakinfo) <- peakinfo$site_name

cellinfo <- data.frame(cells = peak_mat@colData@rownames)
row.names(cellinfo) <- peak_mat@colData@rownames

matrix_in <- peak_assays$PeakMatrix
row.names(matrix_in) <- row.names(peakinfo)

input_cds <-  suppressWarnings(new_cell_data_set(matrix_in,
                                                 cell_metadata = cellinfo,
                                                 gene_metadata = peakinfo))

input_cds <- detect_genes(input_cds)
input_cds <- estimate_size_factors(input_cds)
input_cds <- preprocess_cds(input_cds, method = "LSI")
input_cds <- reduce_dimension(input_cds, reduction_method = 'UMAP', 
                              preprocess_method = "LSI")

umap_coords <- reducedDims(input_cds)$UMAP

cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = umap_coords)

chromosome_length <- read.table("./data/mm10_chromosome_length.txt")

conns <- run_cicero(cicero_cds, chromosome_length) 

head(conns)

all_peaks <- row.names(exprs(input_cds))

write.csv(x = all_peaks, file = paste0(save_folder, "peaks.csv"))
write.csv(x = conns, file = paste0(save_folder, "cicero_connections.csv"))