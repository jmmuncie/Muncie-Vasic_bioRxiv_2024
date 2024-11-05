#This Python script loads the base GRN and adata objects for
#E7.75, E8.5, and E9 and then uses CellOracle to infer GRNs
#for WT and KO CMs or heart field cells at each timepoint

import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns

from scipy.io import mmread
from scipy.sparse import csr_matrix


from celloracle import motif_analysis as ma
import celloracle as co
co.__version__


plt.rcParams['figure.figsize'] = [6, 4.5]
plt.rcParams["savefig.dpi"] = 300

for condition_timepoint in [['WT', '775'], ['KO', '775'], ['WT', '85'], ['KO', '85'], ['WT', '9'], ['KO', '9']]:
    wt_or_ko = condition_timepoint[0]
    timepoint = condition_timepoint[1]
    
    print(f'{wt_or_ko} AND e{timepoint}')
    adata = sc.read_h5ad(f'data/adata_objects/e{timepoint}_subset.h5ad')
    
    # RAW DATA IN CASE YOU NEED IT
    raw_mtx = mmread(f"./data/adata_objects/e{timepoint}_matrix.mtx")
    raw_cells = pd.read_csv(f"./data/adata_objects/e{timepoint}_raw_cells.csv", header=None)
    raw_genes = pd.read_csv(f"./data/adata_objects/e{timepoint}_raw_genes.csv", header=None)
    x = pd.DataFrame(raw_mtx.toarray()).astype('float64')
    x.index = raw_genes.values.T[0]
    
    adata.X = x.loc[adata.var_names, :].values.T

    
    sc.pp.normalize_per_cell(adata)
    adata.raw = adata
    adata.layers["raw_count"] = adata.raw.X.copy()

    #########
    #########
    #########
    # Below, I have to map the numeric that was output when I transformed the seurat object to anndata. The dictionary provides this mapping.
    if timepoint == '775':
        e775_names = ['SHF_WT', 'SHF_KO', 'CMs/FHF_WT', 'CMs/FHF_KO', 'JCF_WT', 'JCF_KO', 'CrM_WT', 
                          'CrM_KO', 'ExM_WT', 'ExM_KO', 'SoM_WT', 'SoM_KO', 'LPM_WT', 'LPM_KO', 'NMPs_WT', 
                          'NMPs_KO', 'PrxM_WT', 'PrxM_KO', 'KPs_WT', 'KPs_KO', 'HSCs_WT', 'HSCs_KO']
        mapping_dict = dict(zip(range(0, len(e775_names)), e775_names))
        adata.obs['celltype_x_genotype'] = adata.obs['cell_type_pool_x_genotype'].map(mapping_dict)

    if timepoint == '85':
        e85_names = ['pSHF_WT','pSHF_KO', 'aSHF_WT', 'aSHF_KO', 'IFT-CMs_WT', 'IFT-CMs_KO', 'V-CMs_WT', 
                         'V-CMs_KO', 'OFT-CMs_WT', 'OFT-CMs_KO', 'PhM_WT', 'PhM_KO', 'LPM_WT', 'LPM_KO', 
                         'PostM_WT', 'PostM_KO', 'MixM_WT', 'MixM_KO', 'C16_WT', 'C16_KO']
        mapping_dict = dict(zip(range(0, len(e85_names)), e85_names))
        adata.obs['celltype_x_genotype'] = adata.obs['cell_type_pool_x_genotype'].map(mapping_dict)

    if timepoint == '9':
        e9_names = ['SHF_WT', 'SHF_KO', 'Pe_WT', 'Pe_KO', 'VP_WT', 'VP_KO', 'CMs-A_WT', 'CMs-A_KO', 
                        'CMs-AVC_WT', 'CMs-AVC_KO', 'CMs-V_WT', 'CMs-V_KO', 'CMs-OFT_WT', 'CMs-OFT_KO', 
                        'PhM_WT', 'PhM_KO', 'C11_WT', 'C11_KO']
        mapping_dict = dict(zip(range(0, len(e9_names)), e9_names))
        adata.obs['celltype_x_genotype'] = adata.obs['cell_type_pool_x_genotype'].map(mapping_dict)

    
    #########
    #########
    #########
    # Here, I subset the cardiomyocyte clusters that I will use for constructing the GRN.
    if timepoint == "775":
        adata_sub = adata[adata.obs['celltype_x_genotype'].isin([
                f'SHF_{wt_or_ko}',
                f'CMs/FHF_{wt_or_ko}'
            ]
            )]
        mapping_dict = {f'SHF_{wt_or_ko}': f'HF_{wt_or_ko}',
                        f'CMs/FHF_{wt_or_ko}': f'HF_{wt_or_ko}'}
        adata_sub.obs['celltype_x_genotype'] = adata.obs['celltype_x_genotype'].map(mapping_dict)


    if timepoint == "85":
        adata_sub = adata[adata.obs['celltype_x_genotype'].isin([
                f'IFT-CMs_{wt_or_ko}',
                f'V-CMs_{wt_or_ko}',
                f'OFT-CMs_{wt_or_ko}',
            ]
            )]

    if timepoint == "9":
        adata_sub = adata[adata.obs['celltype_x_genotype'].isin([
                f'CMs-AVC_{wt_or_ko}',
                f'CMs-A_{wt_or_ko}',
                f'CMs-OFT_{wt_or_ko}',
                f'CMs-V_{wt_or_ko}'
            ]
            )]

    base_GRN = pd.read_parquet(f'data/base_grn_outputs/E{timepoint}/{wt_or_ko}_base_GRN_dataframe.parquet')
    
    ### If you want to prune based on differential accessibility
    # diffexp = pd.read_csv(f'./data/base_grn_diff_access/e{timepoint}{wt_or_ko.lower()}.txt', sep = '\t')
    # removed_indexes = []
    # for index, row in diffexp.iterrows():
    #     if row['Log2FC'] < -1:
    #         removed_indexes.append(f"{row['seqnames']}_{row['start']}_{row['end']}")                
    # base_GRN = base_GRN[~base_GRN.peak_id.isin(removed_indexes)]
    ###
    
    oracle = co.Oracle()
    
    oracle.import_anndata_as_raw_count(adata=adata_sub,
                                   cluster_column_name="celltype_x_genotype",
                                   embedding_name="X_umap")
    
    oracle.import_TF_data(TF_info_matrix=base_GRN)
    
    # Perform PCA
    oracle.perform_PCA()

    # Select important PCs
    n_comps = np.where(np.diff(np.diff(np.cumsum(oracle.pca.explained_variance_ratio_))>0.002))[0][0]
    n_comps = min(n_comps, 50)
    
    n_cell = oracle.adata.shape[0]
    k = int(0.025*n_cell)
    
    oracle.knn_imputation(n_pca_dims=n_comps, k=k, balanced=True, b_sight=k*8,
                      b_maxl=k*4, n_jobs=4)
    
    # Save oracle object.
    oracle.to_hdf5(f"./data/celloracle/e{timepoint}/{wt_or_ko}_cardiac-subset.celloracle.oracle")
    
    # Load file.
    oracle = co.load_hdf5(f"./data/celloracle/e{timepoint}/{wt_or_ko}_cardiac-subset.celloracle.oracle")
    
    # Check clustering data
    sc.pl.umap(oracle.adata, color="celltype_x_genotype")
    oracle.adata
    
    print("Starting to infer network")
    # Calculate GRN for each population in "louvain_annot" clustering unit.
    # This step may take some time.(~30 minutes)
    links = oracle.get_links(cluster_name_for_GRN_unit="celltype_x_genotype", alpha=10,
                             verbose_level=10)
    
    links.filter_links(p=0.001, weight="coef_abs", threshold_number=4000)
    
    #./data/celloracle/ network scores.
    links.get_network_score()
    links.merged_score.head()
    
    links.to_hdf5(file_path=f"./data/celloracle/e{timepoint}/{wt_or_ko}_cardiac-subset-links.celloracle.links")