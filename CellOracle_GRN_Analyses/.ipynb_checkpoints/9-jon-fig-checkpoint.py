# 0. Import
import os
import sys

import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn2
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns

from scipy.io import mmread
from scipy.sparse import csr_matrix



bedtools_gcc_path = '/project/clark-saucerman/bedtools2/bin:/apps/software/standard/core/gcc/12.2.0/bin'
original_path = os.environ.get('PATH')
os.environ['PATH'] = bedtools_gcc_path + ':' + original_path
lib_path = '/apps/software/standard/core/gcc/12.2.0/lib64'  # Adjust the path if necessary
os.environ['LD_LIBRARY_PATH'] = lib_path + ':' + os.environ.get('LD_LIBRARY_PATH', '')


import warnings
warnings.filterwarnings('ignore')
from celloracle import motif_analysis as ma
import celloracle as co
co.__version__


# visualization settings
# %config InlineBackend.figure_format = 'retina'
# %matplotlib inline

plt.rcParams['figure.figsize'] = [6, 4.5]
plt.rcParams["savefig.dpi"] = 300


links_wt_e85 = co.load_hdf5(file_path="./data/celloracle/e85/WT_cardiac-subset-links.celloracle.links")
links_wt_e9 = co.load_hdf5(file_path="./data/celloracle/e9/WT_cardiac-subset-links.celloracle.links")

links_ko_e85 = co.load_hdf5(file_path="./data/celloracle/e85/KO_cardiac-subset-links.celloracle.links")
links_ko_e9 = co.load_hdf5(file_path="./data/celloracle/e9/KO_cardiac-subset-links.celloracle.links")


links_dict = {}

for timepoint, wt_or_ko in [['85', 'WT'], ['85', 'KO'], ['9', 'WT'], ['9', 'KO']]:
    links = co.load_hdf5(file_path=f"./data/celloracle/e{timepoint}/{wt_or_ko}_cardiac-subset-links.celloracle.links")

    links.filter_links(p=0.001, weight="coef_abs", threshold_number=4500)
    filtered_df = links.filtered_links

    for cluster in filtered_df.keys():
        df = filtered_df[cluster]
        
        if wt_or_ko == 'KO':
            df = df.drop(df[df['source'] == 'Mef2c'].index)

        df = df.sort_values('coef_abs', ascending=False).iloc[:4000,:]

        sources = df['source'].unique()
        targets = df['target'].unique()

        # Create a zero matrix
        matrix = pd.DataFrame(np.zeros((len(sources), len(targets))), index=sources, columns=targets)

        # Populate the matrix with coef_mean values
        for index, row in df.iterrows():
            matrix.at[row['source'], row['target']] = row['coef_mean']

        matrix.to_csv(f"./data/celloracle/e{timepoint}/{wt_or_ko}_{cluster}_gephi.csv")
        links.filtered_links[cluster] = df

    links_dict[f'{timepoint}_{wt_or_ko}'] = links
    
    
    
merged_scores = {}
for key, links in links_dict.items():
    links.get_network_score()
    merged_scores[key] = links.merged_score[links.merged_score.cluster.str.contains('IFT-|-A', regex=True)]
    
    
delta_centr = []
for idx in merged_scores['85_WT'].index:
    if idx not in merged_scores['85_KO'].index:
        continue
    else:
        arr = [idx] + list(merged_scores['85_KO'].loc[idx].values[:-1] - merged_scores['85_WT'].loc[idx].values[:-1])

    delta_centr.append(arr)
    
    
x = pd.DataFrame(delta_centr, columns=['Gene', 'All', 'Centr_All', 'In', 'Centr_In', 'Out', 'Centr_Out', 'Bet', 'Eig'])
x.sort_values('Out', ascending=False).head(55)
# x.sort_values('Out').head(40)
# x[x.Gene == 'Esrrg']


def remove_and_filter(links_df, genes_to_remove, num_edges):
    for g in genes_to_remove:
        links_df = links_df.drop(links_df[links_df.source=='Mef2c'].index)
        links_df = links_df.drop(links_df[links_df.target=='Mef2c'].index)
        
    return links_df.sort_values('coef_abs', ascending=False).iloc[:num_edges]


thresh_num = 4000

links_wt_e85 = co.load_hdf5(file_path="./data/celloracle/e85/WT_cardiac-subset-links.celloracle.links")
links_ko_e85 = co.load_hdf5(file_path="./data/celloracle/e85/KO_cardiac-subset-links.celloracle.links")

links_wt_e85.filter_links(p=0.001, weight="coef_abs", threshold_number=thresh_num)
links_ko_e85.filter_links(p=0.001, weight="coef_abs", threshold_number=thresh_num+1000)

links_ko_e85.filtered_links['IFT-CMs_KO'] = remove_and_filter(links_ko_e85.filtered_links['IFT-CMs_KO'], ['Mef2c'], thresh_num)


def return_adata_raw(timepoint):
    adata = sc.read_h5ad(f'data/adata_objects/{timepoint}_subset.h5ad')
    
    if timepoint == 'e85':
        names = ['pSHF_WT','pSHF_KO', 'aSHF_WT', 'aSHF_KO', 'IFT-CMs_WT', 'IFT-CMs_KO', 'V-CMs_WT', 
                 'V-CMs_KO', 'OFT-CMs_WT', 'OFT-CMs_KO', 'PhM_WT', 'PhM_KO', 'LPM_WT', 'LPM_KO', 
                 'PostM_WT', 'PostM_KO', 'MixM_WT', 'MixM_KO', 'C16_WT', 'C16_KO']
    elif timepoint == 'e9':
        names = ['SHF_WT', 'SHF_KO', 'Pe_WT', 'Pe_KO', 'VP_WT', 'VP_KO', 'CMs-A_WT', 'CMs-A_KO', 
                        'CMs-AVC_WT', 'CMs-AVC_KO', 'CMs-V_WT', 'CMs-V_KO', 'CMs-OFT_WT', 'CMs-OFT_KO', 
                        'PhM_WT', 'PhM_KO', 'C11_WT', 'C11_KO']        
    else:
        return
    
    mapping_dict = dict(zip(range(0, len(names)), names))
    adata.obs['celltype_x_genotype'] = adata.obs['cell_type_pool_x_genotype'].map(mapping_dict)    
    
    if timepoint == 'e9':
        adata.obs.loc[adata.obs['celltype_x_genotype'] == 'CMs-AVC_WT', 'celltype_x_genotype'] = 'CMs-A_WT'
    
    raw_mtx = mmread(f"./data/adata_objects/{timepoint}_matrix.mtx")
    raw_cells = pd.read_csv(f"./data/adata_objects/{timepoint}_raw_cells.csv", header=None)
    raw_genes = pd.read_csv(f"./data/adata_objects/{timepoint}_raw_genes.csv", header=None)
    x = pd.DataFrame(raw_mtx.toarray())
    x.index = raw_genes.values.T[0]
    
    x = x.T
    x.index = raw_cells.values.flatten()
    raw_cells.index = raw_cells.values.flatten()
    raw_genes.index = raw_genes.values.flatten()
    
    new_adata_raw = sc.AnnData(
        X=x.values,  # Use the normalized and log-transformed data
        var=raw_genes,  # Use the same genes
        obs=raw_cells   # Use the same cells
        )
    
    sc.pp.normalize_total(new_adata_raw, target_sum=1e4)
    sc.pp.log1p(new_adata_raw)
    adata.raw = new_adata_raw
    return adata



timepoint = 'e9'
adata = return_adata_raw(timepoint)


# Chamber is Atrial because assuming adata_E9 is your AnnData object
chamber = 'A'

# Filter the data to include only the relevant groups
if timepoint == 'e85':
    adata_subset = adata[adata.obs['celltype_x_genotype'].isin([f'{chamber}-CMs_KO', f'{chamber}-CMs_WT'])]
elif timepoint == 'e9':
    adata_subset = adata[adata.obs['celltype_x_genotype'].isin([f'CMs-{chamber}_KO', f'CMs-{chamber}_WT'])]
else:
    print('oops')
    
# adata_subset.X = adata_subset.X + 8

sc.tl.rank_genes_groups(adata_subset, groupby='celltype_x_genotype', method='wilcoxon')
# sc.tl.rank_genes_groups(adata_subset, groupby='celltype_x_genotype', method='wilcoxon', use_raw=False)

# Extract the results into a DataFrame
# import pandas as pd
if timepoint == 'e85':
    de_results = pd.DataFrame(
        {
            'genes': adata_subset.uns['rank_genes_groups']['names'][f'{chamber}-CMs_KO'],
            'logfoldchanges': adata_subset.uns['rank_genes_groups']['logfoldchanges'][f'{chamber}-CMs_KO'],
            'pvals': adata_subset.uns['rank_genes_groups']['pvals'][f'{chamber}-CMs_KO'],
            'pvals_adj': adata_subset.uns['rank_genes_groups']['pvals_adj'][f'{chamber}-CMs_KO']
        }
    )
elif timepoint == 'e9':
    de_results = pd.DataFrame(
        {
            'genes': adata_subset.uns['rank_genes_groups']['names'][f'CMs-{chamber}_KO'],
            'logfoldchanges': adata_subset.uns['rank_genes_groups']['logfoldchanges'][f'CMs-{chamber}_KO'],
            'pvals': adata_subset.uns['rank_genes_groups']['pvals'][f'CMs-{chamber}_KO'],
            'pvals_adj': adata_subset.uns['rank_genes_groups']['pvals_adj'][f'CMs-{chamber}_KO']
        }
    )
else:
    print('oops2')
    
de_results['genes'] = de_results['genes'].apply(lambda x: x[0] if isinstance(x, tuple) else x)


source_list = list(links_wt_e85.filtered_links['IFT-CMs_WT'].source.unique())
target_list = list(links_wt_e85.filtered_links['IFT-CMs_WT'].target.unique())
all_genes = source_list + target_list

x = links_wt_e85.filtered_links['IFT-CMs_WT']
x = x[x['source'] == 'Mef2c']
x = x[x['coef_mean']>0]
y = [print(z) for z in x.target.values]


sc.pl.rank_genes_groups_heatmap(
    adata_subset,
    key='rank_genes_groups', 
    figsize=(10, 8))

def get_coef_matrix(links_df):
    all_genes = np.unique(list(links_df['source'].values) + list(links_df['target'].values))

    all_coefficients = np.zeros((len(all_genes), len(all_genes)))
    coefficients_matrix = pd.DataFrame(all_coefficients, index=all_genes, columns=all_genes)

    for row in links_df.itertuples():
        coefficients_matrix.at[row.source, row.target] = row.coef_mean
    return coefficients_matrix


def get_gene_change(genes_to_change, coef_matrix, iterations):
    delta_x = np.zeros(coef_matrix.shape[0])

    for g in genes_to_change:
        delta_x[coef_matrix.index.get_loc(g)] = -1

    change = np.zeros(coef_matrix.shape[0])

    for i in range(0, iterations):
        change = change + np.matmul(delta_x, coef_matrix.values)
        delta_x = delta_x + change


    change_dict = dict(zip([x for x in coef_matrix.index], [x for x in change]))

    change_df = pd.DataFrame({'name': [x for x in coef_matrix.index],
                             'delta_x': [x for x in change]})


    pos_genes = []
    neg_genes = []

    for g, val in change_dict.items():
        if val > 0:
            pos_genes.append(g)
        if val < 0: 
            neg_genes.append(g)

    return change_dict, pos_genes, neg_genes, change_df



tf = 'Mef2c'

coef_matrix = get_coef_matrix(links_wt_e85.filtered_links['IFT-CMs_WT'])
gene_change, up_genes, down_genes, change_df = get_gene_change([tf], coef_matrix, iterations=1)
de_lost_in_ko = de_results[de_results.logfoldchanges < -1.5]
de_gained_in_ko = de_results[de_results.logfoldchanges >1.5]

i = 0
print('(POSITIVE RESULT) Genes our model predicts will decrease in KO and are actually downregulated in KO.')

for g in change_df[change_df['delta_x'] < -0.01].sort_values('delta_x')['name'].values:
    if g in de_lost_in_ko['genes'].values:
        print(g)
        i += 1
tp = i
print(i)
print('\n')

print('(POSITIVE RESULT) Genes our model predicts will increase in KO and actually increase in KO')
i = 0
for g in change_df[change_df['delta_x'] > 0.01].sort_values('delta_x')['name'].values:
    if g in de_gained_in_ko['genes'].values:
        print(g)
        i += 1
tn = i
print(i)
print('\n')

print('(NEGATIVE RESULT) Genes our model predicts will decrease in KO, BUT actually INCREASE in KO.')
i = 0
for g in change_df[change_df['delta_x'] < -0.01].sort_values('delta_x')['name'].values:
    if g in de_gained_in_ko['genes'].values:
        print(g)
        i += 1
fp = i
print(i)

print('\n')
print('(NEGATIVE RESULT) Genes our model predicts will increase in KO, BUT actually DECREASE in KO.')
i = 0
for g in change_df[change_df['delta_x'] > 0.01].sort_values('delta_x')['name'].values:
    if g in de_lost_in_ko['genes'].values:
        print(g)
        i += 1
fn = i
print(i)


print(f"{tf}: Accuracy of: {(tp + tn) / (tp+tn+fp+fn)}, and num tf + tn = {tp+tn}")



# First, import packages
import pandas as pd
import networkx as nx
import pyvis
import numpy as np
import celloracle as co
import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn2
import graphviz

# visualization settings
# %config InlineBackend.figure_format = 'retina'
# %matplotlib inline

plt.rcParams['figure.figsize'] = [6, 4.5]
plt.rcParams["savefig.dpi"] = 300



def remove_and_filter(links_df, genes_to_remove, num_edges):
    for g in genes_to_remove:
        links_df = links_df.drop(links_df[links_df.source=='Mef2c'].index)
        links_df = links_df.drop(links_df[links_df.target=='Mef2c'].index)
        
    return links_df.sort_values('coef_abs', ascending=False).iloc[:num_edges]


thresh_num=4000

links_wt_e85 = co.load_hdf5(file_path="./data/celloracle/e85/WT_cardiac-subset-links.celloracle.links")
links_ko_e85 = co.load_hdf5(file_path="./data/celloracle/e85/KO_cardiac-subset-links.celloracle.links")

links_wt_e85.filter_links(p=0.001, weight="coef_abs", threshold_number=thresh_num)
links_ko_e85.filter_links(p=0.001, weight="coef_abs", threshold_number=thresh_num+1000)

links_ko_e85.filtered_links['IFT-CMs_KO'] = remove_and_filter(links_ko_e85.filtered_links['IFT-CMs_KO'], ['Mef2c'], thresh_num)



def get_centralities(links, chamber_name):
    links.get_network_score()
    return links.merged_score[links.merged_score.cluster == chamber_name]


def plot_centrality_bar(genes, chamber, links_wt, links_ko, ylim):
    
    wt_cent = get_centralities(links_wt, f'{chamber}-CMs_WT')
    ko_cent = get_centralities(links_ko, f'{chamber}-CMs_KO')
    
    gene_centrality = {}
    # gene_centrality_ko = {}

    for g in genes:
        gene_centrality[g] = []
        if g in wt_cent.index:
            gene_centrality[g].append(wt_cent.loc[g,:].degree_out)
        else:
            gene_centrality[g].append(0)
        
        if g in ko_cent.index:
            gene_centrality[g].append(ko_cent.loc[g,:].degree_out) 
        else:
            gene_centrality[g].append(0)
      
    # import pdb
    # pdb.set_trace()
    labels = ['WT', 'KO']
    df = pd.DataFrame(gene_centrality, index=labels)
    ax = df.T.plot.bar()
    ax.set_ylim(0, ylim)
    ax.set_ylabel('Degree Out')
    
    
    
genes = ['Nr2f2', 'Gata4']
chamber = 'IFT'
plot_centrality_bar(genes, chamber, links_wt_e85, links_ko_e85, ylim=60)


def plot_subnetwork(network, tf_genes, targets):
    graph_df = {'source': [],
                 'target': [],
                 'weight': [],
                 'edge_cols': []}
    
    for curr_gene in tf_genes:
        curr_gene_df = network[network['source'] == curr_gene]
        gene_df_targs = curr_gene_df[curr_gene_df['target'].isin(targets)]

        for index, row in gene_df_targs.iterrows():
            graph_df['source'].append(row['source'])
            graph_df['target'].append(row['target'])
            graph_df['weight'].append(row['coef_mean'])
            if row['coef_mean'] > 0:
                graph_df['edge_cols'].append('green')
            else:
                graph_df['edge_cols'].append('red')

    graph_df = pd.DataFrame(graph_df)

    G = nx.from_pandas_edgelist(graph_df, source='source', target='target', edge_attr='weight', create_using=nx.DiGraph())
    from pyvis.network import Network
    net = Network(notebook=True, directed=True)
    net.from_nx(G)

    for edge in net.edges:
        if edge['width'] > 0:
            edge['color'] = '#FFB000'
            edge['label'] = '+'
        else:
            edge['color'] = '#D81B60'
            edge['label'] = '_'

        edge['value'] = np.abs(edge['width'])
        edge['arrowStrikethrough'] = False

    for node in net.nodes:
        if node['id'] in tf_genes:
            # node['color'] = '#1f77b4' # blue is atrial
            node['shape'] = 'diamond'
            node['label'] = node['label'] + ' (TF)'
        # elif node['id'] in top_eVent:
        #     # node['color'] = '#ff7f0e' # Orange is ventricular
        #     node['shape'] = 'dot' # 
        #     # node['label'] = node['label'] + ' (OTHER)'
        else:
            node['color'] = 'black' # Black is neither

    return net



x = links_wt_e85.filtered_links['IFT-CMs_WT']




x = links_wt_e85.merged_score
y = x[x['cluster'] == 'IFT-CMs_WT'].sort_values('betweenness_centrality', ascending=False).head(15).index


de_results.sort_values('logfoldchanges', ascending=False).head(20)

de_temp = de_results[de_results.logfoldchanges < 0].sort_values('pvals_adj')
de_temp.head(20).genes.values
# y = [print(z) for z in de_temp.head(20).genes.values]
# de_temp.head(20)
de_results.sort_values('logfoldchanges', ascending=False)
de_subset = de_results[de_results.pvals <1E-7]
# de_subset[de_subset.genes.isin(all_genes)].genes.values


tf_genes = ['Nr2f2', 'Gata4']

targets = links_wt_e85.filtered_links['IFT-CMs_WT'].target.unique().tolist()
targets += tf_genes

net = plot_subnetwork(links_wt_e85.filtered_links['IFT-CMs_WT'], tf_genes, targets)
# net = plot_subnetwork(links_ko_e85.filtered_links['IFT-CMs_KO'], tf_genes, targets)

y = [print(z) for z in tf_genes]

net.show('net.html')



def get_gene_targets(tf, links_wt, links_ko):
    wt_targets = links_wt[links_wt.sort_values('coef_mean')['source'].isin([tf])].target.unique()
    ko_targets = links_ko[links_ko.sort_values('coef_mean')['source'].isin([tf])].target.unique()
    
    intersect_targets = [x for x in wt_targets if x in ko_targets]
    
    wt_only_targets = [x for x in wt_targets if x not in ko_targets]

    ko_only_targets = [x for x in ko_targets if x not in wt_targets]
    
    return intersect_targets, wt_only_targets, ko_only_targets

tf = 'Nr2f2'

links_wt = links_wt_e85.filtered_links['IFT-CMs_WT']
links_ko = links_ko_e85.filtered_links['IFT-CMs_KO']

intersect_targets, wt_targets, ko_targets = get_gene_targets(tf, links_wt, links_ko)


import networkx as nx
import numpy as np
import math
import matplotlib

matplotlib.rc('font', family='Arial')
fig, ax = plt.subplots(figsize=(12, 12))

all_nodes = wt_targets + intersect_targets +  ko_targets
# all_nodes = list(set(wt_targets + intersect_targets))
# all_nodes = list(set(intersect_targets))

G_act = nx.DiGraph()
G_act.add_nodes_from(all_nodes)
G_inhib = nx.DiGraph()
G_inhib.add_nodes_from(all_nodes)

pos = nx.circular_layout(G_act.subgraph(all_nodes))

G_act.add_node(tf)
G_inhib.add_node(tf)

tf_links = links_wt[links_wt.source == tf]
tf_links = links_ko[links_ko.source == tf]
    
for x in tf_links.target:
    connection = tf_links[tf_links.target == x]
    if connection.coef_mean.values[0] < 0:
        G_inhib.add_edge(tf, x, weight=connection.coef_abs*40, color='#f542a4')
    else:
        G_act.add_edge(tf, x, weight=connection.coef_abs*40, color='black')
        

for i, G in enumerate([G_act, G_inhib]):
    edges = G.edges()
    colors = [G[u][v]['color'] for u,v in edges]
    weights = [G[u][v]['weight'] for u,v in edges]
    center_node = tf
    pos[center_node] = np.array([0, 0])
    
    pos_higher = nx.rescale_layout_dict(pos, 16)
    nx.draw_networkx_labels(G, pos_higher)
    
    pos_lower = nx.rescale_layout_dict(pos, 13)
    # nx.draw_networkx_nodes(G, pos, node_color="#005D32")

    if i == 0:
        nx.draw_networkx_edges(G, pos_lower, width=weights, edge_color=colors, arrows=True, arrowstyle=matplotlib.patches.ArrowStyle('->', head_length=0.4, head_width=0.5))
    else:
        nx.draw_networkx_edges(G, pos_lower, width=weights, edge_color=colors, arrows=True, arrowstyle=matplotlib.patches.ArrowStyle('|-|', widthA=0, angleA=0, widthB=.5, angleB=0))


pos_mid = nx.rescale_layout_dict(pos, 14)
nx.draw_networkx_nodes(G, pos_mid, node_color="#42f593")

x_values, y_values = zip(*pos_higher.values())
x_max = max(x_values)
x_min = min(x_values)
x_margin = (x_max - x_min) * 0.1
plt.xlim(x_min - x_margin, x_max + x_margin)

y_max = max(y_values)
y_min = min(y_values)
y_margin = (y_max - y_min) * 0.1
plt.ylim(y_min - y_margin, y_max + y_margin)
# plt.show()
plt.savefig("./nr2f2-ko.pdf", format='pdf')