import os,sys,csv,re
import pandas as pd
import numpy as np
import scanpy as sc
import math
import SpaGCN as spg
from scipy.sparse import issparse
import random, torch
import warnings
warnings.filterwarnings("ignore")
import matplotlib.colors as clr
import matplotlib.pyplot as plt
import SpaGCN as spg
import cv2
from scanpy import read_10x_h5
from anndata import AnnData
import pyhere

analysis_name = '05-bayesspace_meta'

out_dir_processed = pyhere.here('spagcn', 'processed-data', analysis_name)
out_dir_plots = pyhere.here('spagcn', 'plots', analysis_name)
input_adata_path = pyhere.here(
    'spagcn', 'processed-data', '04-spatial_dlpfc', 'spe_anndata.h5ad'
)
cluster_dir = '/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/clustering_results'
N_CLUSTERS = list(range(2, 16))
NUM_META_COLUMNS = 5

#  Column names expected to be present in adata.obs for X and Y array and pixel
#  coordinates, respectively
array_names = ("array_row", "array_col")
pixel_names = ("pxl_col_in_fullres", "pxl_row_in_fullres")

###############################################################################
#  Helper functions
###############################################################################

#  Given an AnnData, a column name in 'adata.obs', a color list, plot title, 
#  and output file path, create a PNG file containing a plot
def plot_adata(adata, colname, title, plot_color, out_file):
    ax = sc.pl.scatter(
        adata, alpha=1, x=pixel_names[0], y=pixel_names[1], color=colname, title=title,
        color_map=plot_color, show=False, size=100000/adata.shape[0]
    )
    ax.set_aspect('equal', 'box')
    ax.axes.invert_yaxis()
    plt.savefig(out_file, dpi=600)
    plt.close()

#  Given a raw Anndata, processed AnnData, X and Y array coordinates, a
#  target spatial domain, and a string specifying column name for cluster calls,
#  return a DataFrame of info about SVGs. Under special conditions, None may be
#  returned instead of a DataFrame, if no valid SVGs can be found.
def get_svgs(raw, adata, x_array, y_array, target, cluster_col):
    #Set filtering criterials
    min_in_group_fraction=0.8
    min_in_out_group_ratio=1
    min_fold_change=1.5
    
    #Search radius such that each spot in the target domain has approximately 10 neighbors on average
    adj_2d=spg.calculate_adj_matrix(x=x_array, y=y_array, histology=False)
    start, end= np.quantile(adj_2d[adj_2d!=0],q=0.001), np.quantile(adj_2d[adj_2d!=0],q=0.1)
    r=spg.search_radius(
        target_cluster=target, cell_id=adata.obs.index.tolist(), x=x_array,
        y=y_array, pred=adata.obs[cluster_col].tolist(), start=start, end=end,
        num_min=10, num_max=14,  max_run=100
    )
    
    #Detect neighboring domains
    nbr_domians=spg.find_neighbor_clusters(
        target_cluster=target,
        cell_id=raw.obs.index.tolist(), 
        x=raw.obs[array_names[0]].tolist(), 
        y=raw.obs[array_names[1]].tolist(), 
        pred=raw.obs[cluster_col].tolist(),
        radius=r,
        ratio=1/2
    )
    
    #  I've observed that domains where no neighbors are found are spatially
    #  scattered and are likely not biologically meaningful or interpretable.
    #  Stop here in that case, and don't find SVGs                               
    if nbr_domians is None:
        print('Found no neighbors for domain', target, '. No SVGs will be found for this domain.')
        return None
    
    nbr_domians=nbr_domians[0:3]
    de_genes_info=spg.rank_genes_groups(
        input_adata=raw,
        target_cluster=target,
        nbr_list=nbr_domians, 
        label_col=cluster_col, 
        adj_nbr=True, 
        log=True
    )
    
    #  Take significant genes only                                
    de_genes_info=de_genes_info[(de_genes_info["pvals_adj"]<0.05)]
    if de_genes_info.shape[0] == 0:
        print("No significant SVGs found for domain " + str(target) + "!")
        return None
    
    filtered_info=de_genes_info
    
    #  Take genes that exceed the expression seen in neighboring domains
    if np.count_nonzero(filtered_info["in_out_group_ratio"] > min_in_out_group_ratio) == 0:
        print("Can't find any genes such that more spots in domain " + str(target) + " express the gene than its neighbors!")
        return None
        
    filtered_info = filtered_info[filtered_info["in_out_group_ratio"] > min_in_out_group_ratio]
    
    #  Relax filtering criteria if required to find at least one SVG
    while np.count_nonzero(
        (filtered_info["in_group_fraction"] > min_in_group_fraction) &
        (filtered_info["fold_change"] > min_fold_change)
        ) == 0:
        #  Relax both parameters simultaneously
        min_in_group_fraction -= 0.1
        min_fold_change = (min_fold_change + 1) / 2
        
        if min_in_group_fraction <= 0:
            print('No "min_in_group_fraction" exists that finds any SVGs for domain ' + str(target) + "!")
            return None
        
        if min_fold_change <= 1.01:
            print('No "fold_change" > 1.01 exists that finds any SVGs for domain ' + str(target) + "!")
            return None
        
        print('Warning: lowering "min_in_group_fraction" to ', min_in_group_fraction, 'and "min_fold_change" to ', min_fold_change, 'for domain', target, 'to find at least 1 SVG.')
    
    #  Filter significant genes
    filtered_info=filtered_info[(filtered_info["in_group_fraction"] > min_in_group_fraction) &
                                (filtered_info["fold_change"] > min_fold_change)]
    
    filtered_info=filtered_info.sort_values(by="in_group_fraction", ascending=False)
    filtered_info["target_dmain"]=target
    filtered_info["neighbors"]=str(nbr_domians)
    
    return filtered_info

###############################################################################
#   Main analysis
###############################################################################

#   Read in AnnData
adata_orig = sc.read_h5ad(input_adata_path)

#  Subset to just this sample
sample_name = adata_orig.obs.sample_id.unique()[os.environ['SGE_TASK_ID'] - 1]
adata_orig = adata_orig[adata_orig.obs.sample_id == sample_name, :]
print('Continuing with sample ' + sample_name + '...')

adata_orig.var_names_make_unique()
spg.prefilter_genes(adata_orig,min_cells=3) # avoiding all genes are zeros
spg.prefilter_specialgenes(adata_orig)

#Normalize and take log for UMI
sc.pp.normalize_per_cell(adata_orig)
sc.pp.log1p(adata_orig)

#   Make the output directories specific to this sample and ensure they exist
out_dir_processed = pyhere.here(out_dir_processed, sample_name)
out_dir_plots = pyhere.here(out_dir_plots, sample_name)
if not os.path.exists(out_dir_processed):
    os.mkdir(out_dir_processed)

if not os.path.exists(out_dir_plots):
    os.mkdir(out_dir_plots)

for k in N_CLUSTERS:
    #   Use output directories specific to this number of clusters
    this_out_dir_plots = pyhere.here(out_dir_plots, str(k) + "_clusters")
    this_out_dir_processed = pyhere.here(
        out_dir_processed, str(k) + "_clusters"
    )
    if not os.path.exists(this_out_dir_plots):
        os.mkdir(this_out_dir_plots)

    if not os.path.exists(this_out_dir_processed):
        os.mkdir(this_out_dir_processed)

    #   Read in the 'clusters.csv' file
    cluster_path = cluster_dir + '/bayesSpace_harmony_' + str(k) + \
        '/clusters.csv'
    assert os.path.exists(cluster_path), 'clusters.csv file not found.'

    clusters = pd.read_csv(cluster_path)

    #   Subset to only the clusters for this sample
    clusters['sample_name'] = ['_'.join(x.split('_')[1:]) for x in clusters.key]
    clusters = clusters[clusters['sample_name'] == sample_name]

    #   Subset AnnData to spots where a cluster label exists
    adata = adata_orig[adata_orig.obs.key.isin(clusters.key), :]
    clusters.key.index = [x.split('_')[0] for x in clusters.key]
    clusters.cluster.index = [x.split('_')[0] for x in clusters.key]

    #   Now cluster calls and the AnnData should line up by key, so we add a
    #   cluster column
    assert(all(clusters.key == adata.obs.key))
    adata.obs['cluster'] = clusters.cluster.astype('category')

    #   Plot what the clusters look like, as a sanity check
    out_file = pyhere.here(this_out_dir_plots, "clusters.png")
    plot_adata(
        adata, "cluster", "BayesSpace Clusters", plot_color=None,
        out_file=out_file
    )

    #   Read AnnData in and subset to current sample. 'adata' is filtered,
    #   whereas 'raw' has all observations seen in the clusters dataframe
    raw = sc.read_h5ad(input_adata_path)
    raw = raw[raw.obs.sample_id == sample_name, :]
    raw.var_names_make_unique()

    #   Subset AnnData to spots where a cluster label exists
    raw = raw[raw.obs.key.isin(clusters.key), :]

    raw.obs["cluster"] = adata.obs["cluster"]

    #Convert sparse matrix to non-sparse
    raw.X=(raw.X.A if issparse(raw.X) else raw.X)
    raw.raw=raw
    sc.pp.log1p(raw)

    #   Now cluster calls should match between AnnData and data frame
    assert all(clusters.cluster.values == raw.obs.cluster.values)

    x_array=adata.obs[array_names[0]].astype(np.int32).tolist()
    y_array=adata.obs[array_names[1]].astype(np.int32).tolist()

        #  Initialize empty meta-gene columns
    for i in range(NUM_META_COLUMNS):
        clusters['meta_gene_' + str(i + 1)] = ''

    #  Find meta genes for each domain found
    for target in clusters.cluster.unique():
        #  Determine SVGs
        filtered_info = get_svgs(raw, adata, x_array, y_array, target, "cluster")
        if filtered_info is None:
            continue
        
        print("SVGs for domain ", str(target),":", filtered_info["genes"].tolist())
        
        #  Order SVGs by "in_group_fraction" and within that, "pvals_adj"
        filtered_info=filtered_info.sort_values(by="pvals_adj", ascending=True)
        filtered_info=filtered_info.sort_values(by="in_group_fraction", ascending=False)
        start_gene = filtered_info.genes.values[0]
        
        meta_name, meta_exp=spg.find_meta_gene(input_adata=raw,
                            pred=raw.obs["pred"].tolist(),
                            target_domain=target,
                            start_gene=start_gene,
                            mean_diff=0,
                            early_stop=True,
                            max_iter=3,
                            use_raw=False)
        
        raw.obs["meta"]=meta_exp
        
        #  Fill in meta-gene columns in 'cluster_list' that are associated with
        #  domain 'target'
        fill_meta_column(raw, meta_name, cluster_list, target, NUM_META_COLUMNS)
        
        color_self = clr.LinearSegmentedColormap.from_list('pink_green', ['#3AB370',"#EAE7CC","#FD1593"], N=256)
        
        #  Plot "start gene" in meta gene set
        title = get_symbol(raw, start_gene) + ' (marker for domain ' + str(target) + ')'
        out_file = pyhere.here(this_out_dir_plots, "start_meta_gene_domain_" + str(target) + ".png")
        raw.obs["exp"]=raw.X[:,raw.var.index==start_gene]
        plot_adata(raw, "exp", title, color_self, out_file)
        
        #  Plot entire meta gene set
        raw.obs["exp"]=raw.obs["meta"]
        out_file = pyhere.here(this_out_dir_plots, "all_meta_genes_domain_" + str(target) + ".png")
        plot_adata(raw, "exp", parse_metaname(raw, meta_name), color_self, out_file)


out_file = pyhere.here(this_out_dir_processed, 'clusters.csv')
cluster_list.to_csv(out_file)
