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

###############################################################################
#   Main analysis
###############################################################################

#   Read in AnnData
adata = sc.read_h5ad(input_adata_path)

#  Subset to just this sample
sample_name = adata.obs.sample_id.unique()[os.environ['SGE_TASK_ID'] - 1]
adata = adata[adata.obs.sample_id == sample_name, :]
print('Continuing with sample ' + sample_name + '...')

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
    adata = adata[adata.obs.key.isin(clusters.key), :]
    clusters.key.index = [x.split('_')[0] for x in clusters.key]
    clusters.cluster.index = [x.split('_')[0] for x in clusters.key]

    #   Now cluster calls and the AnnData should line up by key, so we add a
    #   cluster column
    assert(all(clusters.key == adata.obs.key))
    adata.obs['cluster'] = clusters.cluster

    #   Plot what the clusters look like, as a sanity check
    plot_color=[
        "#F56867","#FEB915","#C798EE","#59BE86","#7495D3","#D1D1D1","#6D1A9C",
        "#15821E","#3A84E6","#997273","#787878","#DB4C6C","#9E7A7A","#554236",
        "#AF5F3C","#93796C","#F9BD3F","#DAB370","#877F6C","#268785"
    ]
    adata.uns["cluster_colors"]=list(plot_color[:k])
    out_file = pyhere.here(this_out_dir_plots, "clusters.png")
    #plot_adata(adata, "cluster", "BayesSpace Clusters", plot_color, out_file)
    plot_adata(adata, "cluster", "BayesSpace Clusters", 'Accent', out_file)