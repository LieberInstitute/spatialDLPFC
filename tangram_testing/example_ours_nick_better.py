#  Rather than train tangram on top 200 genes based on top expression
#  (regardless of cell type), this script explores better ways to select
#  genes, incorporating diversity of expression across cell type

import os, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import torch
import tangram as tg
import seaborn as sns
import scipy as sp

out_dir = '/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spython/tangram_testing/example_ours_nick_better_out'
sc_path = '/dcl01/lieber/ajaffe/Nick/spatial/tangram/sce_dlpfc.h5ad'
sp_path = '/dcl01/lieber/ajaffe/Nick/spatial/tangram/visium_dlpfc.h5ad'

#  Given a numpy array "arr" and integer "N", return a list of the indices of
#  the top N largest values in arr
def topN(arr, N):
    arr_copy = arr.copy()
    min_val = np.min(arr)
    
    inds = []
    for i in range(N):
        temp = int(np.argmax(arr_copy))
        inds.append(temp)
        arr_copy[temp] = min_val - 1
    
    return inds

ad_sp = sc.read_h5ad(sp_path)
ad_sc = sc.read_h5ad(sc_path)

sc.pp.normalize_total(ad_sc)

#  Subset by those genes which overlap
sp_genes = np.unique(list(ad_sp.var.gene_name)).tolist()
sc_genes = np.unique(list(ad_sc.var.Symbol)).tolist()
overlap = [x for x in sp_genes if x in sc_genes]
ad_sc, ad_sp = tg.pp_adatas(ad_sc, ad_sp, genes=overlap)

prop_names = ['propNucleiExprs', 'propExprsIn.Astro', 'propExprsIn.Excit.ambig', 'propExprsIn.Excit.L2:3', 'propExprsIn.Excit.L3:4', 'propExprsIn.Excit.L4:5', 'propExprsIn.Excit.L5', 'propExprsIn.Excit.L5:6', 'propExprsIn.Excit.L6.broad', 'propExprsIn.Inhib.1', 'propExprsIn.Inhib.2', 'propExprsIn.Inhib.3', 'propExprsIn.Inhib.4', 'propExprsIn.Inhib.5', 'propExprsIn.Inhib.6', 'propExprsIn.Micro', 'propExprsIn.Oligo', 'propExprsIn.OPC']

#  Get proportion of cell types, for each gene, that have nonzero expression for that gene
props = np.array([[ad_sc.var[x][y] for x in prop_names]
                  for y in range(len(ad_sc.var.Symbol))])
props = np.mean(props > 0, axis=1)

#  Normalize into a metric with mean 0 and std. dev. 1 (across genes)
props = (props - np.mean(props)) / np.sqrt(np.var(props))

#  Get average number of total counts across all voxels
counts = np.asarray(np.mean(ad_sp.X, axis=0)).reshape(-1)

#  Normalize into a metric with mean 0 and std. dev. 1 (across genes)
counts = (counts - np.mean(counts)) / np.sqrt(np.var(counts))

#  Choose genes based on equally weighting the 2 previously defined metrics
#  equally
inds = topN(counts + props, 1000)
sp_genes = list(ad_sp.var.gene_name)
our_markers = [sp_genes[x] for x in inds]

#  Subset based on these genes
ad_sc, ad_sp = tg.pp_adatas(ad_sc, ad_sp, genes=our_markers)

assert ad_sc.var.index.equals(ad_sp.var.index)

#  Mapping step using GPU
ad_map = tg.map_cells_to_space(
    adata_cells=ad_sc,
    adata_space=ad_sp,
    device='cuda: 0' # device='cpu'
)

ad_map.write_h5ad(os.path.join(out_dir, 'ad_map.h5ad'))

#  Reload original data
ad_sp = sc.read_h5ad(sp_path)
ad_sc = sc.read_h5ad(sc_path)

ad_ge = tg.project_genes(adata_map=ad_map, adata_sc=ad_sc)
ad_ge.write_h5ad(os.path.join(out_dir, 'ad_ge.h5ad'))

df_all_genes = tg.compare_spatial_geneexp(ad_ge, ad_sp)

#  Compute average cosine similarity for test genes
np.mean(df_all_genes.score[np.logical_not(df_all_genes.is_training)])
# 0.16692365277114757 (200 training genes)
# 0.15575753851758423 (1000 training genes)

#  Compute average cosine similarity for training genes
np.mean(df_all_genes.score[df_all_genes.is_training])
# 0.906375931435494 (200 training genes)
# 0.7673710390238581 (1000 training genes)
