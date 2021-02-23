import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import torch
import tangram as tg
import seaborn as sns
import scipy as sp
from os import path

out_dir = 'out/'
sc_path = 'out/sce_dlpfc.h5ad'
sp_path = 'out/visium_dlpfc.h5ad'

os.chdir('/home/arta/Documents/GitHub/spython/local_tangram_expdata')


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

xs = ad_sp.obs.row.values
ys = ad_sp.obs.col.values

#  Manually save voxel coords plot
f = plt.figure()
plt.axis('off')
plt.scatter(xs, ys, s=.7)
f.savefig(os.path.join(out_dir, 'voxel_coords.pdf'), bbox_inches='tight')

ad_sc = sc.read_h5ad(sc_path)

sc.pp.normalize_total(ad_sc)

ad_sc.obs.cell_type.value_counts()

#  Compare which genes are present among:
#    visium data "sp_genes"
#    single cell data "sc_genes"

sp_genes = np.unique(list(ad_sp.var.gene_name)).tolist()
sc_genes = np.unique(list(ad_sc.var.Symbol)).tolist()
overlap = [x for x in sp_genes if x in sc_genes]
len(sp_genes)  # 36591
len(sc_genes)  # 33514
len(overlap)  # 31908

ad_sc, ad_sp = tg.pp_adatas(ad_sc, ad_sp, genes=overlap)

# prop_names = ['propNucleiExprs', 'propExprsIn.Astro', 'propExprsIn.Excit.ambig', 'propExprsIn.Excit.L2:3', 'propExprsIn.Excit.L3:4', 'propExprsIn.Excit.L4:5', 'propExprsIn.Excit.L5', 'propExprsIn.Excit.L5:6', 'propExprsIn.Excit.L6.broad', 'propExprsIn.Inhib.1', 'propExprsIn.Inhib.2', 'propExprsIn.Inhib.3', 'propExprsIn.Inhib.4', 'propExprsIn.Inhib.5', 'propExprsIn.Inhib.6', 'propExprsIn.Micro', 'propExprsIn.Oligo', 'propExprsIn.OPC']

# props = np.array([[ad_sc.var[x][y] for x in prop_names]
#                  for y in range(len(ad_sc.var.Symbol))])
# ent = sp.stats.entropy(props, axis=1)
# ent[np.isnan(ent)] = 0

#  Form a (flat) array of counts by gene
temp = np.sum(ad_sp.X, axis=0)
temp = np.asarray(temp).reshape(-1)

inds = topN(temp, 200)
sp_genes = list(ad_sp.var.gene_name)
our_markers = [sp_genes[x] for x in inds]

ad_sc, ad_sp = tg.pp_adatas(ad_sc, ad_sp, genes=our_markers)

assert ad_sc.var.index.equals(ad_sp.var.index)

if path.exists("out/ad_map.h5ad"):
    ad_map = sc.read("out/ad_map.h5ad")
else:
    #  Mapping step using GPU
    ad_map = tg.map_cells_to_space(
        adata_cells=ad_sc,
        adata_space=ad_sp,
        device='cuda: 0'  # device='cpu'
    )
    ad_map.write_h5ad(os.path.join(out_dir, 'ad_map.h5ad'))
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import torch
import tangram as tg
import seaborn as sns
import scipy as sp
from os import path

out_dir = 'out/'
sc_path = 'out/sce_dlpfc.h5ad'
sp_path = 'out/visium_dlpfc.h5ad'

# os.chdir('/home/arta/Documents/GitHub/spython/local_tangram_expdata')

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

xs = ad_sp.obs.row.values
ys = ad_sp.obs.col.values

#  Manually save voxel coords plot
f = plt.figure()
plt.axis('off')
plt.scatter(xs, ys, s=.7)
f.savefig(os.path.join(out_dir, 'voxel_coords.pdf'), bbox_inches='tight')

ad_sc = sc.read_h5ad(sc_path)

sc.pp.normalize_total(ad_sc)

ad_sc.obs.cell_type.value_counts()

#  Compare which genes are present among:
#    visium data "sp_genes"
#    single cell data "sc_genes"

sp_genes = np.unique(list(ad_sp.var.gene_name)).tolist()
sc_genes = np.unique(list(ad_sc.var.Symbol)).tolist()
overlap = [x for x in sp_genes if x in sc_genes]
len(sp_genes) # 36591
len(sc_genes) # 33514
len(overlap) # 31908

ad_sc, ad_sp = tg.pp_adatas(ad_sc, ad_sp, genes=overlap)

#prop_names = ['propNucleiExprs', 'propExprsIn.Astro', 'propExprsIn.Excit.ambig', 'propExprsIn.Excit.L2:3', 'propExprsIn.Excit.L3:4', 'propExprsIn.Excit.L4:5', 'propExprsIn.Excit.L5', 'propExprsIn.Excit.L5:6', 'propExprsIn.Excit.L6.broad', 'propExprsIn.Inhib.1', 'propExprsIn.Inhib.2', 'propExprsIn.Inhib.3', 'propExprsIn.Inhib.4', 'propExprsIn.Inhib.5', 'propExprsIn.Inhib.6', 'propExprsIn.Micro', 'propExprsIn.Oligo', 'propExprsIn.OPC']

#props = np.array([[ad_sc.var[x][y] for x in prop_names]
#                  for y in range(len(ad_sc.var.Symbol))])
#ent = sp.stats.entropy(props, axis=1)
#ent[np.isnan(ent)] = 0

#  Form a (flat) array of counts by gene
temp = np.sum(ad_sp.X, axis=0)
temp = np.asarray(temp).reshape(-1)

inds = topN(temp, 200)
sp_genes = list(ad_sp.var.gene_name)
our_markers = [sp_genes[x] for x in inds]

ad_sc, ad_sp = tg.pp_adatas(ad_sc, ad_sp, genes=our_markers)

assert ad_sc.var.index.equals(ad_sp.var.index)

if path.exists("out/ad_map.h5ad"):
    ad_map = sc.read("out/ad_map.h5ad")
else:
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

tg.plot_cell_annotation(ad_map, annotation='cell_type', x='row', y='col', nrows=5, ncols=4)
f = plt.gcf()
f.savefig(os.path.join(out_dir, 'cell_annotation.png'), bbox_inches='tight')

tg.plot_training_scores(ad_map, bins=50, alpha=.5)
f = plt.gcf()
f.savefig(os.path.join(out_dir, 'train_scores.pdf'), bbox_inches='tight')

ad_map.uns['train_genes_df']

ad_ge = tg.project_genes(adata_map=ad_map, adata_sc=ad_sc)
ad_ge.write_h5ad(os.path.join(out_dir, 'ad_ge.h5ad'))

#genes = ['Cdh12', 'Satb2', 'Dscaml1']
#ad_map.uns['train_genes_df'].loc[genes]
#tg.plot_genes(genes, adata_measured=ad_sp, adata_predicted=ad_ge)

#mask = ad_map.uns['train_genes_df'].train_score.isna()
#genes = ad_map.uns['train_genes_df'][mask].index.values
#tg.plot_genes(genes, adata_measured=ad_sp, adata_predicted=ad_ge)

(ad_ge.var.is_training == False).sum()

del ad_map, ad_sc, f, sp_genes, temp, xs, ys

df_all_genes = tg.compare_spatial_geneexp(ad_ge, ad_sp)

#  Compute average cosine similarity for test genes
np.mean(df_all_genes.score[np.logical_not(df_all_genes.is_training)])
# 0.16738437239477638

#  Compute average cosine similarity for training genes
np.mean(df_all_genes.score[df_all_genes.is_training])
# 0.9083821825226959

#  Reload original data
ad_sp = sc.read_h5ad(sp_path)
ad_sc = sc.read_h5ad(sc_path)

tg.plot_cell_annotation(ad_map, annotation='cell_type', x='row', y='col', nrows=5, ncols=4)
f = plt.gcf()
f.savefig(os.path.join(out_dir, 'cell_annotation.png'), bbox_inches='tight')

tg.plot_training_scores(ad_map, bins=50, alpha=.5)
f = plt.gcf()
f.savefig(os.path.join(out_dir, 'train_scores.pdf'), bbox_inches='tight')

ad_map.uns['train_genes_df']

ad_ge = tg.project_genes(adata_map=ad_map, adata_sc=ad_sc)
ad_ge.write_h5ad(os.path.join(out_dir, 'ad_ge.h5ad'))

# genes = ['Cdh12', 'Satb2', 'Dscaml1']
# ad_map.uns['train_genes_df'].loc[genes]
# tg.plot_genes(genes, adata_measured=ad_sp, adata_predicted=ad_ge)

# mask = ad_map.uns['train_genes_df'].train_score.isna()
# genes = ad_map.uns['train_genes_df'][mask].index.values
# tg.plot_genes(genes, adata_measured=ad_sp, adata_predicted=ad_ge)

(ad_ge.var.is_training == False).sum()

del ad_map, ad_sc, f, sp_genes, temp, xs, ys

df_all_genes = tg.compare_spatial_geneexp(ad_ge, ad_sp)

print("Compute average cosine similarity for test genes", file=sys.stdout)
print(np.mean(df_all_genes.score[np.logical_not(df_all_genes.is_training)]), file=sys.stdout)

print("Compute average cosine similarity for training genes", file=sys.stdout)
print(np.mean(df_all_genes.score[df_all_genes.is_training]), file=sys.stdout)
