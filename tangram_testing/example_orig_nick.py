#  This script is intended to be run from the working directory:
#    tangram_testing/Tangram/example

import os, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import torch
import tangram as tg
import seaborn as sns

out_dir = '/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spython/tangram_testing/example_orig_nick_out'

path = os.path.join('data', 'slideseq_MOp_1217.h5ad')
ad_sp = sc.read_h5ad(path)

xs = ad_sp.obs.x.values
ys = ad_sp.obs.y.values

#  Manually save voxel coords plot
f = plt.figure()
plt.axis('off')
plt.scatter(xs, ys, s=.7);
f.savefig(os.path.join(out_dir, 'voxel_coords.pdf'), bbox_inches='tight')

path = os.path.join('data','mop_sn_tutorial.h5ad')
ad_sc = sc.read_h5ad(path)

sc.pp.normalize_total(ad_sc)

ad_sc.obs.subclass_label.value_counts()

df_genes = pd.read_csv('data/MOp_markers.csv', index_col=0)
markers = np.reshape(df_genes.values, (-1, ))
markers = list(markers)

ad_sc, ad_sp = tg.pp_adatas(ad_sc, ad_sp, genes=markers)

assert ad_sc.var.index.equals(ad_sp.var.index)

ad_sc.write_h5ad(os.path.join(out_dir, 'ad_sc_readytomap.h5ad'))
ad_sp.write_h5ad(os.path.join(out_dir, 'ad_sp_readytomap.h5ad'))

#  Mapping step using GPU
ad_map = tg.map_cells_to_space(
    adata_cells=ad_sc,
    adata_space=ad_sp,
    device='cuda: 0' # device='cpu'
)

ad_map.write_h5ad(os.path.join(out_dir, 'ad_map.h5ad'))

#  Reload original data (not normalized)
path = os.path.join('data', 'slideseq_MOp_1217.h5ad')
ad_sp = sc.read_h5ad(path)

path = os.path.join('data','mop_sn_tutorial.h5ad')
ad_sc = sc.read_h5ad(path)

# ad_map = sc.read_h5ad(os.path.join(out_dir,'ad_map.h5ad'))

tg.plot_cell_annotation(ad_map, annotation='subclass_label', nrows=5, ncols=4)
tg.plot_training_scores(ad_map, bins=50, alpha=.5)

ad_map.uns['train_genes_df']

ad_ge = tg.project_genes(adata_map=ad_map, adata_sc=ad_sc)

genes = ['Cdh12', 'Satb2', 'Dscaml1']
ad_map.uns['train_genes_df'].loc[genes]

tg.plot_genes(genes, adata_measured=ad_sp, adata_predicted=ad_ge)

mask = ad_map.uns['train_genes_df'].train_score.isna()
genes = ad_map.uns['train_genes_df'][mask].index.values
tg.plot_genes(genes, adata_measured=ad_sp, adata_predicted=ad_ge)

(ad_ge.var.is_training == False).sum()

df_all_genes = tg.compare_spatial_geneexp(ad_ge, ad_sp)

#  Compute average cosine similarity for test genes
np.mean(df_all_genes.score[np.logical_not(df_all_genes.is_training)])
# 0.1616565621617533

#  Compute average cosine similarity for training genes
np.mean(df_all_genes.score[df_all_genes.is_training])
# 0.7694395580684324

def score_by_sparsity(df, sparsity):
    perc_genes = round(100 * np.mean(df.sparsity_2 < sparsity), 2)
    print(perc_genes, '% of genes are less than ', round(100 * sparsity, 2), '% sparse.', sep='')
    
    return np.mean(df.score[np.logical_and(np.logical_not(df.is_training),
                                           df.sparsity_2 < sparsity)])


sns_plot = sns.scatterplot(data=df_all_genes, x='score', y='sparsity_2', hue='is_training', alpha=.5)
fig = sns_plot.get_figure()
fig.savefig(os.path.join(out_dir, 'sparsity.pdf'))

genes = ['Snap25', 'Atp1b1', 'Atp1a3', 'Ctgf', 'Nefh', 'Aak1', 'Fa2h', ]
df_all_genes.loc[genes]

tg.plot_genes(genes, adata_measured=ad_sp, adata_predicted=ad_ge)
