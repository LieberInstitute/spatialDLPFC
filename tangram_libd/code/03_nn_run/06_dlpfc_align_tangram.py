import os, sys
import pyhere

import scanpy as sc
import squidpy as sq
import numpy as np
import pandas as pd
from anndata import AnnData
import matplotlib.pyplot as plt
import matplotlib as mpl
import skimage
import seaborn as sns
import tangram as tg

################################################################################
#   Variable definitions
################################################################################

#-------------------------------------------------------------------------------
#   Paths
#-------------------------------------------------------------------------------

plot_dir = pyhere.here('tangram_libd', 'plots', '03_nn_run', 'DLPFC')
processed_dir = pyhere.here(
    'tangram_libd', 'processed-data', '03_nn_run', 'tangram_out_DLPFC'
)

sample_path = pyhere.here(
    'tangram_libd', 'processed-data', '03_nn_run', 'brain_samples.txt'
)

#-------------------------------------------------------------------------------
#   Dataset-specific variables
#-------------------------------------------------------------------------------

#   Variable name in ad_sp.obs to color by in deconvolution-related plots
cluster_var_plots = 'Cluster'

#   Variable name in ad_sc.obs representing cell type
cell_type_var = 'cellType'

spatial_coords_names = ('pxl_row_in_fullres', 'pxl_col_in_fullres')

################################################################################
#   Alignment (spatial registration)
################################################################################

print('Using tangram version:', tg.__version__)

#-------------------------------------------------------------------------------
#   Load AnnDatas
#-------------------------------------------------------------------------------

#  Grab the full list of sample names we will subset from
with open(sample_path, 'r') as f:
    sample_names = f.read().splitlines()

#  Determine this particular sample name
sample_name = sample_names[int(os.environ['SGE_TASK_ID']) - 1]
print('Only using spatial sample {}.'.format(sample_name))

ad_sp = sc.read_h5ad(
    os.path.join(processed_dir, 'ad_sp_orig_{}.h5ad'.format(sample_name)
)

ad_sc = sc.read_h5ad(
    os.path.join(processed_dir, 'ad_sc_{}.h5ad'.format(sample_name)
)

#-------------------------------------------------------------------------------
#   Align
#-------------------------------------------------------------------------------

#  Mapping step using GPU
gpu_index = os.environ['CUDA_VISIBLE_DEVICES']
ad_map = tg.map_cells_to_space(ad_sc, ad_sp,
    mode="cells",
    density_prior='rna_count_based',
    num_epochs=500,
    device = "cuda:" + gpu_index
)

print('Producing exploratory plots...')
tg.project_cell_annotations(ad_map, ad_sp, annotation=cell_type_var)
annotation_list = list(pd.unique(ad_sc.obs[cell_type_var]))

#  Plot spatial expression by cell-type label
tg.plot_cell_annotation_sc(ad_sp, annotation_list, perc = 0.02)
f = plt.gcf()
f.savefig(
    os.path.join(plot_dir, 'cell_annotation_' + sample_name + '.png'),
    bbox_inches='tight'
)

#  Plot training scores
tg.plot_training_scores(ad_map, bins=20, alpha=.5)
f = plt.gcf()
f.savefig(
    os.path.join(plot_dir, 'train_scores_' + sample_name + '.pdf'),
    bbox_inches='tight'
)

#  Project all cells based on trained mapping
ad_ge = tg.project_genes(adata_map=ad_map, adata_sc=ad_sc)

#  Plot 5 lowest-scoring training genes
genes = ad_map.uns['train_genes_df']['train_score'][-5:].index
ad_map.uns['train_genes_df'].loc[genes]
tg.plot_genes_sc(
    genes, adata_measured=ad_sp, adata_predicted=ad_ge, perc=0.02
)
f = plt.gcf()
f.savefig(
    os.path.join(
        plot_dir, 'mapped_low_scoring_train_genes_' + sample_name + '.png'
    ),
    bbox_inches='tight'
)

#  Plot genes not present in spatial data
uniq_sc_genes = ad_sc.var['gene_id'][
    ~ ad_sc.var['gene_id'].isin(ad_sp.var['gene_id'])
].index
tg.plot_genes_sc(
    uniq_sc_genes[:10], adata_measured=ad_sp, adata_predicted=ad_ge, perc=0.02
)
f = plt.gcf()
f.savefig(
    os.path.join(
        plot_dir, 'mapped_missing_visium_genes_' + sample_name + '.png'
    ),
    bbox_inches='tight'
)

#  Compute average cosine similarity for test genes
df_all_genes = tg.compare_spatial_geneexp(ad_ge, ad_sp, ad_sc)
test_score = np.mean(
    df_all_genes.score[np.logical_not(df_all_genes.is_training)]
)
print('Average test score:', round(float(test_score), 4))

#  Compute average cosine similarity for training genes
train_score = np.mean(df_all_genes.score[df_all_genes.is_training])
print('Average training score:', round(float(train_score), 4))

tg.plot_auc(df_all_genes)
f = plt.gcf()
f.savefig(
    os.path.join(plot_dir, 'orig_test_auc_' + sample_name + '.png'),
    bbox_inches='tight'
)

#   Save all AnnDatas that were produced or modified
ad_map.write_h5ad(
    os.path.join(processed_dir, 'ad_map_{}.h5ad'.format(sample_name))
)
ad_ge.write_h5ad(
    os.path.join(processed_dir, 'ad_ge_{}.h5ad'.format(sample_name))
)
ad_sp.write_h5ad(
    os.path.join(processed_dir, 'ad_sp_aligned_{}.h5ad'.format(sample_name))
)
