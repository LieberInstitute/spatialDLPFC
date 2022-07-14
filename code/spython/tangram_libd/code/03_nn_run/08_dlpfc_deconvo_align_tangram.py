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
from PIL import Image

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

#   Variable name in ad_sp.obs containing cell counts
cell_count_var = 'cell_count'

#   Variable name in ad_sc.obs representing cell type
cell_type_var = 'cellType'

plot_file_type = 'png' # 'pdf' is also supported for higher-quality plots

################################################################################
#   Deconvolution
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
    os.path.join(processed_dir, 'ad_sp_orig_{}.h5ad'.format(sample_name))
)

ad_sc = sc.read_h5ad(
    os.path.join(processed_dir, 'ad_sc_{}.h5ad'.format(sample_name))
)

SCALE_FACTOR = ad_sp.uns['spatial'][sample_name]['scalefactors']['tissue_hires_scalef']

#-------------------------------------------------------------------------------
#   Align in "deconvolution mode"
#-------------------------------------------------------------------------------

print('Re-aligning in "deconvolution mode"...')
gpu_index = os.environ['CUDA_VISIBLE_DEVICES']
ad_map = tg.map_cells_to_space(
    ad_sc,
    ad_sp,
    mode = "constrained",
    target_count = ad_sp.obs[cell_count_var].sum(),
    density_prior = np.array(ad_sp.obs[cell_count_var]) / ad_sp.obs[cell_count_var].sum(),
    num_epochs = 1000,
    device = "cuda:" + gpu_index
)

print('Projecting annotation and plotting AUC...')
tg.project_cell_annotations(ad_map, ad_sp, annotation=cell_type_var)
annotation_list = list(pd.unique(ad_sc.obs[cell_type_var]))
tg.plot_cell_annotation_sc(
    ad_sp, annotation_list, perc=0.02, scale_factor = SCALE_FACTOR
)
f = plt.gcf()
f.savefig(
    os.path.join(
        plot_dir,
        'cell_annotation_deconvo_{}.{}'.format(sample_name, plot_file_type)
    ),
    bbox_inches='tight'
)

ad_ge = tg.project_genes(adata_map=ad_map, adata_sc=ad_sc)
df_all_genes = tg.compare_spatial_geneexp(ad_ge, ad_sp, ad_sc)
tg.plot_auc(df_all_genes)
f = plt.gcf()
f.savefig(
    os.path.join(
        plot_dir, 'deconvo_test_auc_{}.{}'.format(sample_name, plot_file_type)
    ),
    bbox_inches='tight'
)

#-------------------------------------------------------------------------------
#   Format segmentation results, form new AnnData, and plot
#-------------------------------------------------------------------------------

print('Formatting segmentation results...')
tg.create_segment_cell_df(ad_sp)

tg.count_cell_annotations(
    ad_map,
    ad_sc,
    ad_sp,
    annotation=cell_type_var,
)

ad_segment = tg.deconvolve_cell_annotations(ad_sp)

#   Produce the main deconvolution plot of interest
print('Producing main deconvolution plot...')
fig, ax = plt.subplots(1, 1, figsize=(20, 20))
sc.pl.spatial(
    ad_segment,
    color="cluster",
    size=0.6,
    show=False,
    frameon=False,
    alpha_img=0.5,
    legend_fontsize=20,
    ax=ax
)
f = plt.gcf()
f.savefig(
    os.path.join(
        plot_dir, 'deconvo_cells_{}.{}'.format(sample_name, plot_file_type)
    ),
    bbox_inches='tight'
)

#   Save all AnnDatas that were produced or modified
print('Saving AnnDatas...')
ad_map.write_h5ad(
    os.path.join(processed_dir, 'ad_map_deconvo_{}.h5ad'.format(sample_name))
)
ad_ge.write_h5ad(
    os.path.join(processed_dir, 'ad_ge_deconvo_{}.h5ad'.format(sample_name))
)
ad_sp.write_h5ad(
    os.path.join(processed_dir, 'ad_sp_aligned_deconvo_{}.h5ad'.format(sample_name))
)
ad_segment.write_h5ad(
    os.path.join(processed_dir, 'ad_segment_{}.h5ad'.format(sample_name))
)

print('Done all tasks.')
