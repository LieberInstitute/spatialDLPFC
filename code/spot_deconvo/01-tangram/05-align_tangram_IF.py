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

#   We use some modified tangram code to circumvent the need for image
#   segmentation, instead just relying on per-spot cell counts
os.chdir(pyhere.here("code", "spot_deconvo", "01-tangram"))
import custom_tg_code as ctg

################################################################################
#   Variable definitions
################################################################################

cell_group = "layer" # "broad" or "layer"

#-------------------------------------------------------------------------------
#   Paths
#-------------------------------------------------------------------------------

plot_dir = pyhere.here("plots", "spot_deconvo", "01-tangram", "IF", cell_group)
processed_dir = pyhere.here(
    "processed-data", "spot_deconvo", "01-tangram", "IF", cell_group
)
sc_path_in = pyhere.here(processed_dir, '{}', 'ad_sc.h5ad')
sp_path_in = pyhere.here(processed_dir, '{}', 'ad_sp_orig.h5ad')
id_path = pyhere.here(
    "processed-data", "spot_deconvo", "05-shared_utilities", "IF",
    "sample_ids.txt"
)

#-------------------------------------------------------------------------------
#   Dataset-specific variables
#-------------------------------------------------------------------------------

#   Variable name in ad_sp.obs containing cell counts
cell_count_var = 'cell_count'

#   Variable name in ad_sc.obs representing cell type
if cell_group == 'broad':
    cell_type_var = 'cellType_broad_hc'
else:
    cell_type_var = 'layer_level'

plot_file_type = 'pdf'

################################################################################
#   Deconvolution
################################################################################

print('Using tangram version:', tg.__version__)

#-------------------------------------------------------------------------------
#   Load AnnDatas
#-------------------------------------------------------------------------------

#  Grab the full list of sample names we will subset from
with open(id_path, 'r') as f:
    sample_names = f.read().splitlines()

#  Determine this particular sample name
sample_name = sample_names[int(os.environ['SGE_TASK_ID']) - 1]
print('Only using spatial sample {}.'.format(sample_name))

ad_sp = sc.read_h5ad(str(sp_path_in).format(sample_name))
ad_sc = sc.read_h5ad(str(sc_path_in).format(sample_name))

SCALE_FACTOR = ad_sp.uns['spatial'][sample_name]['scalefactors']['tissue_hires_scalef']

#-------------------------------------------------------------------------------
#   Align in "deconvolution mode"
#-------------------------------------------------------------------------------

print('Re-aligning in "deconvolution mode"...')
ad_map = tg.map_cells_to_space(
    ad_sc,
    ad_sp,
    mode = "constrained",
    target_count = ad_sp.obs[cell_count_var].sum(),
    density_prior = np.array(ad_sp.obs[cell_count_var]) / ad_sp.obs[cell_count_var].sum(),
    num_epochs = 1000,
    device = "cuda:0"
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
        plot_dir, sample_name,
        'cell_annotation_deconvo.{}'.format(plot_file_type)
    ),
    bbox_inches='tight'
)

ad_ge = tg.project_genes(adata_map=ad_map, adata_sc=ad_sc)
df_all_genes = tg.compare_spatial_geneexp(ad_ge, ad_sp, ad_sc)
tg.plot_auc(df_all_genes)
f = plt.gcf()
f.savefig(
    os.path.join(
        plot_dir, sample_name, 'deconvo_test_auc.{}'.format(plot_file_type)
    ),
    bbox_inches='tight'
)

#-------------------------------------------------------------------------------
#   Write 'clusters.csv' file containing counts of each cell type per spot
#-------------------------------------------------------------------------------

ctg.custom_count_cell_annotations(
    ad_map,
    ad_sc,
    ad_sp,
    annotation=cell_type_var,
    cell_count_var=cell_count_var
)

clusters = ad_sp.obsm['tangram_ct_count'].iloc[:, 2:]
clusters.index += '_' + sample_name
clusters.index.name = 'key'
clusters.rename(columns={'cell_n': cell_count_var}, inplace=True)
clusters.to_csv(os.path.join(processed_dir, sample_name, 'clusters.csv'))

# ad_segment = tg.deconvolve_cell_annotations(ad_sp)

# #   Produce the main deconvolution plot of interest
# print('Producing main deconvolution plot...')
# fig, ax = plt.subplots(1, 1, figsize=(20, 20))
# sc.pl.spatial(
#     ad_segment,
#     color="cluster",
#     size=0.6,
#     show=False,
#     frameon=False,
#     alpha_img=0.5,
#     legend_fontsize=20,
#     ax=ax
# )
# f = plt.gcf()
# f.savefig(
#     os.path.join(
#         plot_dir, 'deconvo_cells_{}.{}'.format(sample_name, plot_file_type)
#     ),
#     bbox_inches='tight'
# )

#   Save all AnnDatas that were produced or modified
print('Saving AnnDatas...')
ad_map.write_h5ad(
    os.path.join(processed_dir, sample_name, 'ad_map.h5ad')
)
ad_ge.write_h5ad(
    os.path.join(processed_dir, sample_name, 'ad_ge.h5ad')
)
ad_sp.write_h5ad(
    os.path.join(processed_dir, sample_name, 'ad_sp_aligned.h5ad')
)
# ad_segment.write_h5ad(
#     os.path.join(processed_dir, 'ad_segment_{}.h5ad'.format(sample_name))
# )

print('Done all tasks.')
