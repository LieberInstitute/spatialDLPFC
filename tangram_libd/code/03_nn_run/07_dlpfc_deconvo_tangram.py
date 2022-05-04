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

#   Variable name in ad_sp.obs to color by in deconvolution-related plots
cluster_var_plots = 'Cluster'

#   Variable name in ad_sc.obs representing cell type
cell_type_var = 'cellType'

spatial_coords_names = ('pxl_row_in_fullres', 'pxl_col_in_fullres')

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
#   Load and preprocess full-res image
#-------------------------------------------------------------------------------

#   Read in full-resolution histology image
img_path = str(
    pyhere.here(
        "spagcn/raw-data/02-our_data_tutorial/" + sample_name + ".tif"
    )
)

#   Read histology image in as numpy array with values in [0, 65536], and
#   convert to squidpy ImageContainer
img_arr = np.array(Image.open(img_path), dtype = np.uint16) * (2 ** 8)
img = sq.im.ImageContainer(img_arr)

#   Apply smoothing and compute segmentation masks
print('Smoothing and segmenting image...')
sq.im.process(img=img, layer="image", method="smooth")
sq.im.segment(
    img=img,
    layer="image_smooth",
    method="watershed",
    channel=0,
)

img.save(os.path.join(processed_dir, sample_name + "_img.zarr"))

#-------------------------------------------------------------------------------
#   Visualize segmentation results
#-------------------------------------------------------------------------------

print('About to plot segmentation masks...')

inset_y = 6000
inset_x = 6000
inset_sy = 400
inset_sx = 500

fig, axs = plt.subplots(1, 3, figsize=(30, 10))
sc.pl.spatial(
    ad_sp, color=cluster_var_plots, alpha=0.7, frameon=False, show=False,
    ax=axs[0], title = ""
)
axs[0].set_title("Clusters", fontdict={"fontsize": 20})

rect = mpl.patches.Rectangle(
    (inset_y * SCALE_FACTOR, inset_x * SCALE_FACTOR),
    width=inset_sx * SCALE_FACTOR,
    height=inset_sy * SCALE_FACTOR,
    ec="yellow",
    lw=4,
    fill=False,
)
axs[0].add_patch(rect)

axs[0].axes.xaxis.label.set_visible(False)
axs[0].axes.yaxis.label.set_visible(False)

axs[1].imshow(
    img["image"][
        inset_y : inset_y + inset_sy, inset_x : inset_x + inset_sx, 0, 0
    ]
    / 65536,
    interpolation="none",
)
axs[1].grid(False)
axs[1].set_xticks([])
axs[1].set_yticks([])
axs[1].set_title("DAPI", fontdict={"fontsize": 20})

crop = img["segmented_watershed"][
    inset_y : inset_y + inset_sy, inset_x : inset_x + inset_sx
].values.squeeze(-1)
crop = skimage.segmentation.relabel_sequential(crop)[0]
cmap = plt.cm.plasma
cmap.set_under(color="black")

#   Why did we need to change this line?
axs[2].imshow(crop[:, : ,0], interpolation="none", cmap=cmap, vmin=0.001)
# axs[2].imshow(crop, interpolation="none", cmap=cmap, vmin=0.001)

axs[2].grid(False)
axs[2].set_xticks([])
axs[2].set_yticks([])
axs[2].set_title("Nucleous segmentation", fontdict={"fontsize": 20});
f = plt.gcf()
f.savefig(
    os.path.join(
        plot_dir, 'segmentation_test_{}.{}'.format(sample_name, plot_file_type)
    ),
    bbox_inches='tight'
)

#-------------------------------------------------------------------------------
#   Extract info about segmented nuclei
#-------------------------------------------------------------------------------

print('Extracting info about segmented nuclei...')

# define image layer to use for segmentation
features_kwargs = {
    "segmentation": {
        "label_layer": "segmented_watershed",
        "props": ["label", "centroid"],
        "channels": [1, 2],
    }
}

# calculate segmentation features
sq.im.calculate_image_features(
    ad_sp,
    img,
    layer="image",
    key_added="image_features",
    features_kwargs=features_kwargs,
    features="segmentation",
    mask_circle=True
)

ad_sp.obs["cell_count"] = ad_sp.obsm["image_features"]["segmentation_label"]
os.chdir(plot_dir)
sc.pl.spatial(
    ad_sp, color=[cluster_var_plots, "cell_count"], frameon=False,
    save = 'per_spot_cell_counts_{}.{}'.format(sample_name, plot_file_type)
)

#-------------------------------------------------------------------------------
#   Re-align in "deconvolution mode"
#-------------------------------------------------------------------------------

print('Re-aligning in "deconvolution mode"...')
gpu_index = os.environ['CUDA_VISIBLE_DEVICES']
ad_map = tg.map_cells_to_space(
    ad_sc,
    ad_sp,
    mode = "constrained",
    target_count = ad_sp.obs.cell_count.sum(),
    density_prior = np.array(ad_sp.obs.cell_count) / ad_sp.obs.cell_count.sum(),
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

print('Attempting to save segmented AnnData...')
ad_segment.write_h5ad(
    os.path.join(processed_dir, 'ad_segment_' + sample_name + '.h5ad')
)

print('Done all tasks.')
