#  Run Tangram on DLPFC subset of Matt's Neuron snRNA-seq samples coupled with
#  the 12 NN Visium samples available from spatialLIBD. This code is roughly
#  based off of the tutorial found here:
#
#  https://github.com/broadinstitute/Tangram/blob/master/tutorial_tangram_with_squidpy.ipynb
#
#  This python script will be invoked by an array of shell scripts. Some plots
#  are generated for the same marker genes mentioned by Kristen, though since
#  this is a "real run", these markers are allowed to be in the set of training
#  genes (to produce the best fit).

import os, sys
import getopt
import pyhere
from pathlib import Path

import scanpy as sc
import squidpy as sq
import numpy as np
import pandas as pd
from anndata import AnnData
import pathlib
import matplotlib.pyplot as plt
import matplotlib as mpl
import skimage
import seaborn as sns
import tangram as tg
from PIL import Image

plot_dir = pyhere.here('tangram_libd', 'plots', '03_nn_run', 'DLPFC')
out_dir = pyhere.here(
    'tangram_libd', 'processed-data', '03_nn_run', 'tangram_out_DLPFC'
)
Path(plot_dir).mkdir(parents=True, exist_ok=True)
Path(out_dir).mkdir(parents=True, exist_ok=True)

sc_path = pyhere.here(
    'tangram_libd', 'processed-data', '03_nn_run', 'sce_DLPFC.h5ad'
)
sp_path = pyhere.here(
    'tangram_libd', 'processed-data', '03_nn_run', 'visium_DLPFC.h5ad'
)
marker_path = pyhere.here(
    'tangram_libd', 'processed-data', '03_nn_run', 'pan_markers.txt'
)
sample_path = pyhere.here(
    'tangram_libd', 'processed-data', '03_nn_run', 'brain_samples.txt'
)

#  Genes we want to plot predicted vs. actual expression for
select_genes_names = [
    'SNAP25', 'MBP', 'PCP4', 'CCK', 'RORB', 'ENC1', 'CARTPT', 'NR4A2', 'RELN'
]

print('Using tangram version:', tg.__version__)

#  Grab the full list of sample names we will subset from
with open(sample_path, 'r') as f:
    sample_names = f.read().splitlines()

#  Determine this particular sample name
sample_name = sample_names[int(os.environ['SGE_TASK_ID']) - 1]
print('Subsetting to just sample {}.'.format(sample_name))

#  Load AnnDatas and list of marker genes
ad_sp = sc.read_h5ad(sp_path)
ad_sp = ad_sp[ad_sp.obs['sample_id'] == sample_name, :]
ad_sc = sc.read_h5ad(sc_path)

with open(marker_path, 'r') as f:
    markers = f.read().splitlines()

#  Use Ensembl IDs as gene names for both AnnDatas
ad_sc.var.index = ad_sc.var.gene_id
    
#  Note when genes of interest are present in the training set  
select_genes = ad_sp.var.gene_id[ad_sp.var.gene_name.isin(select_genes_names)]
assert len(select_genes_names) == len(select_genes)

for i in range(len(select_genes)):
    gene = select_genes[i]
    gene_name = select_genes_names[i]
    
    #  Verify genes of interest were measured in the experiment
    assert gene in ad_sp.var.gene_id.values
    assert gene in ad_sc.var.gene_id.values
    
    if gene in markers:
        print('Gene', gene, '(' + gene_name + ') is in the training set.')
    else:
        print('Gene', gene, '(' + gene_name + ') is in the test set.')

tg.pp_adatas(ad_sc, ad_sp, genes=markers)

#  Mapping step using GPU
gpu_index = os.environ['CUDA_VISIBLE_DEVICES']
ad_map = tg.map_cells_to_space(ad_sc, ad_sp,
    mode="cells",
    density_prior='rna_count_based',
    num_epochs=500,
    device = "cuda:" + gpu_index
)

tg.project_cell_annotations(ad_map, ad_sp, annotation="cellType")
annotation_list = list(pd.unique(ad_sc.obs['cellType']))

#  Plot spatial expression by cell-type label
tg.plot_cell_annotation_sc(
    ad_sp, annotation_list, x='pxl_row_in_fullres', y='pxl_col_in_fullres',
    perc=0.02
)
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

#  Project all cells based on trained mapping, and save the result
ad_ge = tg.project_genes(adata_map=ad_map, adata_sc=ad_sc)
ad_ge.write_h5ad(os.path.join(out_dir, 'ad_ge_' + sample_name + '.h5ad'))

#  Plot 5 lowest-scoring training genes
genes = ad_map.uns['train_genes_df']['train_score'][-5:]
ad_map.uns['train_genes_df'].loc[genes]
tg.plot_genes_sc(genes, adata_measured=ad_sp, adata_predicted=ad_ge, perc=0.02)
f = plt.gcf()
f.savefig(
    os.path.join(
        plot_dir, 'mapped_low_scoring_train_genes_' + sample_name + '.png'
    ),
    bbox_inches='tight'
)

#  Plot genes not present in spatial data
genes=['loc102633833', 'gm5700', 'gm8292'] # adjust for our data
tg.plot_genes_sc(genes, adata_measured=ad_sp, adata_predicted=ad_ge, perc=0.02)
f = plt.gcf()
f.savefig(
    os.path.join(
        plot_dir, 'mapped_missing_visium_genes_' + sample_name + '.png'
    ),
    bbox_inches='tight'
)

#  Compute average cosine similarity for test genes
df_all_genes = tg.compare_spatial_geneexp(ad_ge, ad_sp, ad_sc)
test_score = np.mean(df_all_genes.score[np.logical_not(df_all_genes.is_training)])
print('Average test score:', round(float(test_score), 4))

#  Compute average cosine similarity for training genes
train_score = np.mean(df_all_genes.score[df_all_genes.is_training])
print('Average training score:', round(float(train_score), 4))

tg.plot_auc(df_all_genes)
f = plt.gcf()
f.savefig(
    os.path.join(plot_dir, 'test_auc_' + sample_name + '.png'),
    bbox_inches='tight'
)

#   Save the spatial AnnData, which was modified to include additional data in
#   the above steps
ad_sp.write_h5ad(os.path.join(out_dir, 'ad_sp_' + sample_name + '.h5ad'))

################################################################################
#   Deconvolution
################################################################################

#   Read histology image in as numpy array
img_arr = np.array(
    Image.open(
        str(
            pyhere.here(
                "spagcn/raw-data/02-our_data_tutorial/" + sample_name + ".tif"
            )
        )
    )
)

img_arr = img_arr[:1000, :1000, :]

#   Convert to squidpy ImageContainer
img = sq.im.ImageContainer(img_arr)

#   Apply smoothing and compute segmentation masks
sq.im.process(img=img, layer="image", method="smooth")
sq.im.segment(
    img=img,
    layer="image_smooth",
    method="watershed",
    channel=0,
)

#-------------------------------------------------------------------------------
#   Visualize segmentation results
#-------------------------------------------------------------------------------

# inset_y = 1500
# inset_x = 1700
# inset_sy = 400
# inset_sx = 500

sf = ad_sp.uns['scaleFactor'][0]

inset_y = 1000
inset_x = 1000
inset_sy = 400
inset_sx = 500

fig, axs = plt.subplots(1, 3, figsize=(30, 10))
sc.pl.spatial(
    ad_sp, color="Cluster", alpha=0.7, frameon=False, show=False, ax=axs[0], 
    title="", spot_size = 50, scale_factor = sf, img = img
)
axs[0].set_title("Clusters", fontdict={"fontsize": 20})

rect = mpl.patches.Rectangle(
    (inset_y * sf, inset_x * sf),
    width=inset_sx * sf,
    height=inset_sy * sf,
    ec="yellow",
    lw=4,
    fill=False,
)
axs[0].add_patch(rect)

axs[0].axes.xaxis.label.set_visible(False)
axs[0].axes.yaxis.label.set_visible(False)

axs[1].imshow(
    img["image"][inset_y : inset_y + inset_sy, inset_x : inset_x + inset_sx, 0, 0]
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
axs[2].imshow(crop, interpolation="none", cmap=cmap, vmin=0.001)
axs[2].grid(False)
axs[2].set_xticks([])
axs[2].set_yticks([])
axs[2].set_title("Nucleous segmentation", fontdict={"fontsize": 20});
