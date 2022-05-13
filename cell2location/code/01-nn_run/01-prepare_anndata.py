import sys
import os

import scanpy as sc
import anndata
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.backends.backend_pdf import PdfPages

import cell2location
from cell2location.utils.filtering import filter_genes
from cell2location.models import RegressionModel
import scvi

from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42 # enables correct plotting of text
import seaborn as sns

import pyhere
from pathlib import Path
from PIL import Image
import json

################################################################################
#   Variable definitions
################################################################################

sc_path = pyhere.here(
    'tangram_libd', 'processed-data', '03_nn_run', 'sce_DLPFC.h5ad'
)
sp_path = pyhere.here(
    'tangram_libd', 'processed-data', '03_nn_run', 'visium_DLPFC.h5ad'
)

processed_dir = pyhere.here('cell2location', 'processed-data', '01-nn_run')
plot_dir = pyhere.here('cell2location', 'plots', '01-nn_run')
Path(plot_dir).mkdir(parents=True, exist_ok=True)
Path(processed_dir).mkdir(parents=True, exist_ok=True)

marker_path = pyhere.here(
    'tangram_libd', 'processed-data', '03_nn_run', 'pan_markers.txt'
)

#   Parent directory whose subdirectories include JSON files with spot sizes
#   and scale factors for each sample
json_dir = '/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/10X'

# create paths and names to results folders for reference regression and
# cell2location models
ref_run_name = f'{processed_dir}/reference_signatures'
run_name = f'{processed_dir}/cell2location_map'

#   Naming conventions used for different columns in the spatial AnnData
sample_id_var = 'sample_id'   # in spatial object only
ensembl_id_var = 'gene_id'    # in both spatial and single-cell objects
gene_symbol_var = 'gene_name' # in both spatial and single-cell objects
cell_type_var = 'cellType'    # in single-cell only
spatial_coords_names = ('pxl_row_in_fullres', 'pxl_col_in_fullres')

cell_types_to_drop = ['Macrophage', 'Mural', 'Tcell']

plot_file_type = 'png' # 'pdf' is also supported for higher-quality plots

################################################################################
#   Preprocessing
################################################################################

#  Load AnnDatas
print('Loading AnnDatas...')
adata_vis = sc.read_h5ad(sp_path)
adata_ref = sc.read_h5ad(sc_path)

adata_vis.obs['sample'] = adata_vis.obs[sample_id_var]

# rename genes to ENSEMBL
adata_vis.var['SYMBOL'] = adata_vis.var[gene_symbol_var]
adata_vis.var_names = adata_vis.var[ensembl_id_var]
adata_vis.var_names.name = None

# find mitochondria-encoded (MT) genes
adata_vis.var['MT_gene'] = [
    gene.startswith('MT-') for gene in adata_vis.var['SYMBOL']
]

# remove MT genes for spatial mapping (keeping their counts in the object).
#   Actually, we'll probably use Louise's marker genes, in which case we should
#   just confirm those don't include mitochondrial genes
adata_vis.obsm['MT'] = adata_vis[:, adata_vis.var['MT_gene'].values].X.toarray()
adata_vis = adata_vis[:, ~adata_vis.var['MT_gene'].values]

# Use ENSEMBL as gene IDs to make sure IDs are unique and correctly matched
adata_ref.var['SYMBOL'] = adata_ref.var[gene_symbol_var]
adata_ref.var.index = adata_ref.var[ensembl_id_var]
adata_ref.var_names = adata_ref.var[ensembl_id_var]
adata_ref.var.index.name = None

#   Drop rare cell types from single-cell data
adata_ref = adata_ref[
    ~ adata_ref.obs[cell_type_var].isin(cell_types_to_drop), :
]

#   Subset to specific genes
selected = filter_genes(
    adata_ref, cell_count_cutoff=5, cell_percentage_cutoff2=0.03,
    nonz_mean_cutoff=1.12
)

#   This code to save plot doesn't work! TODO
f = plt.gcf()
f.savefig(
    os.path.join(
        plot_dir, f'cell_annotation.{plot_file_type}'
    ),
    bbox_inches='tight'
)

adata_ref = adata_ref[:, selected].copy()

#-------------------------------------------------------------------------------
#   Attach hi-res images and scaleFactors to spatial AnnData
#-------------------------------------------------------------------------------

adata_vis.uns['spatial'] = {}

for id in adata_vis.obs['sample'].cat.categories:
    #   Path to JSON from spaceranger including spot size for this sample
    json_path = os.path.join(json_dir, id, 'scalefactors_json.json')

    with open(json_path) as f: 
        json_data = json.load(f)

    #   Store scalefactors in AnnData as cell2location expects
    adata_vis.uns['spatial'][id] = {
        'scalefactors': json_data
    }

    #   Read in high-res image as numpy array with values in [0, 1] rather than
    #   [0, 255] (note that it isn't verified the original range is [0, 255]!). 
    #   Then attach to AnnData object
    img_path = str(
        pyhere.here(
            "tangram_libd/raw-data/03_nn_run/hires_histology",
            id + ".png"
        )
    )

    img_arr = np.array(Image.open(img_path), dtype = np.float32) / 256
    assert img_arr.shape == (2000, 2000, 3), img_arr.shape

    adata_vis.uns['spatial'][id]['images'] = { 'hires': img_arr }

#-------------------------------------------------------------------------------
#   Attach spatialCoords to spatial AnnData
#-------------------------------------------------------------------------------

#   Squidpy expects spatial coords in a very specific format
adata_vis.obsm['spatial'] = np.array(
    list(
        zip(
            adata_vis.obs[spatial_coords_names[0]],
            adata_vis.obs[spatial_coords_names[1]]
        )
    )
)

#-------------------------------------------------------------------------------
#   Save AnnDatas
#-------------------------------------------------------------------------------

adata_vis.write_h5ad(
    os.path.join(processed_dir, 'adata_vis_orig.h5ad')
)

adata_ref.write_h5ad(
    os.path.join(processed_dir, 'adata_ref.h5ad')
)

################################################################################
#   Perform regression
################################################################################

# prepare anndata for the regression model
scvi.data.setup_anndata(
    adata = adata_ref,
    batch_key = 'Sample',
    labels_key = cell_type_var
)
scvi.data.view_anndata_setup(adata_ref)

# create and train the regression model
mod = RegressionModel(adata_ref)

# Use all data for training (validation not implemented yet, train_size=1)
mod.train(max_epochs=250, batch_size=2500, train_size=1, lr=0.002, use_gpu=True)

# plot ELBO loss history during training, removing first 20 epochs from the plot
mod.plot_history(20)
