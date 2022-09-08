import sys
import os

import scanpy as sc
import anndata
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

import cell2location
from cell2location.utils.filtering import filter_genes

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

cell_group = "layer" # "broad" or "layer"


sc_path = pyhere.here(
    "processed-data", "spot_deconvo", "05-shared_utilities",
    "sce_" + cell_group + ".h5ad"
)
sp_path = pyhere.here(
    "processed-data", "spot_deconvo", "05-shared_utilities", "IF", "spe.h5ad"
)

processed_dir = pyhere.here(
    "processed-data", "spot_deconvo", "03-cell2location", "IF", cell_group
)
plot_dir = pyhere.here(
    "plots", "spot_deconvo", "03-cell2location", "IF", cell_group
)
Path(plot_dir).mkdir(parents=True, exist_ok=True)
Path(processed_dir).mkdir(parents=True, exist_ok=True)

#   Directory containing hires image and a JSON containing scale factors and
#   spot size for a given sample. Here '{}' will be replaced by a single
#   sample name
spaceranger_dir = pyhere.here(
    'processed-data', '01_spaceranger_IF', '{}', 'outs', 'spatial'
)

marker_path = pyhere.here(
    "processed-data", "spot_deconvo", "05-shared_utilities",
    "markers_C2L_" + cell_group + ".txt"
)

#   Naming conventions used for different columns in the spatial AnnData
sample_id_var = 'sample_id'          # in spatial object only
ensembl_id_var = 'gene_id'           # in both spatial and single-cell objects
gene_symbol_var = 'gene_name'        # in both spatial and single-cell objects

# in single-cell only
if cell_group == 'broad':
    cell_type_var = 'cellType_broad_hc'
else:
    cell_type_var = 'layer_level'

spatial_coords_names = ['pxl_col_in_fullres', 'pxl_row_in_fullres']

#   If True, important per-spot cell counts from a 'clusters.csv' file.
#   Otherwise it is assumed the object already contains cell counts that will be
#   used later
IMPORT_CELL_COUNTS = True
cell_count_var = 'cell_count' # Column in ad_sp.obs containing cell counts
cell_counts_path = pyhere.here(
    'processed-data', 'spot_deconvo', '02-cellpose', '{}', 'clusters.csv'
)

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
adata_vis.obsm['MT'] = adata_vis[:, adata_vis.var['MT_gene'].values].X.toarray()
adata_vis = adata_vis[:, ~adata_vis.var['MT_gene'].values]

#   Spatial AnnData needs unique indices. Rather than using barcode (repeated
#   for every sample), use "key" (barcode + sample ID)
adata_vis.obs_names = adata_vis.obs['key']
adata_vis.obs_names.name = None

# Use ENSEMBL as gene IDs to make sure IDs are unique and correctly matched
adata_ref.var['SYMBOL'] = adata_ref.var[gene_symbol_var]
adata_ref.var.index = adata_ref.var[ensembl_id_var]
adata_ref.var_names = adata_ref.var[ensembl_id_var]
adata_ref.var.index.name = None

#   Subset to marker genes
with open(marker_path, 'r') as f:
    selected = f.read().splitlines()

adata_ref = adata_ref[:, selected].copy()

#-------------------------------------------------------------------------------
#   Attach hi-res images and scaleFactors to spatial AnnData
#-------------------------------------------------------------------------------

adata_vis.uns['spatial'] = {}

for sample_id in adata_vis.obs['sample'].cat.categories:
    #   Path to JSON from spaceranger including spot size for this sample
    json_path = pyhere.here(
        str(spaceranger_dir).format(sample_id), 'scalefactors_json.json'
    )
    
    with open(json_path) as f: 
        json_data = json.load(f)
    
    #   Store scalefactors in AnnData as cell2location expects
    adata_vis.uns['spatial'][sample_id] = {
        'scalefactors': json_data
    }
    
    #   Read in high-res image as numpy array with values in [0, 1] rather than
    #   [0, 255], then attach to AnnData object
    img_path = str(
        pyhere.here(
            str(spaceranger_dir).format(sample_id), 'tissue_hires_image.png'
        )
    )
    img_arr = np.array(Image.open(img_path), dtype = np.float32) / 256
    
    #   Store image and scalefactors in AnnData as squidpy expects
    adata_vis.uns['spatial'][sample_id] = {
        'scalefactors': json_data,
        'images' : { 'hires' : img_arr }
    }

#-------------------------------------------------------------------------------
#   Attach spatialCoords to spatial AnnData
#-------------------------------------------------------------------------------

#   Correct how spatialCoords are stored. Currently, they are a pandas
#   DataFrame, with the columns potentially in the wrong order (depending on the
#   version of SpatialExperiment used in R). We need them as a numpy array.
adata_vis.obsm['spatial'] = np.array(
    adata_vis.obsm['spatial'][spatial_coords_names]
)

#-------------------------------------------------------------------------------
#   Add cell counts from a 'clusters.csv' file, if applicable
#-------------------------------------------------------------------------------

if IMPORT_CELL_COUNTS:
    cell_counts = pd.DataFrame()

    for sample_id in adata_vis.obs['sample'].cat.categories:
        this_path = str(cell_counts_path).format(sample_id)
        
        #   Read in counts, whose associated barcodes should match the spatial
        #   AnnData's
        cell_counts = pd.concat((cell_counts, pd.read_csv(this_path)))

    cell_counts['n_cells'].index = cell_counts['key']
    adata_vis.obs[cell_count_var] = cell_counts['n_cells']
    assert not any(adata_vis.obs[cell_count_var].isna())

#-------------------------------------------------------------------------------
#   Replace special characters in some layer groups
#-------------------------------------------------------------------------------

if cell_group == "layer":
    adata_ref.obs[cell_type_var] = pd.Series(
        [x.replace('/', '_') for x in adata_ref.obs[cell_type_var]],
        dtype = 'category', index = adata_ref.obs_names
    )

#-------------------------------------------------------------------------------
#   Save AnnDatas
#-------------------------------------------------------------------------------

if cell_group == 'broad':
    adata_vis.write_h5ad(
        os.path.join(os.path.dirname(processed_dir), 'adata_vis_orig.h5ad')
    )

adata_ref.write_h5ad(
    os.path.join(processed_dir, 'adata_ref_orig.h5ad')
)
