import sys

import scanpy as sc
import anndata
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

import cell2location
import scvi

from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42 # enables correct plotting of text
import seaborn as sns

import pyhere

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
sample_id_var = 'sample_id'
ensembl_id_var = 'gene_id'    # in both spatial and single-cell objects!
gene_symbol_var = 'gene_name' # in both spatial and single-cell objects!
spatial_coords_names = ('pxl_row_in_fullres', 'pxl_col_in_fullres')

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
adata_ref.var['SYMBOL'] = adata_ref.obs[gene_symbol_var]
adata_ref.var.index = adata_ref.var[ensembl_id_var]
adata_ref.var_names = adata_ref.var[ensembl_id_var]
adata_ref.var.index.name = None
adata_ref.raw.var['SYMBOL'] = adata_ref.obs[gene_symbol_var]
adata_ref.raw.var.index = adata_ref.var[ensembl_id_var]
adata_ref.raw.var.index.name = None
