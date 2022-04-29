#   While 03_sce_to_anndata_dlpfc.R ideally would produce an AnnData totally
#   ready for use with Tangram and squidpy, in reality there are a few details
#   to adjust. This code is separate from 03_sce_to_anndata_dlpfc.R because it
#   is simpler to perform in python than R.

import os, sys
import pyhere
from pathlib import Path

import scanpy as sc
import numpy as np
import pandas as pd
from anndata import AnnData
import skimage
import seaborn as sns
import tangram as tg
from PIL import Image
import json

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
Path(plot_dir).mkdir(parents=True, exist_ok=True)
Path(processed_dir).mkdir(parents=True, exist_ok=True)

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

#   Parent directory whose subdirectories include JSON files with spot sizes
#   for each sample
json_dir = '/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/10X'

#-------------------------------------------------------------------------------
#   Dataset-specific variables
#-------------------------------------------------------------------------------

#   Variable name in ad_sp.obs to color by in deconvolution-related plots
cluster_var_plots = 'Cluster'

#   Variable name in ad_sc.obs representing cell type
cell_type_var = 'cellType'

#  Genes we want to plot predicted vs. actual expression for
select_genes_names = [
    'SNAP25', 'MBP', 'PCP4', 'CCK', 'RORB', 'ENC1', 'CARTPT', 'NR4A2', 'RELN'
]

spatial_coords_names = ('pxl_row_in_fullres', 'pxl_col_in_fullres')

################################################################################
#   Preprocessing
################################################################################

#  Grab the full list of sample names we will subset from
with open(sample_path, 'r') as f:
    sample_names = f.read().splitlines()

#  Determine this particular sample name
sample_name = sample_names[int(os.environ['SGE_TASK_ID']) - 1]
print('Subsetting to just sample {}.'.format(sample_name))

#  Load AnnDatas and list of marker genes
print('Loading AnnDatas...')
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

#   Make sure some variables are categorical, which enables correct coloring of
#   some later plots
ad_sc.obs[cell_type_var] = ad_sc.obs[cell_type_var].astype('category')
ad_sp.obs[cluster_var_plots] = ad_sp.obs[cluster_var_plots].astype('category')

#   Path to JSON (from spaceranger?) including spot size for this sample
json_path = json_dir + '/' + sample_name + '/scalefactors_json.json'

with open(json_path) as f: 
    json_data = json.load(f)

#   Store scalefactors in AnnData as squidpy expects
ad_sp.uns['spatial'] = {
    sample_name: {
        'scalefactors': json_data
    }
}

#   Read in high-res image as numpy array with values in [0, 1] rather than
#   [0, 255] (note that it isn't verified the original range is [0, 255]!). 
#   Then attach to AnnData object
img_path = str(
    pyhere.here(
        "tangram_libd/raw-data/03_nn_run/hires_histology", sample_name + ".png"
    )
)

img_arr = np.array(Image.open(img_path), dtype = np.float32) / 256
assert img_arr.shape == (2000, 2000, 3), img_arr.shape

ad_sp.uns['spatial'][sample_name]['images'] = { 'hires': img_arr }

#   Squidpy expects spatial coords in a very specific format
ad_sp.obsm['spatial'] = np.array(
    list(
        zip(
            ad_sp.obs[spatial_coords_names[0]],
            ad_sp.obs[spatial_coords_names[1]]
        )
    )
)

#-------------------------------------------------------------------------------
#   Save AnnDatas
#-------------------------------------------------------------------------------

ad_sp.write_h5ad(
    os.path.join(processed_dir, 'ad_sp_orig_{}.h5ad'.format(sample_name))
)

#   While the contents of the 'ad_sc' object should be identical regardless of
#   sample, it looks like the order of some variables is random,
#   but alignment and other downstream tasks are dependent on variable
#   ordering. Therefore, while it's a bit wasteful to save many "copies" of
#   'ad_sc' as done here, it simplifies code later by avoiding several order-
#   related complications that would need manual resolution
ad_sc.write_h5ad(
    os.path.join(processed_dir, 'ad_sc_{}.h5ad'.format(sample_name))
)
