from pathlib import Path
from pyhere import here
import json
import os
import scanpy as sc

import numpy as np
import pandas as pd
from rasterio import Affine

from loopy.sample import Sample
from loopy.utils.utils import remove_dupes

deconvo_tools = ['tangram', 'cell2location']
spot_diameter_m = 55e-6 # 55-micrometer diameter for Visium spot

sample_info_path = here("code", "spaceranger", "spaceranger_parameters.txt")

cell_types_broad = [
    "Astro", "EndoMural", "Micro", "Oligo", "OPC", "Excit", "Inhib"
]
cell_types_layer = [
    "Astro", "EndoMural", "Micro", "Oligo", "OPC", "Excit_L2_3", "Excit_L3",
    "Excit_L3_4_5", "Excit_L4", "Excit_L5", "Excit_L5_6", "Excit_L6",
    "Inhib"
]

raw_results_path = here(
    "processed-data", "spot_deconvo", "05-shared_utilities", "nonIF",
    "results_raw_{}.csv"
)

out_dir = here('processed-data', 'spot_deconvo', '07-loopy', 'nonIF', '{}')

spe_path = here(
    "processed-data", "spot_deconvo", "05-shared_utilities", "nonIF", "spe.h5ad"
)
marker_broad_path = here(
    "processed-data", "spot_deconvo", "05-shared_utilities", "markers_broad.txt"
)

#   We'll include expression for the top 5 broad markers per cell type, and this
#   additional list (of classical DLPFC layer markers)
other_genes = [
    "AQP4", "HPCAL1", "CUX2", "RORB", "PCP4", "KRT17", "SNAP25", "MOBP"
]
n_markers_per_type = 5
n_markers_per_type_orig = 25 # used for the markers in 'marker_broad_path'

#   Read in sample info and subset to the 30 IDs used for analysis
sample_info = pd.read_csv(
    sample_info_path,
    names = [
        "sample_id_long", "image_id?", "position", "img_path", "img_json_path",
        "fastq"
    ],
    sep = "\t"
)

sample_info['sample_id_short'] = [
    x.split('_')[1] + '_' + x.split('_')[2]
    for x in sample_info['sample_id_long']
]
sample_info.index = sample_info['sample_id_short']

#   Add spaceranger dir
sample_info['spaceranger_dir'] = here(
    "processed-data", "rerun_spaceranger", sample_info['sample_id_long'],
    "outs", "spatial"
)

sample_id = sample_info['sample_id_short'].iloc[int(os.environ['SGE_TASK_ID']) - 1]

out_dir = Path(str(out_dir).format(sample_id))

#   Verify all directories/ files exist for all 30 samples
assert all(sample_info['spaceranger_dir'].apply(os.path.exists))
assert all(sample_info['img_path'].apply(os.path.exists))

#   Read in the spaceranger JSON, ultimately to calculate meters per pixel for
#   the full-resolution image
json_path = os.path.join(
    sample_info.loc[sample_id, 'spaceranger_dir'],
    'scalefactors_json.json'
)
with open(json_path, 'r') as f:
    spaceranger_json = json.load(f)

m_per_px = spot_diameter_m / spaceranger_json['spot_diameter_fullres']

################################################################################
#   Add software deconvolution results (at broad and layer resolution)
################################################################################

results_list = []
for cell_group in ("broad", "layer"):
    raw_results = pd.read_csv(
        str(raw_results_path).format(cell_group), index_col = "barcode"
    )
    
    if cell_group == "broad":
        cell_types = cell_types_broad
    else:
        cell_types = cell_types_layer
    #
    small_results = raw_results[
            (raw_results['sample_id'] == sample_id) &
            raw_results['deconvo_tool'].isin(deconvo_tools)
        ].melt(
            id_vars = 'deconvo_tool', value_vars = cell_types,
            ignore_index = False, var_name = 'cell_type', value_name = 'count'
        ).pivot(
            columns = ["deconvo_tool", "cell_type"], values = "count"
        )
    
    small_results.columns = [
        '_'.join([cell_group, x[0], x[1]]).lower()
        for x in small_results.columns.values
    ]
    
    results_list.append(small_results)

all_results = pd.concat(results_list, axis = 1)

################################################################################
#   Read in spatial coordinates and match up to the barcodes measured in the
#   deconvolution results
################################################################################

#   Read in the tissue positions file to get spatial coordinates. Index by
#   barcode, and only take spots overlapping tissue
tissue_positions = pd.read_csv(
    os.path.join(
        sample_info.loc[sample_id, 'spaceranger_dir'],
        'tissue_positions_list.csv'
    ),
    header = None,
    names = ["in_tissue", "row", "col", "y", "x"], # Note the switch of x and y
    index_col = 0
)
tissue_positions.index.name = "barcode"

#   Subset spatial coordinates to spots measured in the SPE object, which appear
#   to more closely match the tissue than the spaceranger tissue positions
#   subsetted to "in tissue" (they aren't identical!)
tissue_positions = all_results.merge(
    tissue_positions, how = "left", on = "barcode"
)[['x', 'y']]

################################################################################
#   Gather gene-expression data into a DataFrame to later as a feature
################################################################################

spe = sc.read(spe_path)

#   Read in broad markers, taking the top N for each cell type
with open(marker_broad_path, 'r') as f:
    markers_temp = f.read().splitlines()

markers = []
for piece in range(len(markers_temp) // n_markers_per_type_orig):
    markers += markers_temp[
        (n_markers_per_type_orig * piece):
        (n_markers_per_type_orig * piece + n_markers_per_type)
    ]

#   Convert layer-marker genes to Ensembl ID
other_genes = set(
    spe.var['gene_id'][spe.var['gene_name'].isin(other_genes)]
)

#   Take the union of the layer-marker genes with the top N broad markers
markers = set(markers).union(other_genes)

#   Subset AnnData to this sample and just the markers
spe = spe[spe.obs['sample_id'] == sample_id, spe.var['gene_id'].isin(markers)]
assert spe.shape == (all_results.shape[0], len(markers))
assert all(spe.obs.index ==  all_results.index)

#   Convert the sparse gene-expression matrix to pandas DataFrame, with the
#   gene symbols as column names
marker_df = pd.DataFrame(
    spe.X.toarray(),
    index = all_results.index,
    columns = spe.var['gene_name'][spe.var['gene_id'].isin(markers)]
)

################################################################################
#   Use the Samui API to create the importable directory for this sample
################################################################################

this_sample = Sample(name = sample_id, path = out_dir)

this_sample.add_coords(
    tissue_positions, name="coords", mPerPx=m_per_px, size=spot_diameter_m
)

#   Add spot deconvolution results as a single feature with multiple columns
this_sample.add_csv_feature(
    all_results, name = "Spot deconvolution", coordName = "coords"
)

#   Add gene expression results (multiple columns) as a feature
this_sample.add_csv_feature(
    marker_df, name = "Genes", coordName = "coords"
)

#   Add the H&E image for this sample
this_sample.add_image(
    tiff = Path(sample_info['img_path'].loc[sample_id]),
    channels = 'rgb', scale = m_per_px
)

this_sample.set_default_feature(group = "Genes", feature = "SNAP25")

this_sample.write()
