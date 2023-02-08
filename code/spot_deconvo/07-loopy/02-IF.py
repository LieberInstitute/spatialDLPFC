from pathlib import Path
from pyhere import here
import json
import os

import numpy as np
import pandas as pd
from rasterio import Affine

from loopy.sample import Sample
from loopy.utils.utils import remove_dupes

deconvo_tools = ['tangram', 'cell2location']
spot_diameter_m = 55e-6 # 55-micrometer diameter for Visium spot
img_channels = ['Lipofuscin', 'DAPI', 'GFAP', 'NeuN', 'OLIG2', 'TMEM119']

sample_info_path = here(
    'raw-data', 'sample_info', 'Visium_IF_DLPFC_MasterExcel_01262022.xlsx'
)

img_path = here(
    'raw-data', 'Images', 'VisiumIF', 'VistoSeg', '{}.tif'
)
json_path = here(
    'processed-data', '01_spaceranger_IF', '{}', 'outs', 'spatial',
    'scalefactors_json.json'
)
tissue_path = here(
    'processed-data', '01_spaceranger_IF', '{}', 'outs', 'spatial',
    'tissue_positions_list.csv'
)
out_dir = here('processed-data', 'spot_deconvo', '07-loopy', 'IF', '{}')

cell_types_broad = [
    "Astro", "EndoMural", "Micro", "Oligo", "OPC", "Excit", "Inhib"
]
cell_types_layer = [
    "Astro", "EndoMural", "Micro", "Oligo", "OPC", "Excit_L2_3", "Excit_L3",
    "Excit_L3_4_5", "Excit_L4", "Excit_L5", "Excit_L5_6", "Excit_L6",
    "Inhib"
]
cell_types_cart = ["astro", "micro", "neuron", "oligo", "other"]

raw_results_path = here(
    "processed-data", "spot_deconvo", "05-shared_utilities", "IF",
    "results_raw_{}.csv"
)
collapsed_results_path = here(
    "processed-data", "spot_deconvo", "05-shared_utilities", "IF",
    "results_collapsed_{}.csv"
)

#   Different sample IDs are used for different files associated with each
#   sample. Determine both forms of the sample ID for this sample and update
#   path variables accordingly
sample_info = pd.read_excel(sample_info_path)[:4]
sample_ids_img = sample_info['Slide SN #'] + '_' + sample_info['Array #']
sample_ids_spot = 'Br' + sample_info['BrNumbr'].astype(int).astype(str) + \
    '_' + pd.Series([x.split('_')[1] for x in sample_info['Br_Region']]) + \
    '_IF'

#   Subset both types of IDs to this sample only
sample_id_img = sample_ids_img[int(os.environ['SGE_TASK_ID']) - 1]
sample_id_spot = sample_ids_spot[int(os.environ['SGE_TASK_ID']) - 1]

out_dir = Path(str(out_dir).format(sample_id_spot))
json_path = Path(str(json_path).format(sample_id_spot))
img_path = Path(str(img_path).format(sample_id_img))
tissue_path = Path(str(tissue_path).format(sample_id_spot))

#   Read in the tissue positions file to get spatial coordinates. Index by
#   barcode, and only take spots overlapping tissue
tissue_positions = pd.read_csv(
    tissue_path,
    header = None,
    names = ["in_tissue", "row", "col", "y", "x"], # Note the switch of x and y
    index_col = 0
)
tissue_positions.index.name = "barcode"

tissue_positions = tissue_positions.loc[
    tissue_positions['in_tissue'] == 1, ["x", "y"]
]


#   Read in the spaceranger JSON, ultimately to calculate meters per pixel for
#   the full-resolution image
with open(json_path, 'r') as f:
    spaceranger_json = json.load(f)

m_per_px = spot_diameter_m / spaceranger_json['spot_diameter_fullres']

this_sample = Sample(name = sample_id, path = out_dir)
this_sample.add_coords(
    tissue_positions, name="coords", mPerPx=m_per_px, size=spot_diameter_m
)

for cell_group in ("broad", "layer"):
    raw_results = pd.read_csv(
        str(raw_results_path).format(cell_group), index_col = "barcode"
    )
    collapsed_results = pd.read_csv(
        str(collapsed_results_path).format(cell_group), index_col = "barcode"
    )
    
    #   Add software deconvolution results at normal resolution
    if cell_group == "broad":
        cell_types = cell_types_broad
    else:
        cell_types = cell_types_layer
    
    for cell_type in cell_types:
        for deconvo_tool in deconvo_tools:
            small_results = raw_results[
                (raw_results['sample_id'] == sample_id) &
                (raw_results['deconvo_tool'] == deconvo_tool)
            ][[cell_type]]
            
            feature_name = '_'.join(
                [cell_group, deconvo_tool, cell_type]
            ).lower()
            
            small_results = small_results.merge(
                    tissue_positions, how = "left", on = "barcode"
                ).rename(
                    {cell_type: feature_name}, axis = 1
                )
            
            this_sample.add_csv_feature(
                small_results[[feature_name]],
                name = feature_name, coordName="coords"
            )
    
    #   Add CART results (at collapsed resolution)
    for cell_type in cell_types_cart:
        small_results = collapsed_results[
            (collapsed_results['sample_id'] == sample_id) &
            (collapsed_results['deconvo_tool'] == deconvo_tools[0]) &
            (collapsed_results['obs_type'] == 'actual')
        ][[cell_type]]
        
        feature_name = '_'.join([cell_group, 'cart', cell_type])
        
        small_results = small_results.merge(
                tissue_positions, how = "left", on = "barcode"
            ).rename(
                {cell_type: feature_name}, axis = 1
            )
         
        this_sample.add_csv_feature(
            small_results[[feature_name]],
            name = feature_name, coordName="coords"
        )

#   Add the IF image for this sample
this_sample.add_image(
    tiff = img_path, channels = img_channels, scale = m_per_px
)

this_sample.write()
