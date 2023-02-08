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
img_channels = ['Red', 'Green', 'Blue']

sample_ids_path = here(
    "processed-data", "spot_deconvo", "05-shared_utilities", "nonIF",
    "sample_ids.txt"
)
sample_info_path = here(
    "raw-data", "sample_info", "Visium_dlpfc_mastersheet.xlsx"
)

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

#   Read in sample info and subset to the 30 IDs used for analysis
sample_ids = pd.read_csv(sample_ids_path, header = None)
sample_info = pd.read_excel(sample_info_path)

sample_info.index = sample_info['sample name']
sample_info = sample_info.loc[sample_ids.squeeze()]

sample_info['sample_id_long'] = [
    x.split('/')[-1] for x in sample_info['spaceranger file path']
]

#   Update outdated path info
sample_info['spaceranger file path'] = here(
    "processed-data", "rerun_spaceranger", sample_info['sample_id_long'],
    "outs", "spatial"
)
sample_info['image file path'] = [
    x.replace(
        '/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC',
        str(here('raw-data'))
    )
    for x in sample_info['image file path']
]

sample_id = sample_ids.loc[int(os.environ['SGE_TASK_ID']) - 1][0]

out_dir = Path(str(out_dir).format(sample_id))

#   Read in the tissue positions file to get spatial coordinates. Index by
#   barcode, and only take spots overlapping tissue
tissue_positions = pd.read_csv(
    os.path.join(
        sample_info.loc[sample_id, 'spaceranger file path'],
        'tissue_positions_list.csv'
    ),
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
json_path = os.path.join(
    sample_info.loc[sample_id, 'spaceranger file path'],
    'scalefactors_json.json'
)
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

#   Add the H&E image for this sample
this_sample.add_image(
    tiff = Path(sample_info['image file path'].loc[sample_id]),
    channels = img_channels, scale = m_per_px
)

this_sample.write()
