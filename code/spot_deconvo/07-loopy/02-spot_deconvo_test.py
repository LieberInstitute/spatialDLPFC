from pathlib import Path
from pyhere import here
import json

import numpy as np
import pandas as pd
from rasterio import Affine

from loopy.sample import Sample
from loopy.utils.utils import remove_dupes

deconvo_tools = ['tangram', 'cell2location']
cell_types = ["astro", "micro", "neuron", "oligo", "other"]
spot_diameter_m = 55e-6 # 55-micrometer diameter for Visium spot
img_channels = ['Lipofuscin', 'DAPI', 'GFAP', 'NeuN', 'OLIG2', 'TMEM119']

loopy_results_path = here(
    'processed-data', 'spot_deconvo', '05-shared_utilities', 'IF',
    'loopy_Br2720_Ant_IF.csv'
)
img_path = here(
    'raw-data', 'Images', 'VisiumIF', 'VistoSeg', 'V10B01-087_A1.tif'
)
json_path = here(
    'processed-data', '01_spaceranger_IF', 'Br2720_Ant_IF', 'outs', 'spatial',
    'scalefactors_json.json'
)
out_dir = here('processed-data', 'spot_deconvo', '07-loopy', 'Br2720_Ant_IF')

#   Import spot-deconvo results
loopy_csv = pd.read_csv(loopy_results_path)
loopy_csv.index = loopy_csv['barcode']
loopy_csv.rename(
    {'barcode': 'id', 'x': 'y', 'y': 'x'}, axis = 1, inplace = True
)

#   Define the Sample object and add coordinates
this_sample = Sample(name = "Br2720_Ant_IF", path = out_dir)

coords = loopy_csv.loc[
    (loopy_csv['deconvo_tool'] == deconvo_tools[0]) &
    (loopy_csv['cell_type'] == cell_types[0])
][['x', 'y']]

#   Read in the spaceranger JSON, ultimately to calculate meters per pixel for
#   the full-resolution image
with open(json_path, 'r') as f:
    spaceranger_json = json.load(f)

m_per_px = spot_diameter_m / spaceranger_json['spot_diameter_fullres']

this_sample.add_coords(coords, name="coords", mPerPx=m_per_px, size=spot_diameter_m)

#   For each deconvolution tool and cell type, add a feature (cell counts for
#   each spot)
for deconvo_tool in deconvo_tools:
    for cell_type in cell_types:
        feature_name = deconvo_tool + "_" + cell_type
        
        small_csv = loopy_csv.loc[
            (loopy_csv['deconvo_tool'] == deconvo_tool) &
            (loopy_csv['cell_type'] == cell_type)
        ].rename({'count': feature_name + "_count"}, axis = 1)
        
        this_sample.add_csv_feature(
            small_csv[[feature_name + "_count"]],
            name = feature_name, coordName="coords"
        )


#   Add the IF image for this sample
this_sample.add_image(tiff = img_path, channels = img_channels, scale = m_per_px)

this_sample.write()
