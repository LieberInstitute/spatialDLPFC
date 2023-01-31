from pathlib import Path
from pyhere import here
import json
import os

import numpy as np
import pandas as pd
from rasterio import Affine

from loopy.sample import Sample
from loopy.utils.utils import remove_dupes

deconvo_tools = ['tangram', 'cell2location', 'CART']
cell_types = ["astro", "micro", "neuron", "oligo", "other"]
spot_diameter_m = 55e-6 # 55-micrometer diameter for Visium spot
img_channels = ['Lipofuscin', 'DAPI', 'GFAP', 'NeuN', 'OLIG2', 'TMEM119']

loopy_results_path = here(
    'processed-data', 'spot_deconvo', '05-shared_utilities', 'IF',
    'loopy_{}.csv'
)
img_path = here(
    'raw-data', 'Images', 'VisiumIF', 'VistoSeg', '{}.tif'
)
json_path = here(
    'processed-data', '01_spaceranger_IF', '{}', 'outs', 'spatial',
    'scalefactors_json.json'
)
out_dir = here('processed-data', 'spot_deconvo', '07-loopy', '{}')

sample_info_path = here(
    'raw-data', 'sample_info', 'Visium_IF_DLPFC_MasterExcel_01262022.xlsx'
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

#   Fill in sample ID for the path-related variables
img_path = Path(str(img_path).format(sample_id_img))
json_path = Path(str(json_path).format(sample_id_spot))
loopy_results_path = Path(str(loopy_results_path).format(sample_id_spot))
out_dir = Path(str(out_dir).format(sample_id_spot))

#   Import spot-deconvo results
loopy_csv = pd.read_csv(loopy_results_path)
loopy_csv.index = loopy_csv['barcode']
loopy_csv.rename({'barcode': 'id'}, axis = 1, inplace = True)

#   Define the Sample object and add coordinates
this_sample = Sample(name = sample_id_spot, path = out_dir)

coords = loopy_csv.loc[
    (loopy_csv['deconvo_tool'] == deconvo_tools[0]) &
    (loopy_csv['cell_type'] == cell_types[0])
][['x', 'y']]

#   Read in the spaceranger JSON, ultimately to calculate meters per pixel for
#   the full-resolution image
with open(json_path, 'r') as f:
    spaceranger_json = json.load(f)

m_per_px = spot_diameter_m / spaceranger_json['spot_diameter_fullres']

this_sample.add_coords(
    coords, name="coords", mPerPx=m_per_px, size=spot_diameter_m
)

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
this_sample.add_image(
    tiff = img_path, channels = img_channels, scale = m_per_px
)

this_sample.write()
