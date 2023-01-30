from pathlib import Path
from pyhere import here

import numpy as np
import pandas as pd
from rasterio import Affine

from loopy.sample import Sample
from loopy.utils.utils import remove_dupes

deconvo_tools = ['tangram', 'cell2location']
cell_types = ["astro", "micro", "neuron", "oligo", "other"]
loopy_results_path = here(
    'processed-data', 'spot_deconvo', '05-shared_utilities', 'IF',
    'loopy_Br2720_Ant_IF.csv'
)
out_dir = here('processed-data', 'spot_deconvo', '07-loopy', 'Br2720_Ant_IF')

#   Import spot-deconvo results
loopy_csv = pd.read_csv(loopy_results_path)
loopy_csv.index = loopy_csv['barcode']
loopy_csv.rename({'barcode': 'id'}, axis = 1, inplace = True)

#   Define the Sample object and add coordinates
this_sample = Sample(name = "Br2720_Ant_IF", path = out_dir)

coords = loopy_csv.loc[
    (loopy_csv['deconvo_tool'] == deconvo_tools[0]) &
    (loopy_csv['cell_type'] == cell_types[0])
][['x', 'y']]

this_sample.add_coords(coords, name="coords", mPerPx=1e-4, size=2e-3)

#   For each deconvolution tool and cell type, add a feature (cell counts for
#   each spot)
for deconvo_tool in deconvo_tools:
    for cell_type in cell_types:
        feature_name = deconvo_tool + "_" + cell_type
        #
        small_csv = loopy_csv.loc[
            (loopy_csv['deconvo_tool'] == deconvo_tool) &
            (loopy_csv['cell_type'] == cell_type)
        ].rename({'count': feature_name + "_count"}, axis = 1)
        #
        this_sample.add_csv_feature(
            small_csv[[feature_name + "_count"]],
            name = feature_name, coordName="coords"
        )

this_sample.write()
