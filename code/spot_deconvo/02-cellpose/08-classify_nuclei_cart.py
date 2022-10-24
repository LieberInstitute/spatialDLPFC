#   Given mean fluorescence intensities for each nucleus, classify cell types
#   present in each spot.

import os
import matplotlib.pyplot as plt
import numpy as np
from numpy.random import default_rng
import pandas as pd
import seaborn as sns
import tifffile
from scipy import ndimage
from skimage.measure import regionprops, regionprops_table
import pyhere
from pathlib import Path
import pickle

################################################################################
#   Variable definitions
################################################################################

#-------------------------------------------------------------------------------
#   Paths
#-------------------------------------------------------------------------------

sample_info_path = pyhere.here(
    'raw-data', 'sample_info', 'Visium_IF_DLPFC_MasterExcel_01262022.xlsx'
)
img_path = pyhere.here('raw-data', 'Images', 'VisiumIF', 'VistoSeg', '{}.tif')
expanded_mask_path = pyhere.here(
    'processed-data', 'spot_deconvo', '02-cellpose', '{}', 'expanded_masks.npy'
)
spot_path = pyhere.here(
    'processed-data', '01_spaceranger_IF', '{}', 'outs', 'spatial',
    'tissue_positions_list.csv'
)
df_path = pyhere.here(
    'processed-data', 'spot_deconvo', '02-cellpose', '{}', 'df.csv'
)

#   Trained DecisionTreeClassifier path
model_path = pyhere.here(
    'processed-data', 'spot_deconvo', '02-cellpose', 'decision_tree.pkl'
)

#   Main output: rows are spots and columns are cell types (values are counts)
clusters_path = pyhere.here(
    'processed-data', 'spot_deconvo', '02-cellpose', '{}', 'clusters.csv'
)

#   Secondary output: rows are cells and columns are metrics/info
cells_path = pyhere.here(
    'processed-data', 'spot_deconvo', '02-cellpose', '{}', 'cell_metrics.csv'
)

plot_dir = pyhere.here("plots", "spot_deconvo", "02-cellpose", "cart")

Path(plot_dir).mkdir(parents=True, exist_ok=True)

#-------------------------------------------------------------------------------
#   Dataset-specific variables
#-------------------------------------------------------------------------------

names = {0: "junk", 1: "dapi", 2: "gfap", 3: "neun", 4: "olig2", 5: "tmem119"}
cell_types = {
    "neun": "neuron",
    "olig2": "oligo",
    "tmem119": "micro",
    "gfap": "astro"
}

plot_file_type = 'png' # 'pdf' is also supported for higher-quality plots

################################################################################
#   Functions
################################################################################

def plot_roi(img, props, indices, vmax: int = 128, pad: int = 25):
    #   Set up the plot
    fig, axs = plt.subplots(
        nrows=len(indices), ncols=6, figsize=(18, len(indices) * 3)
    )
    axs = axs.flatten()
    
    #   Loop through each nucleus, each of which will be a row in the final plot
    j = 0
    for idx in indices:
        bbox = props[idx]["bbox"]
        roi = props[idx]["image"]
        
        axs[j].imshow(np.pad(roi, (pad, pad), mode="constant", constant_values=0), aspect="equal")
        axs[j].grid(False)
        if j < 6:
            axs[j].set_title("Dilated mask")
        
        j += 1
        
        for i in range(1, 6):
            axs[j].imshow(
                img[
                    i,
                    max(0, bbox[0] - pad) : min(img.shape[1], bbox[2] + pad),
                    max(0, bbox[1] - pad) : min(img.shape[2], bbox[3] + pad),
                ],
                vmax=vmax,
                aspect="equal",
            )
            
            if j < 6:
                axs[j].set_title(names[i])
            
            axs[j].grid(False)
            j += 1
    
    return fig

################################################################################
#   Analysis
################################################################################

# os.environ['SGE_TASK_ID'] = '1'

rng = default_rng()

img_path = pyhere.here('raw-data', 'Images', 'VisiumIF', 'VistoSeg', '{}.tif')
expanded_mask_path = pyhere.here(
    'processed-data', 'spot_deconvo', '02-cellpose', '{}', 'expanded_masks.npy'
)
spot_path = pyhere.here(
    'processed-data', '01_spaceranger_IF', '{}', 'outs', 'spatial',
    'tissue_positions_list.csv'
)
df_path = pyhere.here(
    'processed-data', 'spot_deconvo', '02-cellpose', '{}', 'df.csv'
)

clusters_path = pyhere.here(
    'processed-data', 'spot_deconvo', '02-cellpose', '{}', 'clusters.csv'
)
cells_path = pyhere.here(
    'processed-data', 'spot_deconvo', '02-cellpose', '{}', 'cell_metrics.csv'
)

#-------------------------------------------------------------------------------
#   Read in sample info and adjust paths for this particular sample ID
#-------------------------------------------------------------------------------

#   Different sample IDs are used for different files associated with each
#   sample. Determine both forms of the sample ID for this sample and update
#   path variables accordingly
sample_info = pd.read_excel(sample_info_path)[:4]
sample_ids_img = sample_info['Slide SN #'] + '_' + sample_info['Array #']
sample_ids_spot = 'Br' + sample_info['BrNumbr'].astype(int).astype(str) + \
    '_' + pd.Series([x.split('_')[1] for x in sample_info['Br_Region']]) + \
    '_IF'

sample_id_img = sample_ids_img[int(os.environ['SGE_TASK_ID']) - 1]
img_path = str(img_path).format(sample_id_img)

sample_id_spot = sample_ids_spot[int(os.environ['SGE_TASK_ID']) - 1]
spot_path = str(spot_path).format(sample_id_spot)
df_path = str(df_path).format(sample_id_spot)
clusters_path = str(clusters_path).format(sample_id_spot)
cells_path = str(cells_path).format(sample_id_spot)
expanded_mask_path = str(expanded_mask_path).format(sample_id_spot)

Path(clusters_path).parents[0].mkdir(parents=True, exist_ok=True)

#-------------------------------------------------------------------------------
#   Read in spot data and other inputs
#-------------------------------------------------------------------------------

#   Read in all spot data
raw = pd.read_csv(
    spot_path,
    header=None,
    names=["barcode", "included", "row", "col", "x", "y"],
)

#   Take only spots that overlap tissue
raw = raw.iloc[raw.included[raw.included == 1].index].reset_index().drop(
    columns=["included", "index"]
)

imgs = tifffile.imread(img_path)

df = pd.read_csv(df_path)
df.rename({'Unnamed: 0': 'id'}, axis = 1, inplace = True)
df.index = df['id']
x = df.drop(['id', 'x', 'y', 'dist', 'idx'], axis = 1)

expanded_masks = np.load(expanded_mask_path)

with open(model_path, 'rb') as f:
    model = pickle.load(f)

#-------------------------------------------------------------------------------
#   Call cell types (classify nuclei)
#-------------------------------------------------------------------------------

df['cell_type'] = model.predict(x)

#   Ensure cell-type names match with what we expect
assert set(df['cell_type']) == set(cell_types.values())

#-------------------------------------------------------------------------------
#   Visually verify cell-type calls
#-------------------------------------------------------------------------------

props = regionprops(expanded_masks)
examples_per_type = 5

for cell_type in df['cell_type'].unique():
    #   Randomly pick 5 distinct rows for this cell type
    indices = rng.choice(
        df[df['cell_type'] == cell_type].index,
        examples_per_type, replace = False
    )
    
    #   Print numeric intensities for these cells
    print(f'Intensities for {examples_per_type} random {cell_type} cells:')
    print(df.loc[indices][cell_types.keys()])
    
    #   Plot intensities
    fig = plot_roi(imgs, props, indices)
    plt.suptitle(f'Cells classified as {cell_type}')
    fig.savefig(
        os.path.join(
            plot_dir, f'{cell_type}_{sample_id_img}.{plot_file_type}'
        )
    )
    plt.close('all')

#-------------------------------------------------------------------------------
#   Count cells per spot and print some related statistics
#-------------------------------------------------------------------------------

#   Count cell types in each spot
for cell_type in df['cell_type'].unique():
    raw[cell_type] = df[
        df['cell_type'] == cell_type
    ].groupby('idx')['idx'].count().astype(int)

#   Some spots will have no cells (and NaN values for this reason). Also count
#   the total number of cells per spot
raw.fillna(0, inplace=True)
raw['n_cells'] = raw[[x for x in df['cell_type'].unique()]].sum(axis=1)

#   Correct the data type
for column_name in list(df['cell_type'].unique()) + ['n_cells']:
    raw[column_name] = raw[column_name].astype(int)

#   Print number of cells of each type
for cell_type in df['cell_type'].unique():
    print(f'Total number of {cell_type} cells: {raw[cell_type].sum()}')


#   Print some quick stats about the cell counts per spot
prop_nonzero = round(100 * np.count_nonzero(raw['n_cells']) / raw.shape[0], 1)
print(f"Percentage of spots with at least one cell: {prop_nonzero}%")
print(f"Mean number of cells per spot: {round(np.mean(raw['n_cells']), 3)}")
print(f"Standard deviation (num cells per spot): {round(np.std(raw['n_cells']), 2)}")
print(f"Max number of cells per spot: {np.max(raw['n_cells'])}")

#-------------------------------------------------------------------------------
#   Export the table of metrics for each cell
#-------------------------------------------------------------------------------

#   Use more informative column names for output
df.rename(
    columns = {
        'area': 'nucleus_area', 'x': 'centroid_x', 'y': 'centroid_y',
        'dist': 'dist_to_nearest_spot'
    },
    inplace = True
)

df.rename(
    columns = {
        x: f'{x}_intensity' for x in cell_types.keys()
    },
    inplace = True
)

df.drop(['idx'], axis = 1, inplace = True)

#   Save
df.to_csv(cells_path, float_format="%.3f")

#-------------------------------------------------------------------------------
#   Export spot-level table as a 'clusters.csv' file
#-------------------------------------------------------------------------------

#   Make compatible with spatialLIBD 'clusters.csv' format
raw.index = raw['barcode'] + '_' + sample_id_spot
raw.index.name = 'key'
raw.drop(['row', 'col', 'x', 'y', 'barcode'], axis = 1, inplace = True)

#   Save
raw.to_csv(clusters_path, float_format="%.3f")
