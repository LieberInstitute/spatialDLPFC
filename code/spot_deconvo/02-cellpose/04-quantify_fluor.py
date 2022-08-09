#   Given masks from cellpose (from segmenting the DAPI channel of Visium IF
#   images), quantify mean fluorescence in each non-DAPI channel within each
#   nucleus (dilated to include a region around each nucleus), and save a pandas
#   DataFrame with these values for each nucleus.

import os
import matplotlib.pyplot as plt
import numpy as np
from numpy.random import default_rng
import pandas as pd
import seaborn as sns
import tifffile
from scipy.spatial import KDTree
from scipy import ndimage
from skimage.measure import regionprops, regionprops_table
import pyhere
from pathlib import Path
import json

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
mask_path = pyhere.here(
    'processed-data', 'spot_deconvo', '02-cellpose', 'masks', '{}_mask.npy'
)
spot_path = pyhere.here(
    'processed-data', '01_spaceranger_IF', '{}', 'outs', 'spatial',
    'tissue_positions_list.csv'
)
scale_path = pyhere.here(
    'processed-data', '01_spaceranger_IF', '{}', 'outs', 'spatial',
    'scalefactors_json.json'
)
out_df_path = pyhere.here(
    'processed-data', 'spot_deconvo', '02-cellpose', '{}', 'df.csv'
)
out_masks_path = pyhere.here(
    'processed-data', 'spot_deconvo', '02-cellpose', '{}', 'expanded_masks.npy'
)
plot_dir = pyhere.here("plots", "spot_deconvo", "02-cellpose", "count_cells_no_GFAP")

Path(plot_dir).mkdir(parents=True, exist_ok=True)

#-------------------------------------------------------------------------------
#   Dataset-specific variables
#-------------------------------------------------------------------------------

names = {0: "junk", 1: "dapi", 2: "gfap", 3: "neun", 4: "olig2", 5: "tmem119"}
cell_types = {
    "neun": "neuron",
    "olig2": "oligo",
    "tmem119": "micro"
}

#   Nucleus area, in number of pixels, below which a cell is ignored. Tested by
#   eye
area_threshold = 60

dilation_radius = 15
dilation_chunk_size = 6

plot_file_type = 'png' # 'pdf' is also supported for higher-quality plots

################################################################################
#   Functions
################################################################################

#   Perform a dilation-like transformation on 'img' by total size
#   'dilation_radius'. The whole transformation actually involves a series
#    of dilations by size [chunk_size / 2] followed by inversion. The idea
#   is to expand each mask equally without bias towards masks with larger-valued
#   labels (in case of overlap, which frequently happens)
def balanced_dilation(img, dilation_radius, chunk_size, verbose = False):
    assert chunk_size % 2 == 0, 'chunk_size must be even'
    assert 2 * dilation_radius % chunk_size == 0, 'cannot break this radius into chunks'
    
    num_chunks =  int(2 * dilation_radius / chunk_size)
    dilation_chunked = int(dilation_radius / num_chunks)
    assert num_chunks % 2 == 1, 'must use an odd number of chunks'
    
    #   We'll use -1 * MAX_VALUE as a placeholder, assumed to be smaller than
    #   all elements in 'img'. Check this assumption
    MAX_VALUE = 2 ** 31 - 1
    assert(np.all(img < MAX_VALUE))
    
    expanded_masks = img.copy().astype(np.int32)
    
    for i in range(num_chunks):
        if verbose:
            print(f'Dilating by {dilation_chunked} pixels...')
        
        #   Make sure zero-valued elements are always treated as the smallest
        #   value possible for dilation
        zero_indices = expanded_masks == 0
        expanded_masks[zero_indices] = -1 * MAX_VALUE
        
        expanded_masks = ndimage.grey_dilation(
            expanded_masks,
            size = (dilation_chunked, dilation_chunked)
        )
        
        #   Return "zero-valued" elements to a true value of 0
        zero_indices = expanded_masks == -1 * MAX_VALUE
        expanded_masks[zero_indices] = 0
        
        if i < num_chunks - 1:
            if verbose:
                print('Inverting...')
            
            expanded_masks *= -1
    
    return expanded_masks.astype(img.dtype)

################################################################################
#   Analysis
################################################################################

# os.environ['SGE_TASK_ID'] = '1'

rng = default_rng()

#-------------------------------------------------------------------------------
#   Read in sample info and adjust paths for this particular sample ID
#-------------------------------------------------------------------------------

#   Different sample IDs are used for different files associated with each
#   sample. Determine both forms of the sample ID for this sample and update
#   path variables accordingly
sample_info = pd.read_excel(sample_info_path, header = 1)[:4]
sample_ids_img = sample_info['Slide SN #'] + '_' + sample_info['Array #']
sample_ids_spot = 'Br' + sample_info['BrNumbr'].astype(int).astype(str) + \
    '_' + pd.Series([x.split('_')[1] for x in sample_info['Br_Region']]) + \
    '_IF'

sample_id_img = sample_ids_img[int(os.environ['SGE_TASK_ID']) - 1]
img_path = str(img_path).format(sample_id_img)
mask_path = str(mask_path).format(sample_id_img)

sample_id_spot = sample_ids_spot[int(os.environ['SGE_TASK_ID']) - 1]
spot_path = str(spot_path).format(sample_id_spot)
out_df_path = str(out_df_path).format(sample_id_spot)
out_masks_path = str(out_masks_path).format(sample_id_spot)

Path(out_df_path).parents[0].mkdir(parents=True, exist_ok=True)

#   Path to JSON from spaceranger including spot size for this sample
json_path = str(scale_path).format(sample_id_spot)
with open(json_path) as f: 
    json_data = json.load(f)

#-------------------------------------------------------------------------------
#   Read in spot data
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

#-------------------------------------------------------------------------------
#   Quantify mean fluorescence for each channel at each nucleus
#-------------------------------------------------------------------------------

#   Load multi-channel image and masks from segmenting DAPI channel
imgs = tifffile.imread(img_path)
masks = np.load(mask_path)

#   Dilate the original masks
print(f'Dilating original masks by radius {dilation_radius} and chunk size {dilation_chunk_size}.')
expanded_masks = balanced_dilation(masks, dilation_radius, dilation_chunk_size)

#   Quantify the mean image fluorescence intensity at each nucleus
#   identified by segmenting the DAPI channel. This is done for each
#   (non-lipofuscin, non-DAPI) channel
its = {
    names[i]: regionprops_table(
        expanded_masks, intensity_image=imgs[i], properties=["intensity_mean"]
    )["intensity_mean"]
    for i in range(2, 6)
}

#   Create a table containing the centroids and areas of each mask
#   (nucleus), and add this info to the intensities table
general = regionprops_table(masks, properties=["centroid", "area"])
its["area"] = general["area"]
its["x"] = general["centroid-0"]
its["y"] = general["centroid-1"]

df = pd.DataFrame(its)

#-------------------------------------------------------------------------------
#   Exploratory plot: show the distribution of masks over spots
#-------------------------------------------------------------------------------

#   Plot mask spatial distribution vs. spot distribution; there should be
#   quite a bit of overlap
plt.clf()
plt.scatter(raw["x"], raw["y"], 2)
plt.scatter(df["x"], df["y"], 2)
plt.savefig(
    os.path.join(
        plot_dir, f'mask_spot_overlap_{sample_id_img}.{plot_file_type}'
    )
)

#-------------------------------------------------------------------------------
#   Filter masks that are too small or not within a spot
#-------------------------------------------------------------------------------

# Build KD tree for nearest neighbor search.
kd = KDTree(raw[["x", "y"]])

#   For each mask, assign a distance to the nearest spot ('dist') and index of
#   that spot in 'df' ('idx'). Add this info to 'df'
dist, idx = kd.query(df[["x", "y"]])
dist = pd.DataFrame({"dist": dist, "idx": idx})
df = pd.concat([df, dist], axis=1)

# Filters out masks that are smaller than a threshold and
# masks whose centroid is farther than the spot radius (aka not inside the
# spot).
frac_kept = round(100 * np.sum(df.dist < json_data['spot_diameter_fullres'] / 2) / len(df.dist), 1)
print(f'Keeping {frac_kept}% of masks, which were within spots covered by tissue.')
df = df[df.dist < json_data['spot_diameter_fullres'] / 2]

frac_kept = round(100 * np.sum(df.area > area_threshold) / len(df.area), 1)
print(f'Keeping {frac_kept}% of remaining masks, which met a minimum area threshold.')
df = df[df.area > area_threshold]

#-------------------------------------------------------------------------------
#   Save relevant data
#-------------------------------------------------------------------------------

df.to_csv(out_df_path)
np.save(out_masks_path, expanded_masks)
