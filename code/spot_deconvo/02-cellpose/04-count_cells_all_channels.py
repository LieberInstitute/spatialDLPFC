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

#   Build upon Richard's method (as seen in 02-count_cells.py) to develop a
#   more robust approach to calling cell types in Visium-IF data. In particular,
#       1. incorporate the GFAP channel and associated cell type
#       2. call cell types for each nucleus, such that the sum of cells of every
#          type is equal to the total number of cells in the spot (not true in
#          the previous approach!)

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
clusters_path = pyhere.here(
    'processed-data', 'spot_deconvo', '02-cellpose', '{}', 'clusters.csv'
)
cells_path = pyhere.here(
    'processed-data', 'spot_deconvo', '02-cellpose', '{}', 'cell_metrics.csv'
)
plot_dir = pyhere.here("plots", "spot_deconvo", "02-cellpose")

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
m_per_px = 0.497e-6
spot_radius = 65e-6
area_threshold = 200

dilation_radius = 15
dilation_chunk_size = 6

#   Mean fluorescence below this value will be considered an indication of an 
#   absence of the given cell type
noise_cutoff = 3

plot_file_type = 'png' # 'pdf' is also supported for higher-quality plots

################################################################################
#   Functions
################################################################################

def plot_roi(img, props, idx: int, vmax: int = 128, nrows: int = 3, ncols: int = 2, pad: int = 5):
    bbox = props[idx]["bbox"]
    roi = props[idx]["image"]
    fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(8, 8))
    axs = axs.flatten()
    axs[0].imshow(np.pad(roi, (pad, pad), mode="constant", constant_values=0), aspect="equal")
    axs[0].set_title("Mask")
    axs[0].grid(False)
    for i in range(1, 6):
        axs[i].imshow(
            img[
                i,
                max(0, bbox[0] - pad) : min(img.shape[1], bbox[2] + pad),
                max(0, bbox[1] - pad) : min(img.shape[2], bbox[3] + pad),
            ],
            vmax=vmax,
            aspect="equal",
        )
        axs[i].set_title(names[i])
        axs[i].grid(False)

    fig.tight_layout()
    return fig

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
clusters_path = str(clusters_path).format(sample_id_img)
cells_path = str(cells_path).format(sample_id_img)
mask_path = str(mask_path).format(sample_id_img)
sample_id_spot = sample_ids_spot[int(os.environ['SGE_TASK_ID']) - 1]
spot_path = str(spot_path).format(sample_id_spot)

Path(clusters_path).parents[0].mkdir(parents=True, exist_ok=True)

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
props = regionprops(expanded_masks)

#-------------------------------------------------------------------------------
#   Exploratory plots: plot an example ROI and the distribution of masks over
#   spots
#-------------------------------------------------------------------------------

fig = plot_roi(imgs, props, 400) # 400 here is an arbitrarily chosen index
fig.savefig(
    os.path.join(plot_dir, 'random_roi_{}.{}'.format(sample_id_img, plot_file_type)),
    bbox_inches='tight'
)

#   Plot mask spatial distribution vs. spot distribution; there should be
#   quite a bit of overlap
plt.clf()
plt.scatter(raw["x"], raw["y"], 2)
plt.scatter(df["y"], df["x"], 2)
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
px_dist = spot_radius / m_per_px  # meter per px.
df = df[df.dist < px_dist]

frac_kept = round(100 * np.sum(df.area > area_threshold) / len(df.area), 1)
print(f'Keeping {frac_kept}% of masks, which met a minimum area threshold.')
df = df[df.area > area_threshold]

#-------------------------------------------------------------------------------
#   Explore how fluorescence of each channel is distributed
#-------------------------------------------------------------------------------

#   Histograms of raw fluorescence values in each channel
plt.clf()
fig, axs = plt.subplots(nrows=4, ncols=1, figsize=(10, 10))
for i, channel in enumerate(cell_types.keys()):
    axs[i].hist(df[channel], bins = 100)
    axs[i].set_title(channel)

# plt.show()
plt.savefig(
    os.path.join(
        plot_dir, f'fluor_histograms_raw_{sample_id_img}.{plot_file_type}'
    )
)

#   Histograms of trimmed fluorescence values in each channel (above a
#   "noise cutoff")
plt.clf()
fig, axs = plt.subplots(nrows=4, ncols=1, figsize=(10, 10))
for i, channel in enumerate(cell_types.keys()):
    axs[i].hist(df[channel][df[channel] > noise_cutoff], bins = 30)
    axs[i].set_title(channel + ' (trimmed)')

plt.savefig(
    os.path.join(
        plot_dir, f'fluor_histograms_trimmed_{sample_id_img}.{plot_file_type}'
    )
)

#   Check what fraction of rows have an unclear cell type (very low
#   intensity in all channels)
temp = (df['neun'] < noise_cutoff) & \
    (df['olig2'] < noise_cutoff) & \
    (df['tmem119'] < noise_cutoff) & \
    (df['gfap'] < noise_cutoff)

#   However, some nuclei don't have a clear cell type with this cutoff
bad_perc = round(100 * np.sum(temp) / len(temp), 1) # ~5%
print(f"{bad_perc}% of nuclei don't have a clear cell type.")

#   Look at a few random examples of "bad cells". Generally, it looks like GFAP
#   and TMEM119 dominate but are still quite low in intensity. This possibly
#   shows that we aren't looking at a large enough region around the nucleus,
#   but looking at the image ROIs show a diverse number of scenarios
num_bad_examples = 5
indices = np.random.randint(len(temp[temp]), size=num_bad_examples)
for i in range(num_bad_examples):
    print(f'Intensities for "bad cell" {i}:')

    bad_index = temp[temp].index[indices[i]]
    print(df.loc[bad_index][cell_types.keys()])

    fig = plot_roi(imgs, props, bad_index)
    fig.savefig(
        os.path.join(
            plot_dir, f'bad_cell_{i+1}_{sample_id_img}.{plot_file_type}'
        )
    )

#-------------------------------------------------------------------------------
#   Call cell types
#-------------------------------------------------------------------------------

#   Filtering out what we think are cells of a non-target type dramatically
#   changes the mean intensity, which is a good sign (we get a clear signal for
#   what intensity looks like in the target cell type)
weights = {}
for channel in cell_types.keys():
    print(f'Mean for {channel} above cutoff: {np.mean(df[channel][df[channel] > noise_cutoff])}')
    print(f'Mean for {channel} otherwise: {np.mean(df[channel])}')
    weights[channel] = 1 / np.mean(df[channel][df[channel] > noise_cutoff])


#   Weight the intensity of each channel inversely by the average intensity seen
#   in cells above the "noise cutoff". Then call cell types by choosing the
#   maximal weighted intensity across the 4 possible types
temp = df[[x for x in cell_types.keys()]].copy()
for channel in cell_types.keys():
    temp[channel] *= weights[channel]

temp['max_channel'] = [temp.columns[np.argmax(temp.loc[i])] for i in temp.index]
df['cell_type'] = [cell_types[x] for x in temp['max_channel']]

#   For now, assume all cells with low fluorescence in all channels are a cell
#   type not measured by a marker/ image channel
df.loc[
    (df['neun'] < noise_cutoff) &
    (df['olig2'] < noise_cutoff) &
    (df['tmem119'] < noise_cutoff) &
    (df['gfap'] < noise_cutoff),
    'cell_type'
] = 'other'

#-------------------------------------------------------------------------------
#   Verify how well the deconvolution method performed
#-------------------------------------------------------------------------------

#   For each classified cell type, show how mean intensities look for the other
#   channels. We'd like to see a large gap between the target cell type (as
#   classified) and all others. Note that some gap is likely to arise simply by
#   the definition of classification itself: e.g. neurons will generally have
#   higher neun intensity (that's what made us call them neurons)
plt.clf()
fig, axs = plt.subplots(nrows=4, ncols=1, figsize=(10, 10))
plt.suptitle(f'Weighted mean intensity comparison')
for i, channel in enumerate(cell_types.keys()):
    a = temp[temp['max_channel'] == channel]
    data_list = [a[x] for x in cell_types.keys()]
    axs[i].boxplot(data_list)
    axs[i].set_title(f'{cell_types[channel]} cells')
    axs[i].set_xticklabels(cell_types.keys())

# plt.show()
plt.savefig(
    os.path.join(
        plot_dir,
        f'classified_intensity_comparison_{sample_id_img}.{plot_file_type}'
    )
)

#   Next, plot 5 nuclei of each called cell type for visual inspection
rng = default_rng()
examples_per_type = 5

for cell_type in df['cell_type'].unique():
    #   Randomly pick 5 distinct rows for this cell type
    indices = rng.choice(
        df[df['cell_type'] == cell_type].index,
        examples_per_type, replace = False
    )

    for i in indices:
        plt.close('all')
        fig = plot_roi(imgs, props, i)
        fig.savefig(
            os.path.join(
                plot_dir, f'{cell_type}_{i}_{sample_id_img}.{plot_file_type}'
            )
        )


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

#   Print number of cells of each type
for cell_type in df['cell_type'].unique():
    print(f'Total number of {cell_type} cells: {raw[cell_type].sum()}')

#   Correct the data type
for column_name in list(df['cell_type'].unique()) + ['n_cells']:
    raw[column_name] = raw[column_name].astype(int)

#   Print some quick stats about the cell counts per spot
prop_nonzero = round(100 * np.count_nonzero(raw['n_cells']) / raw.shape[0], 1)
print(f"Percentage of spots with at least one cell: {prop_nonzero}%")
print(f"Mean number of cells per spot: {round(np.mean(raw['n_cells']), 3)}")
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
raw.index = raw['barcode'] + '_' + sample_id_img
raw.index.name = 'key'
raw.drop(['row', 'col', 'x', 'y', 'barcode'], axis = 1, inplace = True)

#   Save
raw.to_csv(clusters_path, float_format="%.3f")
