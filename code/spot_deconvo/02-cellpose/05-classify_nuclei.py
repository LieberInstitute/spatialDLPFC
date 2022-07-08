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

clusters_path = pyhere.here(
    'processed-data', 'spot_deconvo', '02-cellpose', '{}', 'clusters.csv'
)
cells_path = pyhere.here(
    'processed-data', 'spot_deconvo', '02-cellpose', '{}', 'cell_metrics.csv'
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
m_per_px = 0.497e-6
spot_radius = 65e-6

#   Mean fluorescence below this value will be considered an indication of an 
#   absence of the given cell type
noise_cutoff = 3

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
sample_info = pd.read_excel(sample_info_path, header = 1)[:4]
sample_ids_img = sample_info['Slide SN #'] + '_' + sample_info['Array #']
sample_ids_spot = 'Br' + sample_info['BrNumbr'].astype(int).astype(str) + \
    '_' + pd.Series([x.split('_')[1] for x in sample_info['Br_Region']]) + \
    '_IF'

sample_id_img = sample_ids_img[int(os.environ['SGE_TASK_ID']) - 1]
img_path = str(img_path).format(sample_id_img)
clusters_path = str(clusters_path).format(sample_id_img)
cells_path = str(cells_path).format(sample_id_img)
expanded_mask_path = str(expanded_mask_path).format(sample_id_img)
df_path = str(df_path).format(sample_id_img)

sample_id_spot = sample_ids_spot[int(os.environ['SGE_TASK_ID']) - 1]
spot_path = str(spot_path).format(sample_id_spot)

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
expanded_masks = np.load(expanded_masks_path, expanded_masks)

props = regionprops(expanded_masks)

#-------------------------------------------------------------------------------
#   Exploratory plot: show 5 random nuclei
#-------------------------------------------------------------------------------

indices = rng.choice(df.index, 5, replace = False)

fig = plot_roi(imgs, props, indices)
fig.savefig(
    os.path.join(plot_dir, 'random_rois_{}.{}'.format(sample_id_img, plot_file_type)),
    bbox_inches='tight'
)

#-------------------------------------------------------------------------------
#   Explore how fluorescence of each channel is distributed
#-------------------------------------------------------------------------------

#   Histograms of raw fluorescence values in each channel
plt.clf()
fig, axs = plt.subplots(nrows=len(cell_types), ncols=1, figsize=(10, 10))
for i, channel in enumerate(cell_types.keys()):
    axs[i].hist(df[channel], bins = 100)
    axs[i].set_title(channel)

plt.savefig(
    os.path.join(
        plot_dir, f'fluor_histograms_raw_{sample_id_img}.{plot_file_type}'
    )
)

#   Histograms of trimmed fluorescence values in each channel (above a
#   "noise cutoff")
plt.clf()
fig, axs = plt.subplots(nrows=len(cell_types), ncols=1, figsize=(10, 10))
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
    (df['tmem119'] < noise_cutoff)

#   However, some nuclei don't have a clear cell type with this cutoff
bad_perc = round(100 * np.sum(temp) / len(temp), 1) # ~15%
print(f"{bad_perc}% of nuclei don't have a clear cell type.")

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
    (df['tmem119'] < noise_cutoff),
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
fig, axs = plt.subplots(nrows=len(cell_types), ncols=1, figsize=(10, 10))
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
examples_per_type = 5

plt.close('all')
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
