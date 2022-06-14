import os
import matplotlib.pyplot as plt
import numpy as np
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
out_path = pyhere.here(
    'processed-data', 'spot_deconvo', '02-cellpose', '{}', 'clusters.csv'
)
mask_path = pyhere.here(
    'processed-data', 'spot_deconvo', '02-cellpose', 'masks', '{}_mask.npy'
)
spot_path = pyhere.here(
    'processed-data', '01_spaceranger_IF', '{}', 'outs', 'spatial',
    'tissue_positions_list.csv'
)
plot_dir = pyhere.here("plots", "spot_deconvo", "02-cellpose")

Path(plot_dir).mkdir(parents=True, exist_ok=True)

#-------------------------------------------------------------------------------
#   Dataset-specific variables
#-------------------------------------------------------------------------------

names = {0: "junk", 1: "dapi", 2: "gfap", 3: "neun", 4: "olig2", 5: "tmem119"}
thresholds = {
    "neun": 10,
    "olig2": 10,
    "tmem119": 25,
    "gfap": 100 # need to adjust this!
}
cell_types = {
    "neun": "neuron",
    "olig2": "oligo",
    "tmem119": "micro",
    "gfap": "astro"
}
m_per_px = 0.497e-6
spot_radius = 65e-6
area_threshold = 200

dilation_radius = 3

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
out_path = str(out_path).format(sample_id_img)
mask_path = str(mask_path).format(sample_id_img)
sample_id_spot = sample_ids_spot[int(os.environ['SGE_TASK_ID']) - 1]
spot_path = str(spot_path).format(sample_id_spot)

Path(out_path).parents[0].mkdir(parents=True, exist_ok=True)

#-------------------------------------------------------------------------------
#   Read in spot data
#-------------------------------------------------------------------------------

assert set(thresholds.keys()).issubset(set(names.values()))

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
expanded_masks = ndimage.grey_dilation(
    masks, size = (dilation_radius, dilation_radius)
).astype(masks.dtype)

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
general = regionprops_table(expanded_masks, properties=["centroid", "area"])
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
    os.path.join(plot_dir, 'roi_{}.{}'.format(sample_id_img, plot_file_type)),
    bbox_inches='tight'
)

#   Plot mask spatial distribution vs. spot distribution; there should be
#   quite a bit of overlap
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

plt.clf()
fig, axs = plt.subplots(nrows=4, ncols=1, figsize=(10, 10))
for i, t in enumerate(thresholds):
    axs[i].hist(df[t], bins = 50)
    axs[i].set_title(t)

# plt.show()
plt.savefig(
    os.path.join(
        plot_dir, f'fluor_histograms_{sample_id_img}.{plot_file_type}'
    )
)

#   Check what fraction of rows have an unclear cell type (very low
#   intensity in all channels)
temp = (df['neun'] < noise_cutoff) & \
    (df['olig2'] < noise_cutoff) & \
    (df['tmem119'] < noise_cutoff) & \
    (df['gfap'] < noise_cutoff)

#   However, some nuclei don't have a clear cell type with this cutoff
np.sum(temp) / len(temp) # 0.05242629344958128

#   Let's view an example unclear nucleus. It looks like in this case GFAP and
#   TMEM119 dominate but are still quite low in intensity. This likely shows
#   that we aren't looking at a large enough region around the nucleus
bad_index = temp[temp].index[0]
df.loc[bad_index]
fig = plot_roi(imgs, props, bad_index)
fig.show()

#   Filtering out what we think are cells of a non-target type dramatically
#   changes the mean intensity, which is a good sign (we get a clear signal for
#   what intensity looks like in the target cell type)
weights = {}
for t in thresholds:
    print(f'Mean for {t} above cutoff: {np.mean(df[t][df[t] > noise_cutoff])}')
    print(f'Mean for {t} otherwise: {np.mean(df[t])}')
    weights[t] = 1 / np.mean(df[t][df[t] > noise_cutoff])


#   Weight the intensity of each channel inversely by the average intensity seen
#   in cells above the "noise cutoff". Then call cell types by choosing the
#   maximal weighted intensity across the 4 possible types
temp = df[[t for t in thresholds]].copy()
for t in thresholds:
    temp[t] *= weights[t]

temp['max_channel'] = [temp.columns[np.argmax(temp.loc[i])] for i in temp.index]
temp['cell_type'] = [cell_types[x] for x in temp['max_channel']]

#   Keep the spot identity for each cell
temp['idx'] = df['idx']

