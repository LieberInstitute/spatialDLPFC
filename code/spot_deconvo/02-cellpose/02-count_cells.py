import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import tifffile
from scipy.spatial import KDTree
from skimage.measure import regionprops, regionprops_table
import pyhere
from pathlib import Path

sns.set()

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
}
m_per_px = 0.497e-6
spot_radius = 65e-6
area_threshold = 200

plot_file_type = 'png' # 'pdf' is also supported for higher-quality plots

################################################################################
#   Functions
################################################################################

def setup():
    #   Load multi-channel image and masks from segmenting DAPI channel
    imgs = tifffile.imread(img_path)
    masks = np.load(mask_path)

    #   Quantify the mean image fluorescence intensity at each nucleus
    #   identified by segmenting the DAPI channel. This is done for each
    #   (non-lipofuscin) channel
    its = {
        names[i]: regionprops_table(masks, intensity_image=imgs[i], properties=["intensity_mean"])[
            "intensity_mean"
        ]
        for i in range(2, 6)
    }

    #   Create a table containing the centroids and areas of each mask
    #   (nucleus), and add this info to the intensities table
    general = regionprops_table(masks, properties=["centroid", "area"])
    its["area"] = general["area"]
    its["x"] = general["centroid-0"]
    its["y"] = general["centroid-1"]

    return pd.DataFrame(its), masks, imgs

def plot_roi(idx: int, vmax: int = 128, nrows: int = 3, ncols: int = 2):
    bbox = props[idx]["bbox"]
    roi = props[idx]["image"]
    fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(8, 8))
    axs = axs.flatten()
    axs[0].imshow(np.pad(roi, (pad, pad), mode="constant", constant_values=0), aspect="equal")
    axs[0].set_title("Mask")
    axs[0].grid(False)
    for i in range(1, 6):
        axs[i].imshow(
            target[
                i,
                max(0, bbox[0] - pad) : min(target.shape[1], bbox[2] + pad),
                max(0, bbox[1] - pad) : min(target.shape[2], bbox[3] + pad),
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

df, masks, target = setup()
props = regionprops(masks)

idx = 400
pad = 5

# Plot ROI - sanity check
plot_roi(4)
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

# ### Process ROI properties.
#
# Calculates mask parameters. See available parameters [here](https://scikit-image.org/docs/dev/api/skimage.measure.html#skimage.measure.regionprops).
#
# Builds a $k$-d tree to assign masks to spots.

# Build KD tree for nearest neighbor search.
kd = KDTree(raw[["x", "y"]])

#   For each mask, assign a distance to the nearest spot ('dist') and index of
#   that spot in 'df' ('idx'). Add this info to 'df', now called 'combi'
dist, idx = kd.query(df[["x", "y"]])
dist = pd.DataFrame({"dist": dist, "idx": idx})
combi = pd.concat([df, dist], axis=1)

#   Create boolean columns where the values for each channel name are True
#   if the mean intensity for that channel is sufficiently high
for name, t in thresholds.items():
    combi[f"N_{name}"] = combi[name] > t

#   Plot distribution of cell areas
plt.clf()
sns.histplot(data=df, x="area")
plt.savefig(
    os.path.join(
        plot_dir, 'cell_areas_{}.{}'.format(sample_id_img, plot_file_type)
    ),
    bbox_inches='tight'
)

# Filters out masks that are smaller than a threshold and
# masks whose centroid is farther than the spot radius (aka not inside the
# spot).
px_dist = spot_radius / m_per_px  # meter per px.
filtered = combi[(combi.area > area_threshold) & (combi.dist < px_dist)]

#   For each channel, count how many nuclei per spot have sufficiently high
#   mean intensity. Note that a nucleus can be "counted" for multiple channels,
#   or even for none
summed = filtered[[f"N_{name}" for name in thresholds] + \
    ["idx"]].groupby("idx").sum().astype(int)

# Export
out = pd.concat(
    [
        raw,
        summed,
        filtered[["idx", "dist"]].groupby("idx").count().dist.rename("counts"),
    ],
    axis=1,
)
out.fillna(0, inplace=True)
for name in thresholds:
    out[f"N_{name}"] = out[f"N_{name}"].astype(int)

out.counts = out.counts.astype(int)

#   Make compatible with spatialLIBD 'clusters.csv' format
out.index = out['barcode'] + '_' + sample_id_img
out.index.name = 'key'
out.drop(['row', 'col', 'x', 'y', 'barcode'], axis = 1, inplace = True)

out.to_csv(out_path, float_format="%.3f")

#   Print some quick stats about the cell counts per spot
prop_nonzero = round(100 * np.count_nonzero(out['counts']) / out.shape[0], 1)
print(f"Percentage of spots with at least one cell: {prop_nonzero}%")
print(f"Mean number of cells per spot: {round(np.mean(out['counts']), 3)}")
print(f"Max number of cells per spot: {np.max(out['counts'])}")
