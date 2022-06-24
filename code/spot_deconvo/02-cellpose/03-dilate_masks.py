import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tifffile
from scipy.spatial import KDTree
from scipy import ndimage
from skimage.measure import regionprops, regionprops_table, label
import pyhere

#   This script is a proof of concept, showing how image dilation can be used
#   to enlarge each mask identified by cellpose. The enlarged masks might be
#   better able to overlap fluorescence from the GFAP channel, since this
#   marker doesn't always nicely coincide with nuclei of cells. The goal is to
#   eventually be able to leverage info from the GFAP channel to quantify all
#   4 cell types measured in the IF images (from GFAP and 3 other channels).

sample_info_path = pyhere.here(
    'raw-data', 'sample_info', 'Visium_IF_DLPFC_MasterExcel_01262022.xlsx'
)
mask_path = pyhere.here(
    'processed-data', 'spot_deconvo', '02-cellpose', 'masks', '{}_mask.npy'
)
img_path = pyhere.here('raw-data', 'Images', 'VisiumIF', 'VistoSeg', '{}.tif')
plot_path = pyhere.here(
    "plots", "spot_deconvo", "02-cellpose", "mask_dilation_test_{}pix.png"
)
dilation_radii = range(3, 45, 6)
mask_index = 400
sample_index = 0

#   Determine paths to mask and fluorescence image for the first sample
sample_info = pd.read_excel(sample_info_path, header = 1)[:4]
sample_ids_img = sample_info['Slide SN #'] + '_' + sample_info['Array #']
sample_id_img = sample_ids_img[sample_index]
mask_path = str(mask_path).format(sample_id_img)
img_path = str(img_path).format(sample_id_img)

masks = np.load(mask_path)
imgs = tifffile.imread(img_path)

#   Function to plot what the original ('img') and expanded ('img2') masks
#   look like at a particular mask index ('idx')
def plot_roi(img, img2, props, idx: int, pad: int = 5):
    bbox = props[idx]["bbox"]
    roi = props[idx]["image"]
    fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(8, 8))
    axs = axs.flatten()
    axs[0].imshow(
        np.pad(roi, (pad, pad), mode="constant", constant_values=0),
        aspect="equal"
    )
    axs[0].set_title("Original mask")
    axs[0].grid(False)
    axs[1].imshow(
        img2[
            max(0, bbox[0] - pad) : min(img2.shape[0], bbox[2] + pad),
            max(0, bbox[1] - pad) : min(img2.shape[1], bbox[3] + pad),
        ]
    )
    axs[1].set_title('Expanded mask')
    axs[1].grid(False)
    fig.tight_layout()
    return fig

#   Use 'regionprops' to form a list of masks and their bounding boxes
props_orig = regionprops(masks)

#   Relabel masks. This will be later used to have a more consistent estimate
#   as to how many nuclei merged due to dilation
temp = masks.copy()
temp[temp > 0] = 1
temp = label(temp, connectivity = 2)
num_nuclei_orig = np.max(temp)

for dilation_radius in dilation_radii:
    #   Dilate the original masks
    expanded_masks = ndimage.grey_dilation(
        masks, size = (dilation_radius, dilation_radius)
    ).astype(masks.dtype)

    #   Plot the comparison and save
    fig = plot_roi(masks, expanded_masks, props_orig, mask_index, 5 + dilation_radius)
    fig.savefig(str(plot_path).format(dilation_radius))

    #   Check if dilation merges previously distinct masks. Here we relabel the
    #   masks manually, since dilation doesn't combine labels for masks that
    #   merge
    temp = expanded_masks.copy()
    temp[temp > 0] = 1
    temp = label(temp, connectivity = 2)
    num_nuclei_new = np.max(temp)

    perc_left = round(100 * num_nuclei_new / num_nuclei_orig, 1)
    print(f'Dilation by {dilation_radius} pixels results in ~{perc_left}% of masks kept.')
