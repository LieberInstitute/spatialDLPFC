import os
import matplotlib.pyplot as plt
import numpy as np
from numpy.random import default_rng
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
    "plots", "spot_deconvo", "02-cellpose", "mask_dilation_test_{}pix_{}.png"
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

################################################################################
#   Function definitions
################################################################################

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

    expanded_masks = img.copy().astype(np.int32)

    for i in range(num_chunks):
        if verbose:
            print(f'Dilating by {dilation_chunked} pixels...')

        expanded_masks = ndimage.grey_dilation(
            expanded_masks,
            size = (dilation_chunked, dilation_chunked)
        )

        if i < num_chunks - 1:
            if verbose:
                print('Inverting...')

            expanded_masks *= -1

    return expanded_masks.astype(img.dtype)

#   Perform a grey dilation of size 'dilation_radius' on 'img'. 'chunk_size' is
#   ignored.
def normal_dilation(img, dilation_radius, chunk_size):
    return ndimage.grey_dilation(
        img, size = (dilation_radius, dilation_radius)
    ).astype(img.dtype)

################################################################################
#   Analysis
################################################################################

#   Use 'regionprops' to form a list of masks and their bounding boxes
props_orig = regionprops(masks)

#   Relabel masks. This will be later used to have a more consistent estimate
#   as to how many nuclei merged due to dilation
temp = masks.copy()
temp[temp > 0] = 1
temp = label(temp, connectivity = 2)
num_nuclei_orig = np.max(temp)

rng = default_rng()
indices = rng.choice(num_nuclei_orig, size = 100, replace = False)

for dilation_radius in dilation_radii:
    #   We'll compare two mask-expanding methods. 
    for dilation_fun, description in zip(
        (normal_dilation, balanced_dilation),
        ('Normal dilation', 'Balanced dilation')
    ):
        #   Dilate the original masks by the appropriate function
        expanded_masks = dilation_fun(masks, dilation_radius, 6)
        
        #   Plot the comparison and save
        fig = plot_roi(masks, expanded_masks, props_orig, mask_index, 5 + dilation_radius)
        fig.savefig(
            str(plot_path).format(
                dilation_radius, description.split()[0].lower()
            )
        )
        
        #   Estimate how many nuclei are contained in a typical expanded mask (we
        #   only want 1, but will often get more from dilation!)
        a = []
        for i in indices:
            a.append(
                np.unique(masks[(expanded_masks == i) & (masks != 0)]).shape[0]
            )
        
        a = np.array(a)
        print(f'------------- {description} by {dilation_radius} pixels:')
        print(f'Average number of nuclei covered per expanded mask: ~{np.mean(a)}')
        print(f'~{round(100 * np.count_nonzero(a > 1) / a.size, 1)}% masks have at least 2 nuclei.')
        print(f'~{round(100 * np.count_nonzero(a == 0) / a.size, 1)}% masks have no nuclei (dilation "ate" the nucleus)')
        
        #   Check if dilation merges previously distinct masks. Here we relabel the
        #   masks manually, since dilation doesn't combine labels for masks that
        #   merge
        temp = expanded_masks.copy()
        temp[temp > 0] = 1
        temp = label(temp, connectivity = 2)
        num_nuclei_new = np.max(temp)
        
        perc_left = round(100 * num_nuclei_new / num_nuclei_orig, 1)
        print(f'Dilation-induced overlap of expanded masks resulted in ~{perc_left}% of masks kept.')
