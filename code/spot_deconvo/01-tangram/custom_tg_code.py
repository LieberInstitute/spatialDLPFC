import pandas as pd
import tangram as tg
import logging
import numpy as np

#   tg.count_cell_annotations is limited in that it requires cell centroids,
#   which require quite a bit of work to obtain except via squidpy
#   segmentation, which we've found to be slow and ineffective on H&E histology
#   images. This function is a replacement which doesn't return centroids in
#   adata_sp.obsm['tangram_ct_count'] and therefore only requires cell counts
#   as input, but performs deconvolution exactly the same way.
#
#   This code is based off of Tangram code from version 1.0.2, and more
#   precisely, from:   
#   https://github.com/broadinstitute/Tangram/blob/19c60d17d883b96aa227b6908576f5d1ed49288c/tangram/utils.py
def custom_count_cell_annotations(
    adata_map, adata_sc, adata_sp, annotation="cell_type",
    cell_count_var="cell_count", threshold=0.5,
):
    xs = adata_sp.obsm["spatial"][:, 1]
    ys = adata_sp.obsm["spatial"][:, 0]

    cell_count = adata_sp.obs[cell_count_var]

    # create a dataframe
    df_vox_cells = df_vox_cells = pd.DataFrame(
        data={"x": xs, "y": ys, "cell_n": cell_count},
        index=list(adata_sp.obs.index),
    )

    # get the most probable voxel for each cell
    resulting_voxels = np.argmax(adata_map.X, axis=1)

    # create a list with filtered cells and the voxels where they have been
    # placed with the highest probability a cell i is filtered if
    # F_i > threshold'
    if "F_out" in adata_map.obs.keys():
        filtered_voxels_to_types = [
            (j, adata_sc.obs[annotation][k])
            for i, j, k in zip(
                adata_map.obs["F_out"], resulting_voxels, range(len(adata_sc))
            )
            if i > threshold
        ]
        vox_ct = filtered_voxels_to_types
    else:
        vox_ct = [(resulting_voxels, adata_sc.obs[annotation])]

    df_classes = tg.one_hot_encoding(adata_sc.obs[annotation])
    for index, i in enumerate(df_classes.columns):
        df_vox_cells[i] = 0

    for k, v in vox_ct:
        df_vox_cells.iloc[k, df_vox_cells.columns.get_loc(v)] += 1

    adata_sp.obsm["tangram_ct_count"] = df_vox_cells
    logging.info(
        f"spatial cell count dataframe is saved in `obsm` `tangram_ct_count` of the spatial AnnData."
    )
