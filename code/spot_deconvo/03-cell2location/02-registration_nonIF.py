import sys
import os

import scanpy as sc
import anndata
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

import cell2location
from cell2location.models import RegressionModel
from cell2location.utils import select_slide
from cell2location.plt import plot_spatial
import scvi

from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42 # enables correct plotting of text
import seaborn as sns

import pyhere
import session_info
from pathlib import Path

################################################################################
#   Variable definitions
################################################################################

processed_dir = pyhere.here(
    "processed-data", "spot_deconvo", "03-cell2location", "nonIF"
)
plot_dir = pyhere.here(
    "plots", "spot_deconvo", "03-cell2location", "nonIF"
)

sp_path = os.path.join(processed_dir, 'adata_vis_orig.h5ad')
sc_path = os.path.join(processed_dir, 'adata_ref_orig.h5ad')

#   TODO: adjust these as appropriate!

cell_type_var = 'cellType_broad_hc'    # in single-cell only
cell_count_var = 'count'               # in spatial only
plot_file_type = 'pdf'

#   For spatial mapping model: tutorial recommends 20 as default but to try 200
detection_alpha = 20

#   Number of epochs used to train each portion of the model
epochs_signature = 400 # default: 250
epochs_spatial = 30000

################################################################################
#   Function definitions
################################################################################

def post_training(mod, adata, adata_name, max_epochs, sample_kwargs, plot_name):
    # plot ELBO loss history during training, removing first 10% of epochs from
    # the plot
    mod.plot_history(int(max_epochs / 10))
    f = plt.gcf()
    f.savefig(
        os.path.join(plot_dir, f'{plot_name}.{plot_file_type}'),
        bbox_inches='tight'
    )
    plt.close(f)

    # In this section, we export the estimated cell abundance (summary of the
    # posterior distribution).
    adata = mod.export_posterior(
        adata, sample_kwargs=sample_kwargs
    )

    # Save model
    mod.save(f'{processed_dir}/{adata_name}', overwrite=True)

    # Save anndata object with results
    adata.write_h5ad(
        os.path.join(processed_dir, f'{adata_name}_after.h5ad')
    )

    # Examine reconstruction accuracy to assess if there are any issues with
    # mapping the plot should be roughly diagonal, strong deviations will signal
    # problems
    mod.plot_QC()
    plt.savefig(
        os.path.join(
            plot_dir, f'{adata_name}_QC.{plot_file_type}'
        ),
        bbox_inches='tight'
    )
    plt.close('all')

    return (adata, mod)

################################################################################
#   Load AnnDatas
################################################################################

#  Load AnnDatas
print('Loading AnnDatas...')
adata_vis = sc.read_h5ad(sp_path)
adata_ref = sc.read_h5ad(sc_path)

################################################################################
#   Estimate reference cell-type signatures
################################################################################

#   Prepare anndata for the regression model. The original code from the tutorial
#   was changed, given the info here: https://github.com/BayraktarLab/cell2location/issues/145#issuecomment-1107410480
RegressionModel.setup_anndata(
    adata = adata_ref,
    batch_key = 'subject',
    labels_key = cell_type_var
)

# create and train the regression model
mod = RegressionModel(adata_ref)
RegressionModel.view_anndata_setup(mod)

mod.train(max_epochs=epochs_signature, use_gpu=True)

adata_ref, mod = post_training(
    mod, adata_ref, 'adata_ref', epochs_signature,
    {'num_samples': 1000, 'batch_size': 2500, 'use_gpu': True},
    'cell_signature_training_history'
)

# export estimated expression in each cluster
if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
    inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[
        f'means_per_cluster_mu_fg_{i}'
        for i in adata_ref.uns['mod']['factor_names']
    ]].copy()
else:
    inf_aver = adata_ref.var[[f'means_per_cluster_mu_fg_{i}'
        for i in adata_ref.uns['mod']['factor_names']
    ]].copy()

inf_aver.columns = adata_ref.uns['mod']['factor_names']
inf_aver.iloc[0:5, 0:5]

# find shared genes and subset both anndata and reference signatures
intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
adata_vis = adata_vis[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()

################################################################################
#   Spatial mapping
################################################################################

# prepare anndata for cell2location model
cell2location.models.Cell2location.setup_anndata(
    adata=adata_vis, batch_key="sample"
)

# create and train the model

#   Note that 'detection_alpha' was changed from its default following:
#   https://twitter.com/vitaliikl/status/1526510089514926081.
#
#   Although we have cell counts at each spot, and the tutorial claims to
#   support providing these, the authors clarify use of per-spot counts is not
#   officially supported. After following the instructions here:
#   https://github.com/BayraktarLab/cell2location/issues/103#issuecomment-993583464
#   I encountered an error, and so we'll just use average counts (a scalar).
#
# N_CELLS_PER_SPOT = np.array(
#     adata_vis.obs[cell_count_var], dtype=np.float32
# ).reshape((adata_vis.shape[0], 1))

mod = cell2location.models.Cell2location(
    adata_vis, cell_state_df=inf_aver,
    # the expected average cell abundance: tissue-dependent
    # hyper-prior which can be estimated from paired histology:
    N_cells_per_location = float(np.mean(adata_vis.obs[cell_count_var])), # N_CELLS_PER_SPOT,
    # hyperparameter controlling normalisation of
    # within-experiment variation in RNA detection:
    detection_alpha=detection_alpha
)

cell2location.models.Cell2location.view_anndata_setup(mod)

mod.train(
    max_epochs=epochs_spatial,
    # train using full data (batch_size=None)
    batch_size=None,
    # use all data points in training because
    # we need to estimate cell abundance at all locations
    train_size=1,
    use_gpu=True
)

adata_vis, mod = post_training(
    mod, adata_vis, 'adata_vis', epochs_spatial,
    {'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'use_gpu': True},
    'spatial_mapping_training_history'
)

fig = mod.plot_spatial_QC_across_batches()
fig.savefig(
    os.path.join(
        plot_dir, f'spatial_qc_across_batches.{plot_file_type}'
    ),
    bbox_inches='tight'
)
plt.close(fig)

# add 5% quantile, representing confident cell abundance, 'at least this amount
# is present', to adata.obs with nice names for plotting
adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm[
    'q05_cell_abundance_w_sf'
]

################################################################################
#   Visualization
################################################################################

#   Loop through each sample and produce plots for each
for sample_id in adata_vis.obs['sample'].cat.categories:
        #   Subset to this sample
    slide = select_slide(adata_vis, sample_id)

    cell_types = adata_ref.obs[cell_type_var].cat.categories

    #   Plot cell types individually for this sample
    with mpl.rc_context({'axes.facecolor':  'black', 'figure.figsize': [4.5, 5]}):
        sc.pl.spatial(
            slide, cmap='magma',
            color=cell_types,
            ncols=int(len(cell_types) / 2),
            size=1.3,
            img_key='hires',
            # limit color scale at 99.2% quantile of cell abundance
            vmin=0, vmax='p99.2'
        )
        f = plt.gcf()
        f.savefig(
            os.path.join(
                plot_dir, f'individual_cell_types_{sample_id}.{plot_file_type}'
            ),
            bbox_inches='tight'
        )
        plt.close(f)

    #   Plot all cell types on the same spatial grid for this sample
    with mpl.rc_context({'figure.figsize': (15, 15)}):
        fig = plot_spatial(
            adata=slide,
            # labels to show on a plot
            color=cell_types, labels=cell_types,
            show_img=True,
            # 'fast' (white background) or 'dark_background'
            style='fast',
            # limit color scale at 99.2% quantile of cell abundance
            max_color_quantile=0.992,
            # size of locations (adjust depending on figure size)
            circle_diameter=6,
            colorbar_position='right'
        )
        fig.savefig(
            os.path.join(
                plot_dir, f'multi_cell_types_{sample_id}.{plot_file_type}'
            ),
            bbox_inches='tight'
        )
        plt.close(fig)

################################################################################
#   Export deconvo results
################################################################################

#   For export, use mean cell abundances (not 5% quantile)
adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm[
    'means_cell_abundance_w_sf'
]

clusters = adata_vis.obs[
    ['sample'] + list(adata_vis.uns['mod']['factor_names'])
]
clusters.index.name = 'key'

clusters_by_id = clusters.groupby('sample')

for sample_id in adata_vis.obs['sample'].cat.categories:
    clusters_subset = clusters_by_id.get_group(sample_id)

    Path(os.path.join(processed_dir, sample_id)).mkdir(
        parents=True, exist_ok=True
    )
    clusters_subset.to_csv(os.path.join(processed_dir, sample_id, 'clusters.csv'))

session_info.show(html=False)
