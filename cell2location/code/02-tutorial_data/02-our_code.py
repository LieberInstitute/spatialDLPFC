#   Use the tutorial data but my cell2location code. Compare results from this
#   script with the tutorial code (01-exact_tutorial.py) to figure out if my
#   code was the problem

import sys
import os

import scanpy as sc
import anndata
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

import cell2location
from cell2location.utils.filtering import filter_genes
from cell2location.models import RegressionModel
from cell2location.utils import select_slide
from cell2location.plt import plot_spatial
import scvi

from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42 # enables correct plotting of text
import seaborn as sns

import pyhere
from pathlib import Path
from PIL import Image
import json
import session_info

results_folder = pyhere.here(
    "cell2location", "processed-data", "02-tutorial_data",
    "modified_tutorial_results"
)
plot_dir = pyhere.here(
    'cell2location', 'plots', '02-tutorial_data', 'modified_tutorial_results'
)

cell_type_var = 'Subset'    # in single-cell only
plot_file_type = 'png' # 'pdf' is also supported for higher-quality plots

#   Default is 30 in tutorial, but 5 is recommended as an initial guess for
#   Visium data:
#   https://github.com/BayraktarLab/cell2location/blob/master/docs/images/Note_on_selecting_hyperparameters.pdf
N_CELLS_PER_SPOT = 3

#   For spatial mapping model: tutorial recommends 20 as default but to try 200
detection_alpha = 20

################################################################################
#   Function definitions
################################################################################

def perform_regression(
    mod, adata, adata_name, max_epochs, lr, sample_kwargs, plot_name
):
    # Use all data for training (validation not implemented yet, train_size=1)
    mod.train(
        max_epochs=max_epochs, batch_size=sample_kwargs['batch_size'],
        train_size=1, lr=lr, use_gpu=True
    )

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
    mod.save(f'{results_folder}/{adata_name}', overwrite=True)

    # Save anndata object with results
    adata.write_h5ad(
        os.path.join(results_folder, f'{adata_name}_after.h5ad')
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
#   Actually, we'll use the preprocessing code from the tutorial, not my code
################################################################################

# create paths and names to results folders for reference regression and cell2location models
ref_run_name = f'{results_folder}/reference_signatures'
run_name = f'{results_folder}/cell2location_map'

Path(plot_dir).mkdir(parents=True, exist_ok=True)
Path(ref_run_name).mkdir(parents=True, exist_ok=True)
Path(run_name).mkdir(parents=True, exist_ok=True)

adata_vis = sc.datasets.visium_sge(sample_id="V1_Human_Lymph_Node")
adata_vis.obs['sample'] = list(adata_vis.uns['spatial'].keys())[0]

adata_vis.var['SYMBOL'] = adata_vis.var_names
adata_vis.var.set_index('gene_ids', drop=True, inplace=True)

# find mitochondria-encoded (MT) genes
adata_vis.var['MT_gene'] = [gene.startswith('MT-') for gene in adata_vis.var['SYMBOL']]

# remove MT genes for spatial mapping (keeping their counts in the object)
adata_vis.obsm['MT'] = adata_vis[:, adata_vis.var['MT_gene'].values].X.toarray()
adata_vis = adata_vis[:, ~adata_vis.var['MT_gene'].values]

# Read data
adata_ref = sc.read(
    f'{results_folder}/sc_tutorial.h5ad',
    backup_url='https://cell2location.cog.sanger.ac.uk/paper/integrated_lymphoid_organ_scrna/RegressionNBV4Torch_57covariates_73260cells_10237genes/sc.h5ad'
)

adata_ref.var['SYMBOL'] = adata_ref.var.index
# rename 'GeneID-2' as necessary for your data
adata_ref.var.set_index('GeneID-2', drop=True, inplace=True)

# delete unnecessary raw slot (to be removed in a future version of the tutorial)
del adata_ref.raw

selected = filter_genes(adata_ref, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)

# filter the object
adata_ref = adata_ref[:, selected].copy()

################################################################################
#   Estimate reference cell-type signatures
################################################################################

RegressionModel.setup_anndata(
    adata = adata_ref,
    batch_key='Sample',
    labels_key=cell_type_var,
    categorical_covariate_keys=['Method']
)

# create and train the regression model
mod = RegressionModel(adata_ref)
RegressionModel.view_anndata_setup(mod)

adata_ref, mod = perform_regression(
    mod, adata_ref, 'adata_ref', 250, 0.002,
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
#   https://twitter.com/vitaliikl/status/1526510089514926081
mod = cell2location.models.Cell2location(
    adata_vis, cell_state_df=inf_aver,
    # the expected average cell abundance: tissue-dependent
    # hyper-prior which can be estimated from paired histology:
    N_cells_per_location=N_CELLS_PER_SPOT,
    # hyperparameter controlling normalisation of
    # within-experiment variation in RNA detection:
    detection_alpha=detection_alpha
)

cell2location.models.Cell2location.view_anndata_setup(mod)

adata_vis, mod = perform_regression(
    mod, adata_vis, 'adata_vis', 30000, None,
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

#   Subset to this sample
sample_id = "V1_Human_Lymph_Node"
slide = select_slide(adata_vis, sample_id)

cell_types = adata_ref.obs[cell_type_var].cat.categories

#   Plot 8 cell types separately for this sample
sc.pl.spatial(
    slide, cmap='magma',
    # show first 8 cell types
    color=cell_types[:8],
    ncols=4, size=1.3,
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

fig = plot_spatial(
    adata=slide,
    # labels to show on a plot
    color=cell_types[:6], labels=cell_types[:6],
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
