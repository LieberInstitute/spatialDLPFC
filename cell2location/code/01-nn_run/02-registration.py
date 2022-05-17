import sys
import os

import scanpy as sc
import anndata
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.backends.backend_pdf import PdfPages

import cell2location
from cell2location.models import RegressionModel
import scvi

from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42 # enables correct plotting of text
import seaborn as sns

import pyhere
import session_info

################################################################################
#   Variable definitions
################################################################################

processed_dir = pyhere.here('cell2location', 'processed-data', '01-nn_run')
plot_dir = pyhere.here('cell2location', 'plots', '01-nn_run')

sp_path = os.path.join(processed_dir, 'adata_vis_orig.h5ad')
sc_path = os.path.join(processed_dir, 'adata_ref.h5ad')


# create paths and names to results folders for reference regression and
# cell2location models
ref_run_name = f'{processed_dir}/reference_signatures'
run_name = f'{processed_dir}/cell2location_map'

#   Naming conventions used for different columns in the spatial AnnData
cell_type_var = 'cellType'    # in single-cell only

plot_file_type = 'png' # 'pdf' is also supported for higher-quality plots

#   Default is 30 in tutorial, but 5 is recommended as an initial guess for
#   Visium data:
#   https://github.com/BayraktarLab/cell2location/blob/master/docs/images/Note_on_selecting_hyperparameters.pdf
N_CELLS_PER_SPOT = 5 

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
    batch_key = 'donor', # tried 'processBatch' as well with poor results
    labels_key = cell_type_var
)

# create and train the regression model
mod = RegressionModel(adata_ref)
RegressionModel.view_anndata_setup(mod)

# Use all data for training (validation not implemented yet, train_size=1)
mod.train(max_epochs=250, batch_size=2500, train_size=1, lr=0.002, use_gpu=True)

# plot ELBO loss history during training, removing first 20 epochs from the plot
mod.plot_history(20)
f = plt.gcf()
f.savefig(
    os.path.join(
        plot_dir, f'cell_signature_training_history.{plot_file_type}'
    ),
    bbox_inches='tight'
)

# In this section, we export the estimated cell abundance (summary of the
# posterior distribution).
adata_ref = mod.export_posterior(
    adata_ref, sample_kwargs={
        'num_samples': 1000, 'batch_size': 2500, 'use_gpu': True
    }
)

# Save model
mod.save(f"{ref_run_name}", overwrite=True)

# Save anndata object with results
adata_ref.write_h5ad(
    os.path.join(processed_dir, 'adata_ref_after.h5ad')
)

mod.plot_QC()

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
    # within-experiment variation in RNA detection (using default here):
    detection_alpha=20 # default: 200
)

cell2location.models.Cell2location.view_anndata_setup(mod)

mod.train(
    max_epochs=30000,
    batch_size=None, # train using full data (batch_size=None)
    # use all data points in training because
    # we need to estimate cell abundance at all locations
    train_size=1,
    use_gpu=True
)

# plot ELBO loss history during training, removing first 100 epochs from the plot
mod.plot_history(1000)
plt.legend(labels=['full data training']);
f = plt.gcf()
f.savefig(
    os.path.join(
        plot_dir, f'spatial_mapping_training_history.{plot_file_type}'
    ),
    bbox_inches='tight'
)

# In this section, we export the estimated cell abundance (summary of the
# posterior distribution).
adata_vis = mod.export_posterior(
    adata_vis,
    sample_kwargs={
        'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'use_gpu': True
    }
)

# Save model
mod.save(f"{run_name}", overwrite=True)

# Save anndata object with results
adata_vis.write_h5ad(os.path.join(processed_dir, 'adata_vis_after.h5ad'))

# Examine reconstruction accuracy to assess if there are any issues with
# mapping the plot should be roughly diagonal, strong deviations will signal
# problems
mod.plot_QC()

fig = mod.plot_spatial_QC_across_batches()
f.savefig(
    os.path.join(
        plot_dir, f'spatial_qc_across_batches.{plot_file_type}'
    ),
    bbox_inches='tight'
)

session_info.show(html=False)
