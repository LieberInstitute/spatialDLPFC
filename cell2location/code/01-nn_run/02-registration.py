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

################################################################################
#   Perform regression
################################################################################

#  Load AnnDatas
print('Loading AnnDatas...')
adata_vis = sc.read_h5ad(sp_path)
adata_ref = sc.read_h5ad(sc_path)

#   Prepare anndata for the regression model. The original code from the tutorial
#   was changed, given the info here: https://github.com/BayraktarLab/cell2location/issues/145#issuecomment-1107410480
RegressionModel.setup_anndata(
    adata = adata_ref,
    batch_key = 'donor',
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
        plot_dir, f'training_history.{plot_file_type}'
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