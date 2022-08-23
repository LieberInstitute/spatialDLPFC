import sys
import scanpy as sc
import anndata
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import cell2location
from cell2location.models import RegressionModel
from cell2location.utils.filtering import filter_genes
from cell2location.utils import select_slide
from cell2location.plt import plot_spatial
import scvi
from matplotlib import rcParams
from pathlib import Path
rcParams['pdf.fonttype'] = 42 # enables correct plotting of text for PDFs

import pyhere
import os

results_folder = pyhere.here(
    "cell2location", "processed-data", "02-tutorial_data",
    "exact_tutorial_results"
)
plot_dir = pyhere.here(
    'cell2location', 'plots', '02-tutorial_data', 'exact_tutorial_results'
)

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

# prepare anndata for the regression model
cell2location.models.RegressionModel.setup_anndata(
    adata=adata_ref,
    # 10X reaction / sample / batch
    batch_key='Sample',
    # cell type, covariate used for constructing signatures
    labels_key='Subset',
    # multiplicative technical effects (platform, 3' vs 5', donor effect)
    categorical_covariate_keys=['Method']
)

# create the regression model
mod = RegressionModel(adata_ref)

# view anndata_setup as a sanity check
mod.view_anndata_setup()

mod.train(max_epochs=250, use_gpu=True)

mod.plot_history(20)
f = plt.gcf()
f.savefig(
    os.path.join(plot_dir, 'cell_signature_training_history.png')
)
plt.close(f)

# In this section, we export the estimated cell abundance (summary of the posterior distribution).
adata_ref = mod.export_posterior(
    adata_ref, sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'use_gpu': True}
)

# Save model
mod.save(f"{ref_run_name}", overwrite=True)

# Save anndata object with results
adata_file = f"{ref_run_name}/sc.h5ad"
adata_ref.write(adata_file)

mod.plot_QC()
plt.savefig(os.path.join(plot_dir, 'adata_ref_QC.png'))
plt.close('all')

# export estimated expression in each cluster
if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
    inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
else:
    inf_aver = adata_ref.var[[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()

inf_aver.columns = adata_ref.uns['mod']['factor_names']
inf_aver.iloc[0:5, 0:5]

# find shared genes and subset both anndata and reference signatures
intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
adata_vis = adata_vis[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()

# prepare anndata for cell2location model
cell2location.models.Cell2location.setup_anndata(adata=adata_vis, batch_key="sample")

# create and train the model
mod = cell2location.models.Cell2location(
    adata_vis, cell_state_df=inf_aver,
    # the expected average cell abundance: tissue-dependent
    # hyper-prior which can be estimated from paired histology:
    N_cells_per_location=30,
    # hyperparameter controlling normalisation of
    # within-experiment variation in RNA detection:
    detection_alpha=20
)
mod.view_anndata_setup()

mod.train(max_epochs=30000,
          # train using full data (batch_size=None)
          batch_size=None,
          # use all data points in training because
          # we need to estimate cell abundance at all locations
          train_size=1,
          use_gpu=True)

# plot ELBO loss history during training, removing first 100 epochs from the plot
mod.plot_history(1000)
plt.legend(labels=['full data training'])
f = plt.gcf()
f.savefig(
    os.path.join(plot_dir, 'spatial_mapping_training_history.png')
)
plt.close(f)

# In this section, we export the estimated cell abundance (summary of the posterior distribution).
adata_vis = mod.export_posterior(
    adata_vis, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'use_gpu': True}
)

# Save model
mod.save(f"{run_name}", overwrite=True)

# Save anndata object with results
adata_file = f"{run_name}/sp.h5ad"
adata_vis.write(adata_file)

mod.plot_QC()
plt.savefig(os.path.join(plot_dir, 'adata_vis_QC.png'))
plt.close('all')

fig = mod.plot_spatial_QC_across_batches()
fig.savefig(os.path.join(plot_dir, 'spatial_qc_across_batches.png'))
plt.close(fig)

# add 5% quantile, representing confident cell abundance, 'at least this amount is present',
# to adata.obs with nice names for plotting
adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm['q05_cell_abundance_w_sf']

# select one slide
slide = select_slide(adata_vis, 'V1_Human_Lymph_Node')

# plot in spatial coordinates
sc.pl.spatial(slide, cmap='magma',
              # show first 8 cell types
              color=['B_Cycling', 'B_GC_LZ', 'T_CD4+_TfH_GC', 'FDC',
                     'B_naive', 'T_CD4+_naive', 'B_plasma', 'Endo'],
              ncols=4, size=1.3,
              img_key='hires',
              # limit color scale at 99.2% quantile of cell abundance
              vmin=0, vmax='p99.2'
             )
f = plt.gcf()
f.savefig(os.path.join(plot_dir, 'individual_cell_types.png'))
plt.close(f)

# select up to 6 clusters
clust_labels = ['T_CD4+_naive', 'B_naive', 'FDC']
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'V1_Human_Lymph_Node')

fig = plot_spatial(
    adata=slide,
    # labels to show on a plot
    color=clust_col, labels=clust_labels,
    show_img=True,
    # 'fast' (white background) or 'dark_background'
    style='fast',
    # limit color scale at 99.2% quantile of cell abundance
    max_color_quantile=0.992,
    # size of locations (adjust depending on figure size)
    circle_diameter=6,
    colorbar_position='right'
)
f = plt.gcf()
f.savefig(os.path.join(plot_dir, 'multi_cell_types.png'))
plt.close(f)
