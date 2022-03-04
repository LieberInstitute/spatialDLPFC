#  Run Tangram on DLPFC subset of Matt's Neuron snRNA-seq samples coupled with
#  the 12 NN Visium samples available from spatialLIBD. This code is roughly
#  based off of the tutorial found here:
#
#  https://github.com/broadinstitute/Tangram/blob/master/tutorial_tangram_with_squidpy.ipynb
#
#  This python script will be invoked by an array of shell scripts. Some plots
#  are generated for the same marker genes mentioned by Kristen, though since
#  this is a "real run", these markers are allowed to be in the set of training
#  genes (to produce the best fit).

import os, sys
import getopt
import pyhere
from pathlib import Path

import scanpy as sc
import squidpy as sq
import numpy as np
import pandas as pd
from anndata import AnnData
import pathlib
import matplotlib.pyplot as plt
import matplotlib as mpl
import skimage
import seaborn as sns
import tangram as tg

plot_dir = pyhere.here('tangram_libd', 'plots', '03_nn_run', 'DLPFC')
out_dir = pyhere.here(
    'tangram_libd', 'processed-data', '03_nn_run', 'tangram_out_DLPFC'
)
Path(plot_dir).mkdir(parents=True, exist_ok=True)
Path(out_dir).mkdir(parents=True, exist_ok=True)

sc_path = pyhere.here(
    'tangram_libd', 'processed-data', '03_nn_run', 'sce_DLPFC.h5ad'
)
sp_path = pyhere.here(
    'tangram_libd', 'processed-data', '03_nn_run', 'visium_DLPFC.h5ad'
)
marker_path = pyhere.here(
    'tangram_libd', 'processed-data', '03_nn_run', 'pan_markers.txt'
)
sample_path = pyhere.here(
    'tangram_libd', 'processed-data', '03_nn_run', 'brain_samples.txt'
)

#  Genes we want to plot predicted vs. actual expression for
select_genes_names = [
    'SNAP25', 'MBP', 'PCP4', 'CCK', 'RORB', 'ENC1', 'CARTPT', 'NR4A2', 'RELN'
]

print('Using tangram version:', tg.__version__)

#  Grab the full list of sample names we will subset from
with open(sample_path, 'r') as f:
    sample_names = f.read().splitlines()

#  Recieve the '-i' argument, an integer in [1, 12] corresponding to the index
#  in the sample_names list
try:
    opts, args = getopt.getopt(sys.argv[1:], "i:", ["index="])
except getopt.GetoptError:
    print('test.py -i <sample index in 1-12>')
    sys.exit(2)

for opt, arg in opts:
    assert opt in ('-i', '--index='), opt
    sample_index = int(arg)

#  Determine this particular sample name
sample_name = sample_names[sample_index - 1]

#  Load AnnDatas and list of marker genes
ad_sp = sc.read_h5ad(sp_path)
ad_sc = sc.read_h5ad(sc_path)

with open(marker_path, 'r') as f:
    markers = f.read().splitlines()
    
#  Note when genes of interest are present in the training set  
select_genes = ad_sp.var.gene_id[ad_sp.var.gene_name.isin(select_genes_names)]
assert len(select_genes_names) == len(select_genes)

for i in range(len(select_genes)):
    gene = select_genes[i]
    gene_name = select_genes_names[i]
    
    #  Verify genes of interest were measured in the experiment
    assert gene in ad_sp.var.gene_id.values
    assert gene in ad_sc.var.gene_id.values
    
    if gene in markers:
        print('Gene', gene, '(' + gene_name + ') is in the training set.')
    else:
        print('Gene', gene, '(' + gene_name + ') is in the test set.')

tg.pp_adatas(ad_sc, ad_sp, genes=markers)

#  Mapping step using GPU
ad_map = tg.map_cells_to_space(ad_sc, ad_sp,
    mode="cells",
    density_prior='rna_count_based',
    num_epochs=500,
    device="cuda:0"
)

tg.project_cell_annotations(ad_map, ad_sp, annotation="cellType")
annotation_list = list(pd.unique(ad_sc.obs['cellType']))

#  Plot spatial expression by cell-type label
tg.plot_cell_annotation_sc(ad_sp, annotation_list, perc=0.02)
f = plt.gcf()
f.savefig(
    os.path.join(plot_dir, 'cell_annotation_' + sample_name + '.png'),
    bbox_inches='tight'
)

#  Plot training scores
tg.plot_training_scores(ad_map, bins=20, alpha=.5)
f = plt.gcf()
f.savefig(
    os.path.join(plot_dir, 'train_scores_' + sample_name + '.pdf'),
    bbox_inches='tight'
)

#  Project all cells based on trained mapping, and save the result
ad_ge = tg.project_genes(ad_map=ad_map, ad_sc=ad_sc)
ad_ge.write_h5ad(os.path.join(out_dir, 'ad_ge_' + sample_name + '.h5ad'))

#  Plot low-scoring training genes
genes = ['rragb', 'trim17', 'eno1b'] #  Will need to adjust this for our data
ad_map.uns['train_genes_df'].loc[genes]
tg.plot_genes_sc(genes, adata_measured=ad_sp, adata_predicted=ad_ge, perc=0.02)
f = plt.gcf()
f.savefig(
    os.path.join(
        plot_dir, 'mapped_low_scoring_train_genes_' + sample_name + '.png'
    ),
    bbox_inches='tight'
)

#  Plot genes not present in spatial data
genes=['loc102633833', 'gm5700', 'gm8292'] # adjust for our data
tg.plot_genes_sc(genes, adata_measured=ad_sp, adata_predicted=ad_ge, perc=0.02)
f = plt.gcf()
f.savefig(
    os.path.join(
        plot_dir, 'mapped_missing_visium_genes_' + sample_name + '.png'
    ),
    bbox_inches='tight'
)

#  Compute average cosine similarity for test genes
df_all_genes = tg.compare_spatial_geneexp(ad_ge, ad_sp, ad_sc)
test_score = np.mean(df_all_genes.score[np.logical_not(df_all_genes.is_training)])
print('Average test score:', round(float(test_score), 4))

#  Compute average cosine similarity for training genes
train_score = np.mean(df_all_genes.score[df_all_genes.is_training])
print('Average training score:', round(float(train_score), 4))

tg.plot_auc(df_all_genes)
f = plt.gcf()
f.savefig(
    os.path.join(plot_dir, 'test_auc_' + sample_name + '.png'),
    bbox_inches='tight'
)
