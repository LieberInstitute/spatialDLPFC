import os, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import torch
import tangram as tg
import seaborn as sns
import scipy as sp
import getopt
import pyhere
from pathlib import Path

#  This python script will be invoked by an array of shell scripts, where each
#  index corresponds to a sample in the 12 DLPFC spatial samples published in
#  Nature: Neuroscience. The snRNAseq AnnData contains many brain regions
#  ("pan-brain").
#
#  Some plots are generated for the same marker genes mentioned by Kristen,
#  though since this is a "real run", these markers are allowed to be in the
#  set of training genes (to produce the best fit).

plot_dir = pyhere.here('tangram_libd', 'plots', '03_nn_run')
out_dir = pyhere.here('tangram_libd', 'processed-data', '03_nn_run', 'tangram_out')
Path(plot_dir).mkdir(parents=True, exist_ok=True)
Path(out_dir).mkdir(parents=True, exist_ok=True)

sc_path = pyhere.here('tangram_libd', 'processed-data', '03_nn_run', 'sce_pan.h5ad')
sp_path = pyhere.here('tangram_libd', 'processed-data', '03_nn_run', 'visium_DLPFC.h5ad')
marker_path = pyhere.here('tangram_libd', 'processed-data', '03_nn_run', 'pan_markers.txt')
sample_path = pyhere.here('tangram_libd', 'processed-data', '03_nn_run', 'brain_samples.txt')
VisHigh_path = pyhere.here('tangram_libd', 'processed-data', '03_nn_run', 'VisHigh_overlaps.txt')
VisLow_path = pyhere.here('tangram_libd', 'processed-data', '03_nn_run', 'VisLow_overlaps.txt')

#  Genes we want to plot predicted vs. actual expression for, and their
#  corresponding gene symbols
select_genes_names = ['SNAP25', 'MBP', 'PCP4', 'CCK', 'RORB', 'ENC1', 'CARTPT', 'NR4A2', 'RELN']
select_genes = ["ENSG00000132639", "ENSG00000197971", "ENSG00000183036", 
    "ENSG00000187094", "ENSG00000198963", "ENSG00000171617", "ENSG00000164326",
    "ENSG00000153234", "ENSG00000189056"]
    
print('Using tangram version:', tg.__version__)

#  import Visium high and low expressing genes
with open(VisHigh_path, 'r') as f:
    VisHigh = f.read().splitlines()
with open(VisLow_path, 'r') as f:
    VisLow = f.read().splitlines()
    
#  Grab the full list of sample names we will subset from
with open(sample_path, 'r') as f:
    sample_names = f.read().splitlines()

#  Recieve the '-i' argument, an integer in [1, 4] corresponding to the index in
#  the sample_names list
try:
    opts, args = getopt.getopt(sys.argv[1:], "i:", ["index="])
except getopt.GetoptError:
    print('test.py -i <sample index in 1-4>')
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
for i in range(len(select_genes)):
    gene = select_genes[i]
    gene_name = select_genes_names[i]
    
    #  Verify genes of interest were measured in the experiment
    assert gene in ad_sp.var.gene_id
    assert gene in ad_sc.var.gene_id
    
    if gene in markers:
        print('Gene', gene, '(' + gene_name + ') is in the training set.')
    else:
        print('Gene', gene, '(' + gene_name + ') is in the test set.')

#  Subset and otherwise prepare objects for mapping
sc.pp.normalize_total(ad_sc)

ad_sp = ad_sp[ad_sp.obs['sample_id'] == sample_name, :]

tg.pp_adatas(ad_sc, ad_sp, genes=markers)
assert ad_sc.uns['training_genes'] == ad_sp.uns['training_genes']

#  Mapping step using GPU
ad_map = tg.map_cells_to_space(
    adata_sc=ad_sc,
    adata_sp=ad_sp,
    device='cuda:0'
)

ad_map.write_h5ad(os.path.join(out_dir, 'ad_map_' + sample_name + '.h5ad'))

#  Generate plots
tg.plot_cell_annotation(ad_map, ad_sp, annotation='cellType', x='pxl_row_in_fullres', y='pxl_col_in_fullres', nrows=5, ncols=4)
f = plt.gcf()
f.savefig(os.path.join(plot_dir, 'cell_annotation_' + sample_name + '.png'), bbox_inches='tight')

tg.plot_training_scores(ad_map, bins=50, alpha=.5)
f = plt.gcf()
f.savefig(os.path.join(plot_dir, 'train_scores' + sample_name + '.pdf'), bbox_inches='tight')

#  Project all cells based on trained mapping, and save the result
ad_ge = tg.project_genes(adata_map=ad_map, adata_sc=ad_sc)
ad_ge.write_h5ad(os.path.join(plot_dir, 'ad_ge' + sample_name + '.h5ad'))

#  Plot expected vs. actual expression maps for particular test genes of
#  interest
select_genes = [x.lower() for x in select_genes]
tg.plot_genes(select_genes, adata_measured=ad_sp, adata_predicted=ad_ge, x='pxl_row_in_fullres', y='pxl_col_in_fullres')
f = plt.gcf()
f.savefig(os.path.join(plot_dir, 'mapped_select_genes_' + sample_name + '.pdf'), bbox_inches='tight')

VisHigh = [x.lower() for x in VisHigh]
tg.plot_genes(VisHigh, adata_measured=ad_sp, adata_predicted=ad_ge, x='pxl_row_in_fullres', y='pxl_col_in_fullres')
f = plt.gcf()
f.savefig(os.path.join(plot_dir, 'mapped_VisHigh_genes_' + sample_name + '.pdf'), bbox_inches='tight')

VisHigh = [x.lower() for x in VisLow]
tg.plot_genes(VisLow, adata_measured=ad_sp, adata_predicted=ad_ge, x='pxl_row_in_fullres', y='pxl_col_in_fullres')
f = plt.gcf()
f.savefig(os.path.join(plot_dir, 'mapped_VisLow_genes_' + sample_name + '.pdf'), bbox_inches='tight')

#  Compute average cosine similarity for test genes
df_all_genes = tg.compare_spatial_geneexp(ad_ge, ad_sp)
test_score = np.mean(df_all_genes.score[np.logical_not(df_all_genes.is_training)])
print('Average test score:', round(float(test_score), 4))

#  Compute average cosine similarity for training genes
train_score = np.mean(df_all_genes.score[df_all_genes.is_training])
print('Average training score:', round(float(train_score), 4))
