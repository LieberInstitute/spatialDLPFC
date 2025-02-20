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
#  index corresponds to a sample named by Kristen. She named 4 spatial samples,
#  each of which we could map the set of 12 snRNA-seq against. We would compare
#  the mapped results against known/ strongly expected distributions of gene
#  expression based on manually-determined layer segmentations, for a handful
#  of well-known marker genes.

Path(pyhere.here('tangram_libd', 'processed-data', '02_four_sample_demo', 'four_sample_demo_out')).mkdir(parents=True, exist_ok=True)
sc_path = pyhere.here('tangram_libd', 'processed-data', '02_four_sample_demo', 'sce_dlpfc.h5ad')
sp_path = pyhere.here('tangram_libd', 'processed-data', '02_four_sample_demo', 'visium_dlpfc.h5ad')
marker_path = pyhere.here('tangram_libd', 'processed-data', '02_four_sample_demo', 'markers.txt')
out_dir = pyhere.here('tangram_libd', 'processed-data', '02_four_sample_demo', '02_four_sample_demo_out')
plot_dir = pyhere.here('tangram_libd', 'plots', '02_four_sample_demo')
test_genes = ['SNAP25', 'MBP', 'PCP4', 'CCK', 'RORB', 'ENC1', 'CARTPT', 'NR4A2', 'RELN']
sample_path = pyhere.here('tangram_libd', 'processed-data', '02_four_sample_demo', 'brain_samples.txt')

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

# sample_index = 4
#  Determine this particular sample name
sample_name = sample_names[sample_index - 1]

#  Load AnnDatas and list of marker genes
ad_sp = sc.read_h5ad(sp_path)
ad_sc = sc.read_h5ad(sc_path)

with open(marker_path, 'r') as f:
    markers = f.read().splitlines()

#  Ensure no test genes are used in the training set    
for gene in test_genes:
    if gene in markers:
        print('Removing test gene', gene, 'from the training set.')
        markers.remove(gene)

#  Subset and otherwise prepare objects for mapping
sc.pp.normalize_total(ad_sc)

ad_sp = ad_sp[ad_sp.obs['sample_name'] == sample_name, :].copy()

ad_sc, ad_sp = tg.pp_adatas(ad_sc, ad_sp, genes=markers)
assert ad_sc.var.index.equals(ad_sp.var.index)

#  Mapping step using GPU
ad_map = tg.map_cells_to_space(
    adata_cells=ad_sc,
    adata_space=ad_sp,
    device='cuda: 0'
)

ad_map.write_h5ad(os.path.join(out_dir, 'ad_map_' + sample_name + '.h5ad'))

#  Reload the original objects and subset the spatial object by sample
ad_sp = sc.read_h5ad(sp_path)
ad_sc = sc.read_h5ad(sc_path)
ad_sp = ad_sp[ad_sp.obs['sample_name'] == sample_name, :].copy()

#  Generate plots
tg.plot_cell_annotation(ad_map, annotation='cell_type', x='imagerow', y='imagecol', nrows=5, ncols=4)
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
tg.plot_genes(test_genes, adata_measured=ad_sp, adata_predicted=ad_ge, x='imagerow', y='imagecol')
f = plt.gcf()
f.savefig(os.path.join(plot_dir, 'mapped_test_genes_' + sample_name + '.pdf'), bbox_inches='tight')

#  Compute average cosine similarity for test genes
df_all_genes = tg.compare_spatial_geneexp(ad_ge, ad_sp)
test_score = np.mean(df_all_genes.score[np.logical_not(df_all_genes.is_training)])
print('Average test score (all genes):', round(float(test_score), 4))

test_scores = df_all_genes.score[df_all_genes.index.isin(test_genes)]
print('Test scores for our select genes:\n', test_scores, sep='')
print('Average:', round(float(np.mean(test_scores)), 4))

#  Compute average cosine similarity for training genes
train_score = np.mean(df_all_genes.score[df_all_genes.is_training])
print('Average training score:', round(float(train_score), 4))
