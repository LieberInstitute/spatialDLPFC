#!/usr/bin/env python

###--------------------------------------------LOAD LIBRARIES

import numpy as np
import pandas as pd
import scanpy as sc
import tspex
import seaborn as sns
from matplotlib.pyplot import figure
import matplotlib.pyplot as plt
sns.set(style = "whitegrid", font_scale=1)

###--------------------------------------------LOAD DATA

adata = sc.read_h5ad("../../../processed-data/gene-risk-LR-analysis/04-specificity_LR_targets/sce.h5ad")
adata = adata[adata.obs['cellType_broad_hc']!='drop']
adata.X = adata.X.astype(float)
adata.layers['logcounts'] = adata.layers['logcounts'].astype(float)
adata.obsm['tsne'] = adata.obsm['TSNE'].astype(float)

###--------------------------------------------VISUALISE GENES OF INTEREST

df = pd.read_csv('../../../processed-data/gene-risk-LR-analysis/02-SCZ_LRs/SCZ_top_interactions.csv')
unique = list(np.unique(list(df['genesymbol_intercell_source'])+list(df['genesymbol_intercell_target'])))
sc.pl.tsne(adata, layer = 'logcounts', color = unique)

adata.var_names_make_unique(join='.')

df = pd.DataFrame((adata.X.astype(float)).toarray().T, index = adata.var_names, columns = adata.obs['cellType_layer'])

df_med = df.groupby(by=['cellType_layer'], axis=1).median().loc[unique]

###--------------------------------------------RUN ANALYSIS

sns.set(font_scale=0.6)
tso = tspex.TissueSpecificity(df_med, 'tau', log=False)
sns.heatmap(tso.tissue_specificity.to_frame().sort_values(by=0, ascending = False), annot = True) 
plt.savefig('../../../plots/gene-risk-LR-analysis/04-specificity_LR_targets/heatmap_tau_allnuclei.pdf', dpi = 300)
plt.show()

tso.tissue_specificity.to_frame().to_csv('../../../processed-data/gene-risk-LR-analysis/04-specificity_LR_targets/tau_targets_allnuclei.csv')

sns.set(font_scale=1)
spm = tspex.TissueSpecificity(df_med, 'spm', log = False)
plt.figure(figsize=(30,10))
sns.clustermap(spm.tissue_specificity, cmap='viridis', yticklabels = True, col_cluster = False)
plt.savefig('../../../plots/gene-risk-LR-analysis/04-specificity_LR_targets/clustermap_spm_allnuclei.pdf', dpi = 300)
plt.show()

spm.tissue_specificity.to_csv('../../../processed-data/gene-risk-LR-analysis/04-specificity_LR_targets/spm_targets_allnuclei.csv')

