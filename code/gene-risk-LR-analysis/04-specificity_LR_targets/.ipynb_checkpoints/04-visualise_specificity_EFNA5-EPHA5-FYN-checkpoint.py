#!/usr/bin/env python

###--------------------------------------------LOAD LIBRARIES

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

###--------------------------------------------LOAD AND PLOT TAU STATISTIC

df = pd.read_csv('../../../processed-data/gene-risk-LR-analysis/04-specificity_LR_targets/tau_targets_allnuclei.csv', index_col = 0)
df = df.loc[['EFNA5', 'EPHA5', 'FYN']]

sns.set(font_scale=1.7)
sns.heatmap(df.sort_values(by=['0'], ascending = False), annot = True) 
plt.savefig('../../../plots/gene-risk-LR-analysis/04-specificity_LR_targets/heatmap_tau_allnuclei_FYN-EFNA5-EPHA5.pdf', dpi = 300)
plt.show()

###--------------------------------------------LOAD AND PLOT SPM STATISTIC

df = pd.read_csv('../../../processed-data/gene-risk-LR-analysis/04-specificity_LR_targets/spm_targets_allnuclei.csv', index_col = 0)
df = df.loc[['EFNA5', 'EPHA5', 'FYN']]
plt.figure(figsize=(30,10))
sns.clustermap(df, cmap='viridis', yticklabels = True, col_cluster = False, row_cluster=False)
plt.savefig('../../../plots/gene-risk-LR-analysis/04-specificity_LR_targets/heatmap_tau_allnuclei_FYN-EFNA5-EPHA5.pdf', dpi = 300)
plt.show()