###--------------------------------------------LOAD LIBRARIES

import numpy as np
import pandas as pd
import scanpy as sc
import tspex
import seaborn as sns
from matplotlib.pyplot import figure
import matplotlib.pyplot as plt
sns.set(style = "whitegrid", font_scale=1.5)

###--------------------------------------------LOAD DATA

adata = sc.read_h5ad("../../../processed-data/gene-risk-LR-analysis/04-specificity_LR_targets/sce.h5ad")
adata = adata[adata.obs['cellType_broad_hc']!='drop']
adata.X = adata.X.astype(float)
adata.layers['logcounts'] = adata.layers['logcounts'].astype(float)
adata.obsm['tsne'] = adata.obsm['TSNE'].astype(float)

###--------------------------------------------RUN ANALYSIS
phos = pd.DataFrame(
    {'genesymbol_intercell_source': ['FYN', 'EFNA5', 'FYN'],
     'genesymbol_intercell_target': ['EFNA5', 'EPHA5', 'EPHA5']
    })

for n in range(0, np.shape(phos)[0]):
    print('###-------------------------------------------------------------------', phos.iloc[n,0], ' - ', phos.iloc[n,1], '-------------------------------------------------------------------###')
    df1 = adata[:, adata.var['gene_name'].str.contains(phos.iloc[n,0])].to_df().rename(columns={phos.iloc[n,0] : "Ligand"})
    df2 = adata[:, adata.var['gene_name'].str.contains(phos.iloc[n,1])].to_df().rename(columns={phos.iloc[n,1] : "Receptor"})
    df = df1.merge(df2, left_index = True, right_index = True, how = 'left')
    df['LR'] = (df['Ligand']>1) * (df['Receptor']>1)
    
    print('Total number of nuclei ', np.shape(df)[0])
    
    print('% of nuclei expressing ligand ', phos.iloc[n,0], ':', (df['Ligand']>1).sum()/np.shape(df)[0]*100)
    
    print('% of nuclei expressing receptor ', phos.iloc[n,1], ':', (df['Receptor']>1).sum()/np.shape(df)[0]*100)
    
    print('% of nuclei expressing both ligand and receptor: ', (df['LR']==True).sum()/np.shape(df)[0]*100)

    df['Ligand'].hist(bins = 50)
    plt.show()
    
    df['Receptor'].hist(bins = 50)
    plt.show()
    df = df.merge(adata.obs['cellType_layer'], left_index = True, right_index = True, how = 'left')
    print("% of co-expressing nuclei per layer:")
    print(100*(df[df['LR']==True]['cellType_layer'].value_counts(sort=False)/df['cellType_layer'].value_counts(sort=False)).sort_values(ascending=False))
    if (n==0):
        summary = 100*(df[df['LR']==True]['cellType_layer'].value_counts(sort=False)/df['cellType_layer'].value_counts(sort=False))
        summary = summary.to_frame().rename(columns={'cellType_layer' : phos.iloc[n,0]+' - '+phos.iloc[n,1]})
    else:
        summary = summary.merge(100*(df[df['LR']==True]['cellType_layer'].value_counts(sort=False)/df['cellType_layer'].value_counts(sort=False)).to_frame().rename(columns={'cellType_layer' : phos.iloc[n,0]+' - '+phos.iloc[n,1]}), left_index = True, right_index = True, how = 'left')
        
sns.clustermap(summary, cmap = 'viridis', row_cluster = False, col_cluster = False)
plt.savefig('../../../plots/gene-risk-LR-analysis/05-intracellular_coexpression/smalltgtlist_snRNAseq_coexpression_heatmap.pdf', dpi = 300, bbox_inches = 'tight')
plt.show()


adata.var_names_make_unique(join='.')
markers = ['EFNA5', 'EPHA5', 'FYN']

sns.set(style = "whitegrid", font_scale=0.9)
fig = sc.pl.dotplot(adata, markers, groupby=['cellType_layer'], standard_scale='var', cmap = 'plasma_r', layer = 'logcounts', show = False)
plt.savefig('../../../plots/gene-risk-LR-analysis/05-intracellular_coexpression/dotplot_FYN-EFNA5-EPHA5_expression.pdf', bbox_inches = 'tight')
plt.show()