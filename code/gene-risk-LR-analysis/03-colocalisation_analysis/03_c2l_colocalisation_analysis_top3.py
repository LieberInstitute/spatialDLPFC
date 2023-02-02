#!/usr/bin/env python

###--------------------------------------------LOAD LIBRARIES
import scanpy as sc
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import copy
import sys
sns.set(rc={'figure.figsize':(10,9)}, font_scale=2.2)


###--------------------------------------------LOAD DATA
adata1 = sc.read("../../../processed-data/gene-risk-LR-analysis/03-colocalisation_analysis/nonIF_c2l_anndata_combined.h5ad")
adata2 = sc.read("../../../processed-data/gene-risk-LR-analysis/03-colocalisation_analysis/IF_c2l_anndata_combined.h5ad")

adata = adata1.concatenate(adata2, join = 'inner')
if (adata.var['gene_name-1'].astype(str).equals(adata.var['gene_name-0'].astype(str))==True):
    adata.var['gene_name'] = adata.var['gene_name-1']
    adata.var.drop(columns = ['gene_name-1', 'gene_name-0'], inplace = True)
    

###--------------------------------------------FUNCTIONS
def coloc_LR(original, ligand, receptor):
    adata = copy.copy(original)
    gene1 = adata.var['gene_name'].str.contains(ligand)
    gene2 = adata.var['gene_name'].str.contains(receptor)

    if (gene1.sum()<1) or (gene2.sum()<1): 
        print('Gene is not present')
        return
    df1 = adata.to_df()[adata.var[gene1].index].rename(columns={adata.var[gene1].index[0] : "Ligand"})
    df2 = adata.to_df()[adata.var[gene2].index].rename(columns={adata.var[gene2].index[0] : "Receptor"})

    adata.obs = adata.obs.merge(df1, left_index = True, right_index = True, how = 'left')
    adata.obs = adata.obs.merge(df2, left_index = True, right_index = True, how = 'left')

    adata.obs['LR'] = (adata.obs['Ligand']>1) * (adata.obs['Receptor']>1)
    adata.obs['control'] = ((adata.obs['Ligand']==0) & (adata.obs['Receptor']==0))
    LR = adata[adata.obs['LR']==True]
    ctrl = adata[adata.obs['control']==True]
    print('adata', np.shape(adata.obs))
    print('LR', np.shape(LR.obs))
    print('universe', np.shape(ctrl.obs))
    return LR, ctrl


def network_top3(adata, tgt):
    probs = adata.obs[['Astro', 'EndoMural', 'Micro', 'Oligo', 'OPC', 'Excit_L2_3', 'Excit_L3',
       'Excit_L3_4_5', 'Excit_L4', 'Excit_L5', 'Excit_L5_6', 'Excit_L6',
       'Inhib']]

    print("Building adjacency matrix for top 3 cells per spot...")
    hot1 = copy.copy(probs)
    hot1.iloc[:,:] = 0
    for m in range (0,np.shape(probs)[0]):
        hot1.iloc[m,:][probs.T[probs.index[m]].nlargest(3).index]=1

    adj = hot1.T.dot(hot1)
    np.fill_diagonal(adj.values, 0)
    adj = adj = pd.DataFrame(np.tril(adj), index = adj.index, columns = adj.columns)
    adjmat_notnorm = adj
    adjmat = (adj/adj.sum().sum())
    print("Done.")
    print(adjmat.sum().sum())
    
    sns.heatmap(pd.DataFrame(adjmat), cmap = 'inferno_r')
    plt.savefig("../../../plots/gene-risk-LR-analysis/03-colocalisation_analysis/%s_c2l_3cells_heatmap.pdf" % (tgt), dpi = 150, bbox_inches = 'tight')
    plt.show()
    plt.clf()
    adjmat_notnorm.to_csv("../../../processed-data/gene-risk-LR-analysis/03-colocalisation_analysis/%s_c2l_3cells_adjacencymatrix_notnormalised.csv" % (tgt))
    adjmat.to_csv("../../../processed-data/gene-risk-LR-analysis/03-colocalisation_analysis/%s_c2l_3cells_adjacencymatrix.csv" % (tgt))
    return adjmat


###--------------------------------------------RUN ANALYSIS FOR TOP 3 CELLS
print('Identifying spots with co-localising ', sys.argv[1], ' and ', sys.argv[2])
LR, ctrl = coloc_LR(adata, sys.argv[1], sys.argv[2])
if LR is None:
    print('No spots')
    
if (np.shape(LR.obs)[0]<5):
    print('Less than 5 co expressing spots')
    
print('Calculating LR neighbourhood...')
LRadjmat = network_top3(LR, 'LR_'+sys.argv[1]+'-'+sys.argv[2])
print('Calculating control neighbourhood...')
ctrladjmat = network_top3(ctrl, 'Ctrl_'+sys.argv[1]+'-'+sys.argv[2])

print('Taking ratio...')
sns.heatmap(pd.DataFrame(np.tril(LRadjmat/ctrladjmat), index = LRadjmat.index, columns = LRadjmat.columns), cmap = 'inferno_r')
plt.xticks(rotation=90)
plt.yticks(rotation=0)
plt.savefig("../../../plots/gene-risk-LR-analysis/03-colocalisation_analysis/ratio_%s_c2l_3cells_network.pdf" % ('Ctrl'+sys.argv[1]+'-'+sys.argv[2]), dpi = 150, bbox_inches = 'tight')
plt.show()
plt.clf()

lut = dict(zip(LRadjmat.columns, ["#3BB273","#FF56AF","#663894","#E07000","#D2B037","#D0D1E6","#A6BDDB","#74A9CF","#3690C0","#0570B0","#045A8D","#023858","#E94F37"]))
row_colors = LRadjmat.columns.map(lut)

sns.clustermap(LRadjmat/ctrladjmat, cmap = 'inferno_r', row_colors=row_colors, col_colors=row_colors, col_cluster=False, row_cluster=False)

plt.savefig("../../../plots/gene-risk-LR-analysis/03-colocalisation_analysis/ratio_%s_c2l_3cells_network_celllabels.pdf" % ('Ctrl'+sys.argv[1]+'-'+sys.argv[2]), dpi = 150, bbox_inches = 'tight')
plt.show()
plt.clf()

print('Done')
