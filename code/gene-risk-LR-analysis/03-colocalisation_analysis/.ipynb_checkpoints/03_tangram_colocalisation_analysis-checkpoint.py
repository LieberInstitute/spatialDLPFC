#!/usr/bin/env python

###--------------------------------------------LOAD LIBRARIES
import scanpy as sc
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import copy
sns.set(rc={'figure.figsize':(10,9)}, font_scale=2.2)

###--------------------------------------------LOAD DATA
# --Load top interactions of interest
L = ['EFNA5', 'EFNA5']
R = ['EPHA5', 'FYN']

# --Load and combine IF and nonIF data with tangram deconvolution
adata1 = sc.read("../../../processed-data/spot_deconvo/05-shared_utilities/IF/spe.h5ad")
adata1.obs['key'] = adata1.obs.index+'_'+adata1.obs['sample_id'].astype(str)

adata2 = sc.read("../../../processed-data/spot_deconvo/05-shared_utilities/nonIF/spe.h5ad")
adata2.obs['key'] = adata2.obs.index+'_'+adata2.obs['sample_id'].astype(str)

adata = adata1.concatenate(adata2)

df1 = pd.read_csv("../../../processed-data/01-tangram/IF/layer/raw_results/clusters.csv")
df2 = pd.read_csv("../../../processed-data/01-tangram/nonIF/layer/raw_results/clusters.csv")

df = pd.concat([df1, df2])

adata.obs.reset_index(inplace = True)
adata.obs = adata.obs.merge(df, left_on = 'key', right_on = 'key', how = 'left')
adata.obs.set_index('key', inplace = True)

print(np.shape(adata))
if (adata.var['gene_name-1'].astype(str).equals(adata.var['gene_name-0'].astype(str))==True):
    adata.var['gene_name'] = adata.var['gene_name-1']
    adata.var.drop(columns = ['gene_name-1', 'gene_name-0'], inplace = True)
    
###--------------------------------------------FUNCTIONS

# --Identify co-localising spots (raw counts > 1 for genes of interest)
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
    LR = adata[adata.obs['LR']==True]
    ctrl = adata[adata.obs['LR']==False]
    return LR, ctrl

# --Construct adjacency matrix, plot co-localisation network and heatmap
def network_tangram(original, tgt):
    print('###----------------------------------------------------', tgt,'----------------------------------------------------##')
    adata = copy.copy(original)
    probs = adata.obs[['Astro', 'EndoMural', 'Micro', 'Oligo', 'OPC',
           'Excit_L2_3', 'Excit_L3', 'Excit_L3_4_5', 'Excit_L4', 'Excit_L5',
           'Excit_L5_6', 'Excit_L6', 'Inhib']].astype(int)
    probs = probs.fillna(0)
    adj = probs.T.dot(probs)
    np.fill_diagonal(adj.values, 0)
    adj = pd.DataFrame(np.tril(adj), index = adj.index, columns = adj.columns)
    adjmat_notnorm = adj
    adjmat = (adj/adj.sum().sum())
    print("Done.")
    print(adjmat.sum().sum())
    
    sns.heatmap(pd.DataFrame(np.tril(adjmat), index = adjmat.index, columns = adjmat.columns), cmap = 'inferno_r')
    plt.title(tgt)
    plt.savefig("../../../plots/gene-risk-LR-analysis/03-colocalisation_analysis/%s_tangram_heatmap.pdf" % (tgt), dpi = 150, bbox_inches = 'tight')
    plt.show()
    plt.clf()
    
    adjmat_notnorm.to_csv("../../../plots/gene-risk-LR-analysis/03-colocalisation_analysis/%s_tangram_adjacencymatrix_notnormalised.csv" % (tgt))
    adjmat.to_csv("../../../processed-data/gene-risk-LR-analysis/03-colocalisation_analysis/%s_tangram_adjacencymatrix.csv" % (tgt))
    
    return adjmat

###--------------------------------------------RUN ANALYSIS FOR TANGRAM DECONVOLUTED DATA
for n in range(0,len(L)):
    print(n)
    LR, ctrl = coloc_LR(adata, L[n], R[n])
    if LR is None:
        print('No spots')
        continue
    if (np.shape(LR.obs)[0]<5):
        print('Less than 5 co expressing spots')
        continue
    LRadjmat = network_tangram(LR, 'LR_'+L[n]+'-'+R[n])
    ctrladjmat = network_tangram(ctrl, 'Ctrl_'+L[n]+'-'+R[n])
    
    sns.heatmap(pd.DataFrame(np.tril(LRadjmat.replace(0.0, 0.000001)/ctrladjmat.replace(0.0, 0.000001)), index = LRadjmat.index, columns = LRadjmat.columns), cmap = 'inferno_r')
    plt.title(L[n]+'->'+R[n])
    plt.savefig("../../../plots/gene-risk-LR-analysis/03-colocalisation_analysis/ratio_%s_tangram_network.pdf" % ('Ctrl'+L[n]+'-'+R[n]), dpi = 150, bbox_inches = 'tight')
    plt.show()
    plt.clf()
    
    lut = dict(zip(LRadjmat.columns,["#3BB273","#FF56AF","#663894","#E07000","#D2B037","#D0D1E6","#A6BDDB","#74A9CF","#3690C0","#0570B0","#045A8D","#023858","#E94F37"]))
    row_colors = LRadjmat.columns.map(lut)

    sns.clustermap(LRadjmat/ctrladjmat, cmap = 'inferno_r', row_colors=row_colors, col_colors=row_colors, col_cluster=False, row_cluster=False)

    plt.savefig("../../../plots/gene-risk-LR-analysis/03-colocalisation_analysis/ratio_%s_tangram_network_celllabels.pdf" % ('Ctrl'+L[n]+'-'+R[n]), dpi = 150, bbox_inches = 'tight')
    plt.show()
    plt.clf()