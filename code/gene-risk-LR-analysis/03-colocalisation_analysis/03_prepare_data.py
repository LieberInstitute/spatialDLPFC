#!/usr/bin/env python
###--------------------------------------------LOAD LIBRARIES

import scanpy as sc
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import glob
import anndata

###--------------------------------------------COMBINE IMMUNOFLUORESCENCE DATA TO C2L RESULTS

count = 0
for sample in glob.glob('../../../../DLPFC_Visium_LIBD/processed-data/01-tangram/IF/layer/*/ad_sp_orig.h5ad'):
    print(sample.split('/')[9])
    df = pd.read_csv('../../../../DLPFC_Visium_LIBD/processed-data/spot_deconvo/03-cell2location/03-cell2location/IF/layer/'+sample.split('/')[9]+'/clusters.csv')
    adata = sc.read(sample)
    adata.obs['key'] = adata.obs.index+'_'+adata.obs['sample_id'].astype(str)
    print(adata.obs['sample_id'].astype(str))
    adata.obs['key'].astype('category')
    
    adata.obs.reset_index(inplace = True)
    adata.var.drop(columns = ['n_cells'], inplace = True)
    adata.obs = adata.obs.merge(df, left_on = 'key', right_on = 'key', how = 'left')
    adata.obs.set_index('key', inplace = True)
    adata.obs

    if (count == 0):
        combined = adata
    else:
        combined = combined.concatenate(adata)
        print(combined.var.columns)
        if (combined.var['gene_type-1'].astype(str).equals(combined.var['gene_type-0'].astype(str))==True):
            combined.var['gene_type'] = combined.var['gene_type-1']
            combined.var.drop(columns = ['gene_type-1', 'gene_type-0'], inplace = True)
        if (combined.var['gene_name-1'].astype(str).equals(combined.var['gene_name-0'].astype(str))==True):
            combined.var['gene_name'] = combined.var['gene_name-1']
            combined.var.drop(columns = ['gene_name-1', 'gene_name-0'], inplace = True)
        print(combined.var.columns)
    count+=1
    
combined.write('../../../processed-data/gene-risk-LR-analysis/03-colocalisation_analysis/IF_c2l_anndata_combined.h5ad')

adata = sc.read_h5ad('../../../processed-data/gene-risk-LR-analysis/03-colocalisation_analysis/IF_c2l_anndata_combined.h5ad')
adata.obs['index'] = adata.obs['index'].astype(str)+'_'+adata.obs['sample_id'].astype(str)
adata.obs[['index', 'Astro', 'EndoMural', 'Excit_L2_3', 'Excit_L3', 'Excit_L3_4_5', 'Excit_L4', 'Excit_L5', 'Excit_L5_6', 'Excit_L6', 'Inhib', 'Micro', 'OPC', 'Oligo']].to_csv('../../../processed-data/gene-risk-LR-analysis/03-colocalisation_analysis/IF_c2l_annotation.csv', index = False)

###--------------------------------------------COMBINE NON-IMMUNOFLUORESCENCE DATA TO C2L RESULTS

count = 0
for sample in glob.glob('../../../../DLPFC_Visium_LIBD/processed-data/01-tangram/nonIF/layer/*/ad_sp_orig.h5ad'):
    print(sample)
    df = pd.read_csv('../../../../DLPFC_Visium_LIBD/processed-data/spot_deconvo/03-cell2location/03-cell2location/nonIF/layer/'+sample.split('/')[9]+'/clusters.csv')
    adata = sc.read(sample)
    adata.obs['key'] = adata.obs.index+'_'+adata.obs['sample_id'].astype(str)
    print(adata.obs['sample_id'].astype(str))
    adata.obs['key'].astype('category')
    
    adata.obs.reset_index(inplace = True)
    adata.var.drop(columns = ['n_cells'], inplace = True)
    adata.obs = adata.obs.merge(df, left_on = 'key', right_on = 'key', how = 'left')
    adata.obs.set_index('key', inplace = True)
    adata.obs
    if (count == 0):
        combined = adata
    else:
        combined = combined.concatenate(adata)
        print(combined.var.columns)
        if 'gene_type-1' in combined.var.columns:
            if (combined.var['gene_type-1'].astype(str).equals(combined.var['gene_type-0'].astype(str))==True):
                combined.var['gene_type'] = combined.var['gene_type-1']
                combined.var.drop(columns = ['gene_type-1', 'gene_type-0'], inplace = True)
        if 'gene_name-1' in combined.var.columns:
            if (combined.var['gene_name-1'].astype(str).equals(combined.var['gene_name-0'].astype(str))==True):
                combined.var['gene_name'] = combined.var['gene_name-1']
                combined.var.drop(columns = ['gene_name-1', 'gene_name-0'], inplace = True)
        print(combined.var.columns)
    count+=1
    
combined.write('../../../processed-data/gene-risk-LR-analysis/00-OpenTargets_SCZ_risk_genes/nonIF_c2l_anndata_combined.h5ad')

adata = sc.read_h5ad('../../../processed-data/gene-risk-LR-analysis/03-colocalisation_analysis/nonIF_c2l_anndata_combined.h5ad')
adata.obs['index'] = adata.obs['index'].astype(str)+'_'+adata.obs['sample_id'].astype(str)
adata.obs[['index', 'Astro', 'EndoMural', 'Excit_L2_3', 'Excit_L3', 'Excit_L3_4_5', 'Excit_L4', 'Excit_L5', 'Excit_L5_6', 'Excit_L6', 'Inhib', 'Micro', 'OPC', 'Oligo']].to_csv('../../../processed-data/gene-risk-LR-analysis/03-colocalisation_analysis/nonIF_c2l_annotation.csv', index = False)
