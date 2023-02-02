#!/usr/bin/env python

###--------------------------------------------LOAD LIBRARIES
import numpy as np
import pandas as pd

###--------------------------------------------LOAD DATA
#Gene list from OpenTargets platform, downloaded on 12-08-2022
# Filtered to only includes genes with positive genetic association

dis = pd.read_csv('../../../raw-data/gene-risk-LR-analysis/Schizophrenia_MONDO_0005090-associated-diseases.tsv', delimiter = '\t')    

disease = 'Schizophrenia'

###--------------------------------------------FORMAT AND FILTER FOR GENETIC ASSOCIATION > 0.1
dis = dis[['symbol', 'objectObject']]
dis.rename(columns={'objectObject': 'genetic_association'}, inplace = True)
dis = dis[dis['genetic_association']>0.1]
if np.shape(dis)[0]>10:
    genetic_association = dis.set_index('symbol')
    genetic_association.to_csv('../../../processed-data/gene-risk-LR-analysis/00-OpenTargets_SCZ_risk_genes/SCZ_risk_genes_01thr.csv')
    
    dis['score']=int(1)
    dis = dis.drop(columns=['genetic_association']).set_index(['symbol'])
    #creating instance of one-hot-encoder
    encoder = OneHotEncoder(handle_unknown='ignore')
    encoder_df = pd.DataFrame(encoder.fit_transform(dis).toarray())
    encoder_df[disease] = encoder_df[0]
    encoder_df = encoder_df.drop(columns=[0])
    encoder_df['genes'] = list(dis.index)
    encoder_df.set_index('genes', inplace=True)
    encoder_df.to_csv('../../../processed-data/gene-risk-LR-analysis/00-OpenTargets_SCZ_risk_genes/SCZ_risk_genes_01thr_hot1.csv')
    print(encoder_df)
else:
    print("Disease doesn't have enough genes with high genetic scores.")