#!/usr/bin/env python

###--------------------------------------------LOAD LIBRARIES
import pandas as pd
import numpy as np
import omnipath as op

###--------------------------------------------READ DATA
secondary = pd.read_csv('../../../processed-data/gene-risk-LR-analysis/02-SCZ_LRs/SCZ_LorR_interactions.csv')
liana = pd.read_csv('../../../processed-data/gene-risk-LR-analysis/02-SCZ_LRs/TableS3_liana_lr.csv')

ligands = []
receptors = []
for n in range(0,np.shape(secondary)[0]):
    if ((liana['ligand.complex']+liana['receptor.complex']).eq((secondary['genesymbol_intercell_source']+secondary['genesymbol_intercell_target']).loc[n]).sum()>0):
        ligands.append(secondary['genesymbol_intercell_source'].loc[n])
        receptors.append(secondary['genesymbol_intercell_target'].loc[n])
intersect = pd.DataFrame(np.concatenate((np.asarray(ligands).reshape(1,-1), np.asarray(receptors).reshape(1,-1)), axis = 0).T, columns = ['Receptors', 'Ligands'])
intersect.to_csv('../../../processed-data/gene-risk-LR-analysis/02-SCZ_LRs/SCZ_LorR_intersect_LIANA.csv')

db = op.interactions.import_intercell_network(transmitter_params = {"categories":"ligand"}, receiver_params = {"categories": "receptor"})
db = db[np.logical_not(db['genesymbol_intercell_source'].str.startswith('HLA'))]
db = db[np.logical_not(db['genesymbol_intercell_target'].str.startswith('HLA'))]

tgts = pd.read_csv('../../../processed-data/gene-risk-LR-analysis/00-OpenTargets_SCZ_risk_genes/SCZ_risk_genes_01thr_hot1.csv')

###--------------------------------------------KEEP ALL INTERACTIONS WHERE AT LEAST ONE PART IS ASSOCIATED WITH SCZ RISK
ligand_SCZ = []
receptor_SCZ = []
for m in range(0,np.shape(intersect)[0]):
    if (m==0):
        op_LRs = pd.DataFrame(db[(db['genesymbol_intercell_source']==intersect.iloc[m][0]) & (db['genesymbol_intercell_target']==intersect.iloc[m][1])]).reset_index()
        if intersect.iloc[m][0] in list(tgts['genes']):
            ligand_SCZ.append(True)
        else:
            ligand_SCZ.append(False)
        if intersect.iloc[m][1] in list(tgts['genes']):
            receptor_SCZ.append(True)
        else:
            receptor_SCZ.append(False)
    else:
        new = pd.DataFrame(db[(db['genesymbol_intercell_source']==intersect.iloc[m][0]) & (db['genesymbol_intercell_target']==intersect.iloc[m][1])]).reset_index()
        op_LRs = pd.concat([op_LRs,new], axis = 0)
        if intersect.iloc[m][0] in list(tgts['genes']):
            ligand_SCZ.append(True)
        else:
            ligand_SCZ.append(False)
        if intersect.iloc[m][1] in list(tgts['genes']):
            receptor_SCZ.append(True)
        else:
            receptor_SCZ.append(False)
op_LRs['receptor_SCZ_risk'] = receptor_SCZ
op_LRs['ligand_SCZ_risk'] = ligand_SCZ

###--------------------------------------------WRITE REORDERED FILE
op_LRs[['genesymbol_intercell_source','genesymbol_intercell_target','receptor_SCZ_risk',
       'ligand_SCZ_risk',
       'index', 'source', 'target', 'is_stimulation', 'is_inhibition',
       'consensus_direction', 'consensus_stimulation', 'consensus_inhibition',
       'curation_effort', 'references', 'sources', 'type',
       'references_stripped', 'n_references', 'n_sources', 'n_primary_sources',
       'category_intercell_source', 'parent_intercell_source',
       'database_intercell_source', 'scope_intercell_source',
       'aspect_intercell_source', 'category_source_intercell_source',
       'uniprot_intercell_source', 
       'entity_type_intercell_source', 'consensus_score_intercell_source',
       'transmitter_intercell_source', 'receiver_intercell_source',
       'secreted_intercell_source',
       'plasma_membrane_transmembrane_intercell_source',
       'plasma_membrane_peripheral_intercell_source',
       'category_intercell_target', 'parent_intercell_target',
       'database_intercell_target', 'scope_intercell_target',
       'aspect_intercell_target', 'category_source_intercell_target',
       'uniprot_intercell_target', 
       'entity_type_intercell_target', 'consensus_score_intercell_target',
       'transmitter_intercell_target', 'receiver_intercell_target',
       'secreted_intercell_target',
       'plasma_membrane_transmembrane_intercell_target',
       'plasma_membrane_peripheral_intercell_target']].drop(columns = 'index').to_csv('../../../processed-data/gene-risk-LR-analysis/02-SCZ_LRs/SCZ_LorR_intersect_LIANA_op_annotation.csv', index = False)