#   Slightly modify the 'df.csv' files produced by 04-quantify_fluor.py for
#   compatibility with the loopy browser, which will allow manual annotation of
#   cell types, then used to train a decision tree classifier.
#
#   Despite the numbering ("09"-prepare_loopy), this script is intended to be
#   run after 04-quantify_fluor.py and before 06-evaluate_method.py.

import os
import numpy as np
import pandas as pd
import pyhere

################################################################################
#   Variable definitions
################################################################################

#-------------------------------------------------------------------------------
#   Paths
#-------------------------------------------------------------------------------

sample_info_path = pyhere.here(
    'raw-data', 'sample_info', 'Visium_IF_DLPFC_MasterExcel_01262022.xlsx'
)

df_in_path = pyhere.here(
    'processed-data', 'spot_deconvo', '02-cellpose', '{}', 'df.csv'
)

df_out_path = pyhere.here(
    'processed-data', 'spot_deconvo', '02-cellpose', '{}', 'df_loopy.csv'
)

################################################################################
#   Analysis
################################################################################

# os.environ['SGE_TASK_ID'] = '1'

#-------------------------------------------------------------------------------
#   Read in sample info and adjust paths for this particular sample ID
#-------------------------------------------------------------------------------

#   Different sample IDs are used for different files associated with each
#   sample. Determine both forms of the sample ID for this sample and update
#   path variables accordingly
sample_info = pd.read_excel(sample_info_path)[:4]
sample_ids_spot = 'Br' + sample_info['BrNumbr'].astype(int).astype(str) + \
    '_' + pd.Series([x.split('_')[1] for x in sample_info['Br_Region']]) + \
    '_IF'

sample_id_spot = sample_ids_spot[int(os.environ['SGE_TASK_ID']) - 1]
df_in_path = str(df_in_path).format(sample_id_spot)
df_out_path = str(df_out_path).format(sample_id_spot)

#   Read in the CSV, modify some column names, and re-write
df = pd.read_csv(df_in_path)
df.rename(
    {
        'x': 'y',
        'y': 'x',
        'Unnamed: 0': 'id'
    },
    axis = 1, inplace = True
)
df.to_csv(df_out_path, index = False)
