#   Standardize and clean manual cell-type labels and re-write CSVs for
#   immediate downstream use
#
#   1. Drop NA labels and use consistent names for cell types
#   2. Ensure exactly 30 labels per cell type per sample

import pandas as pd
import numpy as np
import pyhere
from pathlib import Path

################################################################################
#   Variable definitions
################################################################################

df_path = pyhere.here(
    'processed-data', 'spot_deconvo', '02-cellpose', '{}', 'df_unfiltered.csv'
)

manual_label_path_in = pyhere.here(
    'processed-data', 'spot_deconvo', '02-cellpose', '{}',
    'manual_labels_raw.csv'
)

manual_label_path_out = pyhere.here(
    'processed-data', 'spot_deconvo', '02-cellpose', '{}',
    'manual_labels_clean.csv'
)

sample_info_path = pyhere.here(
    'raw-data', 'sample_info', 'Visium_IF_DLPFC_MasterExcel_01262022.xlsx'
)

#   Define expected cell-type labels (after cleaning), expected columns
#   in the fluorescence-intensity data frame, and expected number of counts
#   for each cell-type label
expected_labels = ["astro", "micro", "neuron", "oligo"]
expected_df_cols = [
    'id', 'gfap', 'neun', 'olig2', 'tmem119', 'area', 'y', 'x', 'dist', 'idx'
]
expected_label_counts = 30

#   Define original labels and what they should be replaced with
replace_from = ["Astrocytes", "Astrocyte", "Microglia", "Neurons", "Neuron", "Oligo"]
replace_to = ["astro", "astro", "micro", "neuron", "neuron", "oligo"]

################################################################################
#   Clean data
################################################################################

#   Read in the list of sample IDs for the ID data
sample_info = pd.read_excel(sample_info_path)[:4]
sample_ids = 'Br' + sample_info['BrNumbr'].astype(int).astype(str) + \
    '_' + pd.Series([x.split('_')[1] for x in sample_info['Br_Region']]) + \
    '_IF'

for i in range(len(sample_ids)):
    #   Determine paths for this sample
    this_df_path = str(df_path).format(sample_ids[i])
    this_manual_label_path_in = str(manual_label_path_in).format(sample_ids[i])
    this_manual_label_path_out = str(manual_label_path_out).format(sample_ids[i])
    
    this_df = pd.read_csv(this_df_path)
    this_manual_labels = pd.read_csv(this_manual_label_path_in)
    
    #   Check the columns are all as expected
    assert all(
        [x == y for x, y in zip(this_df.columns.tolist(), expected_df_cols)]
    )
    
    #   Fix indices
    this_manual_labels.index = this_manual_labels['id']
    this_df.index = this_df['id']
    
    #   Drop non-labelled spots and check that labels are expected, after
    #   replacing some known values
    print(f"Dropping {sum(this_manual_labels['label'].isna())} NA labels from sample {sample_ids[i]}.")
    this_manual_labels.dropna(inplace = True)
    this_manual_labels.replace(replace_from, replace_to, inplace=True)
    assert set(expected_labels) == set(this_manual_labels['label'])
    
    #   One sample has 2 extra 'astro' labels. Randomly drop 2 'astro' labels to
    #   get the count for each cell type to 30
    if sample_ids[i] == 'Br2720_Ant_IF':
        indices = np.random.choice(
            this_manual_labels['label'].loc[this_manual_labels['label'] == 'astro'].index,
            size=2, replace=False
        )
        this_manual_labels.drop(index = indices, inplace = True)
    
    #   Verify the correct number of labels per cell type is present, and all cell
    #   IDs line up to the original data frame of fluorescence intensities
    assert all([x == expected_label_counts for x in this_manual_labels['label'].value_counts()])
    assert(all(this_manual_labels['id'].isin(this_df['id'])))
    
    #   Write a clean copy of the manual labels
    this_manual_labels.to_csv(this_manual_label_path_out, index = False)
