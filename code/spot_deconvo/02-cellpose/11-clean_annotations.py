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

manual_label_path_orig_in = pyhere.here(
    'processed-data', 'spot_deconvo', '02-cellpose', '{}',
    'manual_labels_raw.csv'
)
manual_label_path_conf_in = pyhere.here(
    'processed-data', 'spot_deconvo', '02-cellpose', '{}',
    'manual_labels_confidence_raw.csv'
)

manual_label_path_orig_out = pyhere.here(
    'processed-data', 'spot_deconvo', '02-cellpose', '{}',
    'manual_labels_clean.csv'
)
manual_label_path_conf_out = pyhere.here(
    'processed-data', 'spot_deconvo', '02-cellpose', '{}',
    'manual_labels_confidence_clean.csv'
)

confidence_path = pyhere.here(
    'processed-data', 'spot_deconvo', '02-cellpose', '{}', 'confidence_unfiltered.csv'
)

sample_info_path = pyhere.here(
    'raw-data', 'sample_info', 'Visium_IF_DLPFC_MasterExcel_01262022.xlsx'
)

#   Define expected cell-type labels (after cleaning), expected columns
#   in the fluorescence-intensity data frame, and expected number of counts
#   for each cell-type label
expected_labels = ["astro", "micro", "neuron", "oligo", "other"]
expected_df_cols = [
    'id', 'gfap', 'neun', 'olig2', 'tmem119', 'area', 'y', 'x', 'dist', 'idx'
]
expected_label_counts = 30

#   Define original labels and what they should be replaced with
replace_from = ["Astrocytes", "Astrocyte", "Microglia", "Neurons", "Neuron", "Oligo", "Other"]
replace_to = ["astro", "astro", "micro", "neuron", "neuron", "oligo", "other"]

################################################################################
#   Clean data
################################################################################

#   Read in the list of sample IDs for the ID data
sample_info = pd.read_excel(sample_info_path)[:4]
sample_ids = 'Br' + sample_info['BrNumbr'].astype(int).astype(str) + \
    '_' + pd.Series([x.split('_')[1] for x in sample_info['Br_Region']]) + \
    '_IF'

for sample_id in sample_ids:
    #   Determine paths for this sample
    this_df_path = str(df_path).format(sample_id)
    this_manual_label_path_orig_in = str(manual_label_path_orig_in).format(sample_id)
    this_manual_label_path_orig_out = str(manual_label_path_orig_out).format(sample_id)
    this_manual_label_path_conf_in = str(manual_label_path_conf_in).format(sample_id)
    this_manual_label_path_conf_out = str(manual_label_path_conf_out).format(sample_id)
    this_confidence_path = str(confidence_path).format(sample_id)
    
    #   Read in all required CSVs
    this_df = pd.read_csv(this_df_path)
    this_manual_labels_orig = pd.read_csv(this_manual_label_path_orig_in)
    this_manual_labels_conf = pd.read_csv(this_manual_label_path_conf_in)
    this_confidence = pd.read_csv(this_confidence_path)
    
    #   Check the columns are all as expected
    assert all(
        [x == y for x, y in zip(this_df.columns.tolist(), expected_df_cols)]
    )
    
    #   Fix indices (index rows by cell ID)
    this_manual_labels_orig.index = this_manual_labels_orig['id']
    this_manual_labels_conf.index = this_manual_labels_conf['id']
    this_df.index = this_df['id']
    this_confidence.index = this_confidence['id']
    
    #   Clean up "confidence" manual labels
    this_manual_labels_conf.dropna(inplace = True)
    this_manual_labels_conf['label'] = [x.split('_')[0] for x in this_manual_labels_conf['label']]
    assert set(this_manual_labels_conf['label']) == set(expected_labels)
    print(f'New cell-type counts (confidence) for sample {sample_id}: \n{this_manual_labels_conf["label"].value_counts()}')
    
    #   Add confidence quantile (and old cell-type call) back from the original CSV
    #   to annotate
    assert set(this_manual_labels_conf.index) == set(this_confidence.index)
    this_manual_labels_conf['quantile'] = this_confidence['quantile']
    this_manual_labels_conf['label_old'] = this_confidence['label']
    
    #   Drop non-labelled spots and check that labels are expected, after
    #   replacing some known values
    print(f"Dropping {sum(this_manual_labels_orig['label'].isna())} NA labels from sample {sample_id}.")
    this_manual_labels_orig.dropna(inplace = True)
    this_manual_labels_orig.replace(replace_from, replace_to, inplace = True)
    assert set(expected_labels) == set(this_manual_labels_orig['label'])
    
    #   Verify the correct number of labels per cell type is present, and all cell
    #   IDs line up to the original data frame of fluorescence intensities
    assert all([x == expected_label_counts for x in this_manual_labels_orig['label'].value_counts()])
    assert(all(this_manual_labels_orig['id'].isin(this_df['id'])))
    assert(all(this_manual_labels_conf['id'].isin(this_df['id'])))
    
    #   Write a clean copy of both sets of manual labels
    this_manual_labels_orig.to_csv(this_manual_label_path_orig_out, index = False)
    this_manual_labels_conf.to_csv(this_manual_label_path_conf_out, index = False)
