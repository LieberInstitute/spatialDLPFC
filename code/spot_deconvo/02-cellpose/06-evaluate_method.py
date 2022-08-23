import os
import matplotlib.pyplot as plt
import numpy as np
from numpy.random import default_rng
import pandas as pd
import seaborn as sns
import pyhere
from pathlib import Path

sample_name = 'Br2720_Ant_IF'

df_path = pyhere.here(
    'processed-data', 'spot_deconvo', '02-cellpose', sample_name, 'df.csv'
)

manual_label_path = pyhere.here(
    'processed-data', 'spot_deconvo', '02-cellpose', sample_name, 'manual_labels.csv'
)

cell_types = {
    "neun": "neuron",
    "olig2": "oligo",
    "tmem119": "micro"
}

#   Mean fluorescence below this value will be considered an indication of an 
#   absence of the given cell type
noise_cutoff = 3

def classify_cells(df, cell_types, verbose = False):
    weights = {}
    for channel in cell_types.keys():
        if verbose:
            print(f'Mean for {channel} above cutoff: {np.mean(df[channel][df[channel] > noise_cutoff])}')
            print(f'Mean for {channel} otherwise: {np.mean(df[channel])}')
        
        weights[channel] = 1 / np.mean(df[channel][df[channel] > noise_cutoff])
    
    
    #   Weight the intensity of each channel inversely by the average intensity seen
    #   in cells above the "noise cutoff". Then call cell types by choosing the
    #   maximal weighted intensity across the 4 possible types
    temp = df[[x for x in cell_types.keys()]].copy()
    for channel in cell_types.keys():
        temp[channel] *= weights[channel]
    
    temp['max_channel'] = [temp.columns[np.argmax(temp.loc[i])] for i in temp.index]
    df['cell_type'] = [cell_types[x] for x in temp['max_channel']]
    
    #   For now, assume all cells with low fluorescence in all channels are a cell
    #   type not measured by a marker/ image channel
    df.loc[
        (df['neun'] < noise_cutoff) &
        (df['olig2'] < noise_cutoff) &
        (df['tmem119'] < noise_cutoff),
        'cell_type'
    ] = 'other'
    
    return temp

#   Read in list of cells, format, and classify cell types
df = pd.read_csv(df_path)
df.rename({'Unnamed: 0': 'id'}, axis = 1, inplace = True)
df.index = df['id']
_ = classify_cells(df, cell_types, True)

#   Read in manual labels and subset list of cells to the ones labelled
manual_labels = pd.read_csv(manual_label_path)
manual_labels.index = manual_labels['id']
df = df.loc[manual_labels['id']]

#   Evaluate classification accuracy
assert set(df['cell_type']) == set(manual_labels['value'])

accuracy = 100 * np.count_nonzero(df['cell_type'] == manual_labels['value']) / df.shape[0]
print(f'Classification accuracy: {accuracy}%.')
