import pandas as pd
import numpy as np
import copy
import pickle

from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn import tree

import pyhere
from pathlib import Path

import matplotlib.pyplot as plt
import graphviz

df_path = pyhere.here(
    'processed-data', 'spot_deconvo', '02-cellpose', '{}', 'df_unfiltered.csv'
)

predictions_path = pyhere.here(
    'processed-data', 'spot_deconvo', '02-cellpose', '{}', 'model_predictions.csv'
)

dataset_path = pyhere.here(
    'processed-data', 'spot_deconvo', '02-cellpose', 'annotation_dataset.pkl'
)

tree_path = pyhere.here(
    'plots', 'spot_deconvo', '02-cellpose', 'cart', 'decision_tree.pdf'
)

model_out_path = pyhere.here(
    'processed-data', 'spot_deconvo', '02-cellpose', 'decision_tree.pkl'
)

sample_info_path = pyhere.here(
    'raw-data', 'sample_info', 'Visium_IF_DLPFC_MasterExcel_01262022.xlsx'
)

random_seed = 0

#   Number of cells per type to predict and output to Loopy-compatible CSV for
#   inspection
num_examples_check = 10

################################################################################
#   Functions
################################################################################

#   Remove splits that don't differentiate classes. Credit to:
#   https://github.com/scikit-learn/scikit-learn/issues/10810#issuecomment-409028102
def prune(tree):
    tree = copy.deepcopy(tree)
    dat = tree.tree_
    nodes = range(0, dat.node_count)
    ls = dat.children_left
    rs = dat.children_right
    classes = [[list(e).index(max(e)) for e in v] for v in dat.value]
    
    leaves = [(ls[i] == rs[i]) for i in nodes]
    
    LEAF = -1
    for i in reversed(nodes):
        if leaves[i]:
            continue
        if leaves[ls[i]] and leaves[rs[i]] and classes[ls[i]] == classes[rs[i]]:
            ls[i] = rs[i] = LEAF
            leaves[i] = True
    return tree

################################################################################
#   Analysis
################################################################################

#-------------------------------------------------------------------------------
#   Read in sample IDs and manual-annotation dataset
#-------------------------------------------------------------------------------

sample_info = pd.read_excel(sample_info_path)[:4]
sample_ids = 'Br' + sample_info['BrNumbr'].astype(int).astype(str) + \
    '_' + pd.Series([x.split('_')[1] for x in sample_info['Br_Region']]) + \
    '_IF'

with open(dataset_path, 'rb') as f:
    x_train, x_test, y_train, y_test = pickle.load(f)

#-------------------------------------------------------------------------------
#   Train the DecisionTreeClassifier
#-------------------------------------------------------------------------------

#   Instantiate, fit, and save the CART. Good hyperparameters were determined in
#   10-explore_models.py (via cross-validation, with final test accuracy of
#   99.0%) and are used here. Note that we still reserve part of the data as
#   test data, since we want to be able to make statements about the test
#   performance of the exact tree used for inference, and trees can
#   qualitatively change with the introduction of new data (if we used all
#   manual labels for inference).
model = tree.DecisionTreeClassifier(
    criterion = 'gini', 
    max_depth = 4,
    min_samples_leaf = 1, 
    random_state = random_seed,
    ccp_alpha = 0
)

model.fit(x_train, y_train)

with open(model_out_path, 'wb') as f:
    pickle.dump(model, f)

#   Compute training and test accuracy on the model
acc_train = round(100 * model.score(x_train, y_train), 1)
acc_test = round(100 * model.score(x_test, y_test), 1)
print(f'CART training accuracy: {acc_train}%.')
print(f'CART test accuracy: {acc_test}%.')

#   Print a more thorough report about training and test scores, making sure
#   performance is good across all classes (since we optimized for accuracy)
labels_train = model.predict(x_train)
labels_test = model.predict(x_test)
print('Training report:\n', classification_report(y_train, labels_train))
print('Test report:\n', classification_report(y_test, labels_test))

#-------------------------------------------------------------------------------
#   Write a small CSV of model predictions for import to loopy to visually
#   verify model performance
#-------------------------------------------------------------------------------

#   Read in the original unfiltered cells for a single sample
sample_id = sample_ids[0]
this_df = pd.read_csv(str(df_path).format(sample_id), index_col = 'id')
this_df = this_df.loc[:, ['gfap', 'neun', 'olig2', 'tmem119', 'area']]

#   Row IDs were not preserved during concatenation earlier, so we'll manually
#   drop rows in 'this_df' that are present in 'x_train' or 'x_test'
x = pd.concat([x_train, x_test])
this_df = this_df.merge(
    x, how='left', indicator=True,
    on = ['gfap', 'neun', 'olig2', 'tmem119', 'area']
)
this_df = this_df[this_df['_merge'] == 'left_only'].drop('_merge', axis=1)
this_df['label'] = model.predict(this_df)

#   For each cell type, pick [num_examples_check] cells randomly. Then form a
#   dataframe with these cells
small_df = pd.DataFrame()
for cell_type in this_df['label'].unique():
    indices = np.random.choice(
        this_df[this_df['label'] == cell_type].index,
        size = num_examples_check, replace = False
    )
    small_df = pd.concat([small_df, this_df.loc[indices]])

#   Make compatible with import to Loopy and write to CSV
small_df['id'] = small_df.index
small_df = small_df.loc[:, ['id', 'label']]
small_df.to_csv(str(predictions_path).format(sample_id), index = False)

#-------------------------------------------------------------------------------
#   Plot the (simplified) decision tree visually (save to PDF)
#-------------------------------------------------------------------------------

model = prune(model)
_ = tree.plot_tree(
    model, class_names = model.classes_, feature_names = x.columns,
    filled = True, rounded = True
)
plt.savefig(tree_path)
