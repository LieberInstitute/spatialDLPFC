import pandas as pd
import numpy as np
import copy

from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report
from sklearn import tree

import pyhere
from pathlib import Path

import matplotlib.pyplot as plt
import graphviz

sample_name = 'Br2720_Ant_IF'

df_path = pyhere.here(
    'processed-data', 'spot_deconvo', '02-cellpose', sample_name, 'df.csv'
)

manual_label_path = pyhere.here(
    'processed-data', 'spot_deconvo', '02-cellpose', sample_name, 'manual_labels.csv'
)

tree_path = pyhere.here(
    'plots', 'spot_deconvo', '02-cellpose', sample_name + '_decision_tree.pdf'
)

cell_types = {
    "neun": "neuron",
    "olig2": "oligo",
    "tmem119": "micro"
}

#   Hyperparameters for CART
random_seed = 0
min_leaf = 10
max_depth = 4
criterion = 'entropy' # 'gini'

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
    #
    leaves = [(ls[i] == rs[i]) for i in nodes]
    #
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

#   Read in list of cells and format
df = pd.read_csv(df_path)
df.rename({'Unnamed: 0': 'id'}, axis = 1, inplace = True)
df.index = df['id']

#   Read in manual labels and subset list of cells to the ones labelled
manual_labels = pd.read_csv(manual_label_path)
manual_labels.index = manual_labels['id']
manual_labels.loc[manual_labels['value'] == 'microglia', 'value'] = 'micro'
df = df.loc[manual_labels['id']]

#   Subset to the features we want the CART to access
# x = df.loc[:, ['gfap', 'neun', 'olig2', 'tmem119', 'area']]
x = df.drop(['id', 'x', 'y', 'dist', 'idx'], axis = 1)
y = manual_labels['value']

#   Split data into training and test sets (20% test), evenly stratified
#   across classes
x_train, x_test, y_train, y_test = train_test_split(
    x, y, test_size = 0.2, random_state = random_seed, stratify = y
)

#   Instantiate and fit the CART
model = tree.DecisionTreeClassifier(
    criterion = criterion, 
    splitter = 'best', 
    max_depth = max_depth,
    class_weight = None,
    min_samples_leaf = min_leaf, 
    random_state = random_seed,
    ccp_alpha = 0.02
)

clf = model.fit(x_train, y_train)

#   Compute training and test accuracy
acc_train = round(100 * model.score(x_train, y_train), 1)
acc_test = round(100 * model.score(x_test, y_test), 1)
print(f'Training accuracy: {acc_train}%.')
print(f'Test accuracy: {acc_test}%.')

#   Print a more thorough report about training and test scores
labels_train = model.predict(x_train)
labels_test = model.predict(x_test)
print(classification_report(y_test, labels_test))
print(classification_report(y_train, labels_train))

#   Plot the (simplified) decision tree visually (save to PDF)
prunedTree = prune(model)
tree.plot_tree(
    prunedTree, class_names = clf.classes_, feature_names = x.columns,
    filled = True, rounded = True
)
plt.savefig(tree_path)
