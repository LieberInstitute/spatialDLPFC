#   We'll use the Decision Tree Classifier trained in 07-cart.py for the paper.
#   This script more thoroughly explores alternative models (such as logistic
#   regression and various SVM approaches), to make sure we aren't missing out
#   on a vastly superior model

import pandas as pd
import numpy as np

from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.metrics import classification_report
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn import tree, svm

import pyhere
from pathlib import Path

sample_name = 'Br2720_Ant_IF'

df_path = pyhere.here(
    'processed-data', 'spot_deconvo', '02-cellpose', sample_name, 'df.csv'
)

manual_label_path = pyhere.here(
    'processed-data', 'spot_deconvo', '02-cellpose', sample_name, 'manual_labels.csv'
)

tree_path = pyhere.here(
    'plots', 'spot_deconvo', '02-cellpose', 'cart', sample_name + '_decision_tree.pdf'
)

model_out_path = pyhere.here(
    'processed-data', 'spot_deconvo', '02-cellpose', sample_name, 'decision_tree.pkl'
)

cell_types = {
    "neun": "neuron",
    "olig2": "oligo",
    "tmem119": "micro"
}

random_seed = 0

################################################################################
#   Analysis
################################################################################

#-------------------------------------------------------------------------------
#   Preprocess and clean fluorescence data + manual cell-type labels
#-------------------------------------------------------------------------------

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

#   Split data into training, test, and validation sets (70%: 15%: 15%), evenly
#   stratified across classes. A validation set is necessary because we'll do
#   a couple grid searchs over hyperparameters, optimizing for test performance
x_train, x_test, y_train, y_test = train_test_split(
    x, y, test_size = 0.3, random_state = random_seed, stratify = y
)

x_test, x_val, y_test, y_val = train_test_split(
    x_test, y_test, test_size = 0.5, random_state = random_seed, stratify = y_test
)

#-------------------------------------------------------------------------------
#   Try a decision tree classifier
#-------------------------------------------------------------------------------

tuned_parameters = [
    {
        'criterion': ['gini', 'entropy'],
        'max_depth': [2, 3, 4, 5],
        'min_samples_leaf': [1, 5, 10, 15, 20, 25],
        'ccp_alpha': [0.001, 0.01, 0.1, 1]
    }
]

grid = GridSearchCV(
    tree.DecisionTreeClassifier(
        random_state = random_seed, splitter = 'best', class_weight = None
    ),
    tuned_parameters, cv=5, scoring = 'accuracy'
)
grid.fit(x_train, y_train)

#   Train the decision tree using the best params and evaluate on validation set
model = tree.DecisionTreeClassifier(**grid.best_params_)
model.fit(x_train, y_train)

#   Compute training and validation accuracy
acc_train = round(100 * model.score(x_train, y_train), 1)
acc_test = round(100 * model.score(x_val, y_val), 1)
print(f'CART training accuracy: {acc_train}%.')
print(f'CART validation accuracy: {acc_test}%.')

#-------------------------------------------------------------------------------
#   Try logistic regression
#-------------------------------------------------------------------------------

tuned_parameters = [
    {
        'logisticregression__C': np.logspace(-2, 4, 20)
    }
]

pipe = make_pipeline(
    StandardScaler(),
    LogisticRegression(random_state = random_seed, max_iter=10000)
)

grid = GridSearchCV(pipe, tuned_parameters, cv=5, scoring = 'accuracy')
grid.fit(x_train, y_train)

pipe = make_pipeline(
    StandardScaler(),
    LogisticRegression(
        random_state = random_seed,
        C = grid.best_params_['logisticregression__C']
    )
)
pipe.fit(x_train, y_train)

acc_train = round(100 * pipe.score(x_train, y_train), 1)
acc_test = round(100 * pipe.score(x_val, y_val), 1)
print(f'Logistic regression training accuracy: {acc_train}%.')
print(f'Logistic regression validation accuracy: {acc_test}%.')

#-------------------------------------------------------------------------------
#   Try SVM (linear and non-linear kernels)
#-------------------------------------------------------------------------------

tuned_parameters = [
    {
        'svc__kernel': ['linear', 'rbf'],
        'svc__gamma': np.logspace(-5, -1, 10),
        'svc__C': np.logspace(0, 4, 5)
    }
]

pipe = make_pipeline(
    StandardScaler(),
    svm.SVC(random_state = random_seed)
)

grid = GridSearchCV(pipe, tuned_parameters, cv=5, scoring = 'accuracy')
grid.fit(x_train, y_train)

pipe = make_pipeline(
    StandardScaler(),
    svm.SVC(
        random_state = random_seed,
        C = grid.best_params_['svc__C'],
        kernel = grid.best_params_['svc__kernel'],
        gamma = grid.best_params_['svc__gamma']
    )
)
pipe.fit(x_train, y_train)

acc_train = round(100 * pipe.score(x_train, y_train), 1)
acc_test = round(100 * pipe.score(x_val, y_val), 1)
print(f'SVM training accuracy: {acc_train}%.')
print(f'SVM validation accuracy: {acc_test}%.')

