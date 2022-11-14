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
import pickle

df_path = pyhere.here(
    'processed-data', 'spot_deconvo', '02-cellpose', '{}', 'df_unfiltered.csv'
)

dataset_path = pyhere.here(
    'processed-data', 'spot_deconvo', '02-cellpose', 'annotation_dataset.pkl'
)

manual_label_orig_path = pyhere.here(
    'processed-data', 'spot_deconvo', '02-cellpose', '{}', 'manual_labels_clean.csv'
)

manual_label_conf_path = pyhere.here(
    'processed-data', 'spot_deconvo', '02-cellpose', '{}', 'manual_labels_confidence_clean.csv'
)

sample_info_path = pyhere.here(
    'raw-data', 'sample_info', 'Visium_IF_DLPFC_MasterExcel_01262022.xlsx'
)

expected_num_labels = 30
num_cell_types = 5

test_proportion = 0.2 # for training/test split

random_seed = 0

################################################################################
#   Analysis
################################################################################

#-------------------------------------------------------------------------------
#   Preprocess and gather fluorescence data + original manual cell-type labels
#-------------------------------------------------------------------------------

sample_info = pd.read_excel(sample_info_path)[:4]
sample_ids = 'Br' + sample_info['BrNumbr'].astype(int).astype(str) + \
    '_' + pd.Series([x.split('_')[1] for x in sample_info['Br_Region']]) + \
    '_IF'

#   Loop through labels and flourescence intensity tables for each sample and
#   combine into a single data frame
df = pd.DataFrame()
for sample_id in sample_ids:
    this_df_path = str(df_path).format(sample_id)
    this_manual_label_orig_path = str(manual_label_orig_path).format(sample_id)
    
    this_df = pd.read_csv(this_df_path, index_col = 'id')
    this_manual_labels_orig = pd.read_csv(this_manual_label_orig_path, index_col = 'id')
    
    this_df['label'] = this_manual_labels_orig['label']
    this_df['label_sample'] = this_df['label'] + '_' + sample_id
    df = pd.concat([df, this_df.dropna()])

#   Verify we have the correct amount of cells
assert(df.shape[0] == len(sample_ids) * expected_num_labels * num_cell_types)

#   Define the inputs (features we want the model to access) and outputs to the
#   model
x = df.loc[:, ['gfap', 'neun', 'olig2', 'tmem119', 'area', 'label_sample']]
y = df['label']

#   Split data into training and test sets (80%: 20%), evenly stratified across
#   classes and sample ID
x_train, x_test, y_train, y_test = train_test_split(
    x, y, test_size = test_proportion, random_state = random_seed,
    stratify = x['label_sample']
)

#   Verify the stratification worked (note the exact equality only works because
#   of the nice divisibility of our data size)
assert all(
    [
        x == expected_num_labels * (1 - test_proportion)
            for x in x_train['label_sample'].value_counts()
    ]
)
assert all(
    [
        x == expected_num_labels * test_proportion
            for x in x_test['label_sample'].value_counts()
    ]
)

#   Remove the column we only used to properly stratify the data
x_train.drop('label_sample', axis = 1, inplace = True)
x_test.drop('label_sample', axis = 1, inplace = True)

#-------------------------------------------------------------------------------
#   Preprocess and gather fluorescence data + additional manual annotation
#   based on a variety of previously unlabelled examples across variable model
#   confidence
#-------------------------------------------------------------------------------

df = pd.DataFrame()
for sample_id in sample_ids:
    this_df_path = str(df_path).format(sample_id)
    this_manual_label_conf_path = str(manual_label_conf_path).format(sample_id)
    
    this_df = pd.read_csv(this_df_path, index_col = 'id')
    this_manual_labels_conf = pd.read_csv(
        this_manual_label_conf_path, index_col = 'id'
    )
    
    #   Form a column combining old cell-type label, confidence quantile, and
    #   sample ID
    this_df['group'] = this_manual_labels_conf['label_old'] + '_' + this_manual_labels_conf['quantile'].astype(str) + '_' + sample_id
    this_df['label'] = this_manual_labels_conf['label']
    df = pd.concat([df, this_df.dropna()])

#   There should be 4 annotated cells per sample per cell type per quantile
assert all([x == 4 for x in df['group'].value_counts()])

#   Define the inputs (features we want the model to access) and outputs to the
#   model
x = df.loc[:, ['gfap', 'neun', 'olig2', 'tmem119', 'area', 'group']]
y = df['label']

x_train2, x_test2, y_train2, y_test2 = train_test_split(
    x, y, test_size = 0.25, random_state = random_seed,
    stratify = x['group']
)

#   Verify the stratification worked (note the exact equality only works because
#   of the nice divisibility of our data size)
assert all([x == 3 for x in x_train2['group'].value_counts()])
assert all([x == 1 for x in x_test2['group'].value_counts()])

#   Remove the column we only used to properly stratify the data
x_train2.drop('group', axis = 1, inplace = True)
x_test2.drop('group', axis = 1, inplace = True)

#-------------------------------------------------------------------------------
#   Combine the two different sets of manual annotation into one dataset and
#   write to disk
#-------------------------------------------------------------------------------

x_train = pd.concat([x_train, x_train2])
x_test = pd.concat([x_test, x_test2])
y_train = pd.concat([y_train, y_train2])
y_test = pd.concat([y_test, y_test2])

perc_train = round(100 * x_train.shape[0] / (x_train.shape[0] + x_test.shape[0]), 2)
print(f'Using {x_train.shape[0]} training and {x_test.shape[0]} test examples ({perc_train}% training).')

#   Write the dataset to disk for later use
with open(dataset_path_out, 'wb') as f:
    pickle.dump((x_train, x_test, y_train, y_test), f)

#-------------------------------------------------------------------------------
#   Try a decision tree classifier
#-------------------------------------------------------------------------------

tuned_parameters = [
    {
        'criterion': ['gini', 'entropy'],
        'max_depth': [2, 3, 4, 5],
        'min_samples_leaf': [1, 5, 10, 15, 20, 25],
        'ccp_alpha': [0, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3]
    }
]

grid = GridSearchCV(
    tree.DecisionTreeClassifier(
        random_state = random_seed, splitter = 'best', class_weight = None
    ),
    tuned_parameters, cv=5, scoring = 'accuracy'
)
grid.fit(x_train, y_train)

#   Compute training and test accuracy on the best model
acc_train = round(100 * grid.best_estimator_.score(x_train, y_train), 1)
acc_test = round(100 * grid.best_estimator_.score(x_test, y_test), 1)
print(f'CART training accuracy: {acc_train}%.')
print(f'CART test accuracy: {acc_test}%.')
print(f'Best params: {grid.best_params_}')

#   Print a more thorough report about training and test scores, making sure
#   performance is good across all classes (since we optimized for accuracy)
labels_train = grid.best_estimator_.predict(x_train)
labels_test = grid.best_estimator_.predict(x_test)
print('Training report:\n', classification_report(y_train, labels_train))
print('Test report:\n', classification_report(y_test, labels_test))

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

acc_train = round(100 * grid.best_estimator_.score(x_train, y_train), 1)
acc_test = round(100 * grid.best_estimator_.score(x_test, y_test), 1)
print(f'Logistic regression training accuracy: {acc_train}%.')
print(f'Logistic regression test accuracy: {acc_test}%.')
print(f'Best params: {grid.best_params_}')

#   Print a more thorough report about training and test scores, making sure
#   performance is good across all classes (since we optimized for accuracy)
labels_train = grid.best_estimator_.predict(x_train)
labels_test = grid.best_estimator_.predict(x_test)
print('Training report:\n', classification_report(y_train, labels_train))
print('Test report:\n', classification_report(y_test, labels_test))

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

acc_train = round(100 * grid.best_estimator_.score(x_train, y_train), 1)
acc_test = round(100 * grid.best_estimator_.score(x_test, y_test), 1)
print(f'SVM training accuracy: {acc_train}%.')
print(f'SVM test accuracy: {acc_test}%.')
print(f'Best params: {grid.best_params_}')

#   Print a more thorough report about training and test scores, making sure
#   performance is good across all classes (since we optimized for accuracy)
labels_train = grid.best_estimator_.predict(x_train)
labels_test = grid.best_estimator_.predict(x_test)
print('Training report:\n', classification_report(y_train, labels_train))
print('Test report:\n', classification_report(y_test, labels_test))
