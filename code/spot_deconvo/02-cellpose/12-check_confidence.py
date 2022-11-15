import pandas as pd
import numpy as np

from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression

import pyhere
from pathlib import Path
import pickle

df_path = pyhere.here(
    'processed-data', 'spot_deconvo', '02-cellpose', '{}', 'df_unfiltered.csv'
)

confidence_path = pyhere.here(
    'processed-data', 'spot_deconvo', '02-cellpose', '{}', 'confidence_unfiltered.csv'
)

dataset_path = pyhere.here(
    'processed-data', 'spot_deconvo', '02-cellpose', 'annotation_dataset.pkl'
)

manual_label_path = pyhere.here(
    'processed-data', 'spot_deconvo', '02-cellpose', '{}', 'manual_labels_clean.csv'
)

sample_info_path = pyhere.here(
    'raw-data', 'sample_info', 'Visium_IF_DLPFC_MasterExcel_01262022.xlsx'
)

random_seed = 0

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

#   The decision tree didn't use area as a feature, so the LogisticRegression
#   model shouldn't either
x_train.drop('area', axis = 1, inplace = True)
x_test.drop('area', axis = 1, inplace = True)

#-------------------------------------------------------------------------------
#   Train logistic-regression model and verify performance
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

#-------------------------------------------------------------------------------
#   Compare confidences on test set vs. unlabelled examples
#-------------------------------------------------------------------------------

#   Read in all the original unfiltered cells for all samples
this_df = pd.DataFrame()
for sample_id in sample_ids:
    this_df_small = pd.read_csv(str(df_path).format(sample_id), index_col = 'id')
    this_df_small['cell_id'] = this_df_small.index
    this_df_small['sample_id'] = sample_id
    this_df_small = this_df_small.loc[
        :, ['gfap', 'neun', 'olig2', 'tmem119', 'sample_id', 'cell_id']
    ]
    this_df = pd.concat([this_df, this_df_small])

#   Merge the DataFrame of all cells with that of just the test cells. This is
#   just a trick to categorize which of the unfiltered cells belong to the test
#   set vs. were never manually annotated
x_test['cell_id'] = x_test.index
this_df = this_df.merge(
    x_test, how='left', indicator=True,
    on = ['gfap', 'neun', 'olig2', 'tmem119', 'cell_id']
)

#   Use more informative value names
this_df['_merge'].cat.remove_unused_categories(inplace=True)
this_df['_merge'].cat.rename_categories(
    {'left_only': 'unlabeled', 'both': 'test'},
    inplace=True
)

#   Ensure we get the correct number of examples in the test set
assert this_df['_merge'].value_counts()['test'] == x_test.shape[0]

#   Calculate confidences (maximal "probability" across classes for each
#   example) by predicted cell type and by group (test set or unlabelled)
conf_df = pd.DataFrame(
    {
        'label': grid.best_estimator_.predict(
            this_df.drop(['_merge', 'sample_id', 'cell_id'], axis = 1)
        ),
        'confidence': np.max(
            grid.best_estimator_.predict_proba(
                this_df.drop(['_merge', 'sample_id', 'cell_id'], axis = 1)
            ),
            axis = 1
        ),
        'group': this_df['_merge'],
        'sample_id': this_df['sample_id'],
        'cell_id': this_df['cell_id']
    }
)

#   While not necessarily perfectly representative of the robustness of the
#   DecisionTreeClassifier, the LogisticRegression model is notably more
#   confident on test-set examples than unlabelled examples (for all cell
#   types).
print('Average model confidences on unlabeled examples vs. labelled examples in the test set:')
print(conf_df.groupby(['group', 'label'])['confidence'].mean())

#-------------------------------------------------------------------------------
#   Generate a new CSV of cells to annotate based on their confidences
#-------------------------------------------------------------------------------

#   Take only the unlabeled cells
conf_df = conf_df.loc[conf_df.group == 'unlabeled'].drop('group', axis = 1)

#   Compute the boundaries of each quantile of the confidences on cells, by
#   group (combination of cell type and sample ID), and add to the DF
quant_df = conf_df \
    .groupby(['label', 'sample_id'])['confidence'] \
    .apply(np.quantile, q = [0.25, 0.5, 0.75]) \
    .reset_index(name = 'q')

conf_df = conf_df.merge(quant_df, how='left', on = ['label', 'sample_id'])

#   Assign which quantile each cell belongs to (in {1, 2, 3, 4})
conf_df['quantile'] = 0
conf_df.loc[conf_df.apply(lambda x: x[1] < x[4][0], axis = 1), 'quantile'] = 1
conf_df.loc[
    conf_df.apply(lambda x: (x[1] >= x[4][0]) & (x[1] < x[4][1]), axis = 1),
    'quantile'
] = 2
conf_df.loc[
    conf_df.apply(lambda x: (x[1] >= x[4][1]) & (x[1] < x[4][2]), axis = 1),
    'quantile'
] = 3
conf_df.loc[conf_df.apply(lambda x: x[1] >= x[4][2], axis = 1), 'quantile'] = 4

#   Sanity check: verify each quantile is nearly equal in size
print('Number of examples in each quantile (sanity check):')
print(conf_df['quantile'].value_counts())

#   Randomly take 4 cells from each group and quantile
conf_df = conf_df \
    .groupby(['label', 'sample_id', 'quantile']) \
    .sample(n = 4) \
    .drop('q', axis = 1)

#   Write a CSV for import to loopy for each sample. These cells can be
#   relabeled if necessary, and the model accordingly retrained
# for sample_id in sample_ids:
#     small_df = conf_df \
#         .loc[conf_df['sample_id'] == sample_id, :] \
#         .rename({'cell_id': 'id'}, axis = 1) \
#         .drop('sample_id', axis = 1)

#     small_df.to_csv(str(confidence_path).format(sample_id), index = False)
