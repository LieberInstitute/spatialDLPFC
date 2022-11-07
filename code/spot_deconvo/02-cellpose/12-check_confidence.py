import pandas as pd
import numpy as np

from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression

import pyhere
from pathlib import Path

df_path = pyhere.here(
    'processed-data', 'spot_deconvo', '02-cellpose', '{}', 'df_unfiltered.csv'
)

manual_label_path = pyhere.here(
    'processed-data', 'spot_deconvo', '02-cellpose', '{}', 'manual_labels_clean.csv'
)

sample_info_path = pyhere.here(
    'raw-data', 'sample_info', 'Visium_IF_DLPFC_MasterExcel_01262022.xlsx'
)

expected_num_labels = 30
num_cell_types = 5

random_seed = 0

################################################################################
#   Analysis
################################################################################

#-------------------------------------------------------------------------------
#   Preprocess and gather fluorescence data + manual cell-type labels
#-------------------------------------------------------------------------------

sample_info = pd.read_excel(sample_info_path)[:4]
sample_ids = 'Br' + sample_info['BrNumbr'].astype(int).astype(str) + \
    '_' + pd.Series([x.split('_')[1] for x in sample_info['Br_Region']]) + \
    '_IF'

#   Loop through labels and flourescence intensity tables for each sample and
#   combine into a single data frame
df = pd.DataFrame()
for i in range(len(sample_ids)):
    this_df_path = str(df_path).format(sample_ids[i])
    this_manual_label_path = str(manual_label_path).format(sample_ids[i])
    
    this_df = pd.read_csv(this_df_path, index_col = 'id')
    this_manual_labels = pd.read_csv(this_manual_label_path, index_col = 'id')
    
    this_df['label'] = this_manual_labels['label']
    df = pd.concat([df, this_df.dropna()])

#   Verify we have the correct amount of cells
assert(df.shape[0] == len(sample_ids) * expected_num_labels * num_cell_types)

#   Define the inputs (features we want the model to access) and outputs to the
#   model
x = df.loc[:, ['gfap', 'neun', 'olig2', 'tmem119']]
y = df['label']

#   Split data into training and test sets (80%: 20%), evenly stratified across
#   classes
x_train, x_test, y_train, y_test = train_test_split(
    x, y, test_size = 0.2, random_state = random_seed, stratify = y
)

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

#   TODO: once training/test split is stratified by sample, take a random subset
#   of the unlabelled examples to also stratify that by sample

#   Read in all the original unfiltered cells for all samples
this_df = pd.DataFrame()
for sample_id in sample_ids:
    this_df_small = pd.read_csv(str(df_path).format(sample_id), index_col = 'id')
    this_df_small = this_df_small.loc[:, ['gfap', 'neun', 'olig2', 'tmem119']]
    this_df = pd.concat([this_df, this_df_small])

#   Merge the DataFrame of all cells with that of just the test cells. This is
#   just a trick to categorize which of the unfiltered cells belong to the test
#   set vs. were never manually annotated
this_df = this_df.merge(
    x_test, how='left', indicator=True,
    on = ['gfap', 'neun', 'olig2', 'tmem119']
)

#   Use more informative value names
this_df['_merge'].cat.remove_unused_categories(inplace=True)
this_df['_merge'].cat.rename_categories(
    {'left_only': 'unlabeled', 'both': 'test'},
    inplace=True
)

#   Calculate confidences (maximal "probability" across classes for each
#   example) by predicted cell type and by group (test set or unlabelled)
conf_df = pd.DataFrame(
    {
        'label': grid.best_estimator_.predict(this_df.drop('_merge', axis = 1)),
        'confidence': np.max(
            grid.best_estimator_.predict_proba(
                this_df.drop('_merge', axis = 1)
            ),
            axis = 1
        ),
        'group': this_df['_merge']
    }
)

#   While not necessarily perfectly representative of the robustness of the
#   DecisionTreeClassifier, the LogisticRegression model is notably more
#   confident on test-set examples than unlabelled examples (for all cell
#   types).
print('Average model confidences on unlabeled examples vs. labelled examples in the test set:')
conf_df.groupby(['group', 'label'])['confidence'].mean()
