library(SpatialExperiment)
library(sessioninfo)
library(MOFAcellulaR)
library(MOFA2)
library(here)
library(tidyverse)

dataset = 'combined'
num_factors = 5
specific_factor = 'Factor1'

model_path = here(
    'processed-data', 'MFA', 'models',
    sprintf('%s_%d.rds', dataset, num_factors)
)
sce_path = here(
    'processed-data', 'MFA', 'combined_rds', sprintf('%s.rds', dataset)
)

model = load_model(model_path)

sce = readRDS(sce_path)$sce
