library(SpatialExperiment)
library(sessioninfo)
library(here)
library(tidyverse)

spe_path = here('processed-data', 'MFA', 'pseudobulk_rds', 'visium_pb.rds')
sce_path = here('processed-data', 'MFA', 'pseudobulk_rds', 'sn_fine_pb.rds')
out_path = here('processed-data', 'MFA', 'combined_rds', 'DLPFC.rds')

dir.create(dirname(out_path), showWarnings = FALSE)

#   Load objects
spe = readRDS(spe_path)
sce = readRDS(sce_path)

#   Ensure genes match
shared_genes = intersect(rownames(spe), rownames(sce))
spe = spe[shared_genes, ]
sce = sce[shared_genes, ]
rowRanges(sce) = rowRanges(spe)

#   Simplify SPE object in preparation for combining with SCE
temp_colnames = colnames(spe)
colData(spe) = colData(spe) |>
    as_tibble() |>
    dplyr::rename(donor = subject, cluster = BayesSpace_harmony_09) |>
    select(donor, cluster, age, sex, ncells) |>
    DataFrame()
colnames(spe) = temp_colnames
reducedDims(spe) = list()
metadata(spe) = list()
int_colData(spe)$spatialCoords = NULL

#   Simplify SCE object in preparation for combining with SPE
temp_colnames = colnames(sce)
colData(sce) = colData(sce) |>
    as_tibble() |>
    dplyr::rename(donor = BrNum, cluster = cellType_layer) |>
    select(donor, cluster, age, sex, ncells) |>
    DataFrame()
colnames(sce) = temp_colnames
reducedDims(sce) = list()
metadata(sce) = list()

#   Gather some donor-level info for later
pd = spe[, match(unique(spe$donor), spe$donor)] |>
    colData() |>
    as_tibble() |>
    select(donor, age, sex)

#   Merge Visium and snRNA-seq objects and save
sce = cbind(spe, sce)
saveRDS(sce, out_path)

session_info()
