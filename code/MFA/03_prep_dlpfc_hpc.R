library(SpatialExperiment)
library(sessioninfo)
library(here)
library(tidyverse)

dlpfc_spe_path = here('processed-data', 'MFA', 'pseudobulk_rds', 'visium_pb.rds')
dlpfc_sce_path = here(
    'processed-data', 'MFA', 'pseudobulk_rds', 'sn_fine_pb.rds'
)
hpc_spe_path = here(
    'processed-data', 'MFA', 'pseudobulk_rds', 'hpc_visium_pb.rds'
)
hpc_sce_path = here(
    'processed-data', 'MFA', 'pseudobulk_rds', 'hpc_sn_pb.rds'
)
combined_out_path = here(
    'processed-data', 'MFA', 'combined_rds', 'combined.rds'
)

dir.create(dirname(combined_out_path), showWarnings = FALSE)

################################################################################
#   Functions
################################################################################

#   Reduce a SingleCellExperiment object into a smaller version that may be
#   combined with a similar reduced object
simplify_sce = function(sce, donor_var, cluster_var, shared_cols) {
    temp_colnames = colnames(sce)
    colData(sce) = colData(sce) |>
        as_tibble() |>
        dplyr::rename(donor = {{ donor_var }}, cluster = {{ cluster_var }}) |>
        select(donor, cluster, ncells, all_of(shared_cols)) |>
        DataFrame()
    colnames(sce) = temp_colnames
    reducedDims(sce) = list()
    metadata(sce) = list()
    int_colData(sce)$spatialCoords = NULL
    int_metadata(sce)$imgData = NULL

    return(sce)
}

#   Harmonize and combine (cbind) two SingleCellExperiment objects that were
#   previously pseudobulked by donor and cluster, but may have different column
#   names for these colData variables. Return a list of inputs for MOFA
combine_sce = function(
        sce1, sce2, donor_var1, donor_var2, cluster_var1, cluster_var2,
        shared_cols
    ) {
    #   Ensure genes match
    shared_genes = intersect(rownames(sce1), rownames(sce2))
    sce1 = sce1[shared_genes, ]
    sce2 = sce2[shared_genes, ]
    rowRanges(sce1) = rowRanges(sce2)

    sce1 = simplify_sce(sce1, donor_var1, cluster_var1, shared_cols)
    sce2 = simplify_sce(sce2, donor_var2, cluster_var2, shared_cols)

    #   Gather some donor-level info for later
    pd = sce1[, match(unique(sce1$donor), sce1$donor)] |>
        colData() |>
        as_tibble() |>
        select(donor, all_of(shared_cols))

    return(list(sce = cbind(sce1, sce2), pd = pd))
}

################################################################################
#   Main
################################################################################

#   Combine spatialDLPFC Visium and fine snRNA-seq objects
dlpfc_mofa_list = combine_sce(
    sce1 = readRDS(dlpfc_sce_path),
    sce2 = readRDS(dlpfc_spe_path),
    donor_var1 = 'BrNum',
    donor_var2 = 'subject',
    cluster_var1 = 'cellType_layer',
    cluster_var2 = 'BayesSpace_harmony_09',
    shared_cols = c('age', 'sex')
)

#   Combine HPC Visium and snRNA-seq objects
hpc_spe = readRDS(hpc_spe_path)
hpc_spe$cluster = NULL
hpc_mofa_list = combine_sce(
    sce1 = readRDS(hpc_sce_path),
    sce2 = hpc_spe,
    donor_var1 = 'brnum',
    donor_var2 = 'brnum',
    cluster_var1 = 'cell.type',
    cluster_var2 = 'domain',
    shared_cols = c()
)

#   Combine DLPFC and HPC
stopifnot(setequal(dlpfc_mofa_list$sce$donor, hpc_mofa_list$sce$donor))
combined_mofa_list = combine_sce(
    sce1 = dlpfc_mofa_list$sce,
    sce2 = hpc_mofa_list$sce,
    donor_var1 = 'donor',
    donor_var2 = 'donor',
    cluster_var1 = 'cluster',
    cluster_var2 = 'cluster',
    shared_cols = c()
)
combined_mofa_list$pd = dlpfc_mofa_list$pd
saveRDS(combined_mofa_list, combined_out_path)

session_info()
