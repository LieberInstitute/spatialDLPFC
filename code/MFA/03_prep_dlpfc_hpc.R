library(SpatialExperiment)
library(sessioninfo)
library(here)
library(tidyverse)

spe_path = here('processed-data', 'MFA', 'pseudobulk_rds', 'visium_pb.rds')
sce_path = here('processed-data', 'MFA', 'pseudobulk_rds', 'sn_fine_pb.rds')
out_path = here('processed-data', 'MFA', 'combined_rds', 'DLPFC.rds')

dir.create(dirname(out_path), showWarnings = FALSE)

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

    return(sce)
}

#   Harmonize and combine (cbind) two SingleCellExperiment objects that were
#   previously pseudobulked by donor and cluster, but may have different column
#   names for these colData variables. Return a list of inputs for MOFA
combine_sce = function(
        sce_path1, sce_path2, donor_var1, donor_var2, cluster_var1, cluster_var2,
        shared_cols
    ) {
    #   Load objects
    sce1 = readRDS(sce_path1)
    sce2 = readRDS(sce_path2)

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
mofa_list = combine_sce(
    sce_path1 = sce_path,
    sce_path2 = spe_path,
    donor_var1 = 'BrNum',
    donor_var2 = 'subject',
    cluster_var1 = 'cellType_layer',
    cluster_var2 = 'BayesSpace_harmony_09',
    shared_cols = c('age', 'sex')
)
saveRDS(mofa_list, out_path)

session_info()
