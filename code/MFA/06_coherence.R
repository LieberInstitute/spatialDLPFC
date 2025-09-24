library(getopt)
library(sessioninfo)
library(SpatialExperiment)
library(here)
library(edgeR)
library(jaffelab)
library(tidyverse)

# Import command-line parameters
spec <- matrix(
    c(
        c("DLPFC_clus", "HPC_clus", "filter_by_region", "clean_expression"),
        c("D", "H", "f", "c"),
        rep("1", 4),
        c("character", "character", "logical", "logical"),
        rep("Add variable description here", 4)
    ),
    ncol = 5
)
opt <- getopt(spec)

sce_path = here('processed-data', 'MFA', 'combined_rds', 'combined.rds')
out_dir = here('processed-data', 'MFA', 'temp_coherence')
covariates = 'ncells'

message("Using the following parameters:")
print(opt)

dir.create(out_dir, showWarnings = FALSE)

sce = readRDS(sce_path)$sce
sce = sce[, sce$cluster %in% c(opt$DLPFC_clus, opt$HPC_clus)]

#   Filter to sufficiently expressed genes
if (opt$filter_by_region) {
    exp_genes = (
        filterByExpr(assays(sce[, sce$cluster == opt$DLPFC_clus])$counts) &
        filterByExpr(assays(sce[, sce$cluster == opt$HPC_clus])$counts)
    )
} else {
    exp_genes = filterByExpr(assays(sce)$counts)
}
stopifnot(length(which(exp_genes)) > 0)
sce = sce[exp_genes, ]

sce_dlpfc = sce[, sce$cluster == opt$DLPFC_clus]
sce_hpc = sce[, sce$cluster == opt$HPC_clus]

#   Either clean expression (regress out effect of covariates) or directly
#   use logcounts, depending on opt$clean_expression
if (opt$clean_expression) {
    dlpfc_mod = model.matrix(
        as.formula(paste('~', paste(covariates, collapse = " + "))),
        data = as.data.frame(colData(sce_dlpfc))
    )[, 2:(1 + length(covariates)), drop = FALSE]
    hpc_mod = model.matrix(
        as.formula(paste('~', paste(covariates, collapse = " + "))),
        data = as.data.frame(colData(sce_hpc))
    )[, 2:(1 + length(covariates)), drop = FALSE]
    dlpfc_exp = cleaningY(assays(sce_dlpfc)$logcounts, dlpfc_mod, P = 1)
    hpc_exp = cleaningY(assays(sce_hpc)$logcounts, hpc_mod, P = 1)
} else {
    dlpfc_exp = assays(sce_dlpfc)$logcounts
    hpc_exp = assays(sce_hpc)$logcounts
}

#   Compute correlation of expression between regions for each donor
co_df_list = list()
for (donor in intersect(sce_dlpfc$donor, sce_hpc$donor)) {
    co_df_list[[donor]] = tibble(
        donor = donor,
        coherence = cor(
            dlpfc_exp[, sce_dlpfc$donor == donor],
            hpc_exp[, sce_hpc$donor == donor]
        )
    )
}

#   Save as one tibble, adding APOE-related variables and input options
do.call(rbind, co_df_list) |>
    mutate(
        DLPFC_cluster = opt$DLPFC_clus,
        HPC_cluster = opt$HPC_clus,
        filter_by_region = opt$filter_by_region,
        clean_expression = opt$clean_expression
    ) |>
    saveRDS(
        file.path(
            out_dir,
            sprintf(
                '%s_%s_filter%s_clean%s.rds',
                opt$DLPFC_clus, opt$HPC_clus, opt$filter_by_region,
                opt$clean_expression
            )
        )
    )

message("Memory usage:")
gc()

session_info()

## This script was made using slurmjobs version 1.3.0
## available from http://research.libd.org/slurmjobs/
