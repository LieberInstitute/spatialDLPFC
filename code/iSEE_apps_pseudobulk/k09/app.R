
library("SingleCellExperiment")
library("iSEE")
library("shiny")
library("paletteer")
library("scuttle")
library("SpatialExperiment")

# load("sce_for_iSEE_LS.rda", verbose = TRUE)# load the pseudobulked object spe_pseudo
spe_pseudo <- readRDS("spe_pseudobulk_bayesSpace_normalized_filtered_region_k9.RDS")

## Make unique gene names
rownames(spe_pseudo) <-
    uniquifyFeatureNames(rowData(spe_pseudo)$gene_id, rowData(spe_pseudo)$gene_name)

# stopifnot(all(unique(spe_pseudo$BayesSpace) %in% names(cell_cols.clean)))

## Don't run this on app.R since we don't want to run this every single time
# lobstr::obj_size(spe_pseudo)
# 876.33 MB

source("initial.R", print.eval = TRUE)

## From https://github.com/LieberInstitute/10xPilot_snRNAseq-human/blob/810b47364af4c8afe426bd2a6b559bd6a9f1cc98/shiny_apps/tran2021_AMY/app.R#L10-L14
## Related to https://github.com/iSEE/iSEE/issues/568
colData(spe_pseudo) <- cbind(
    colData(spe_pseudo)[, !colnames(colData(spe_pseudo)) %in% c("subject", "BayesSpace")],
    colData(spe_pseudo)[, c("BayesSpace", "subject")]
)

spe_pseudo$subject <- as.factor(spe_pseudo$subject)

spe_pseudo <- registerAppOptions(spe_pseudo, color.maxlevels = length(colData(spe_pseudo)$BayesSpace_colors))

iSEE(
    spe_pseudo,
    appTitle = "Spangler2022_pseudobulk_cluster_k09",
    initial = initial,
    colormap = ExperimentColorMap(colData = list(
        # subject = function(n) {
        #     cols <- paletteer::paletteer_d(
        #         palette = "RColorBrewer::Dark2",
        #         n = length(unique(spe_pseudo$subject))
        #     )
        #     cols <- as.vector(cols)
        #     names(cols) <- levels(spe_pseudo$subject)
        #     return(cols)
        # },
        BayesSpace = function(n) {
            return(colData(spe_pseudo)$BayesSpace_colors)
        }
    ))
)
