library("SpatialExperiment")
library("iSEE")
library("shiny")
library("scuttle")

## Load the pseudobulked object sce_pseudo
sce_pseudo <- readRDS("sce_pseudo_BayesSpace_k09.rds")

## Make unique gene names
rownames(sce_pseudo) <-
    scuttle::uniquifyFeatureNames(rowData(sce_pseudo)$gene_id, rowData(sce_pseudo)$gene_name)

## Don't run this on app.R since we don't want to run this every single time
# lobstr::obj_size(sce_pseudo)
# 56.41 MB

source("initial.R", echo = TRUE, max.deparse.length = 500)

## From https://github.com/LieberInstitute/10xPilot_snRNAseq-human/blob/810b47364af4c8afe426bd2a6b559bd6a9f1cc98/shiny_apps/tran2021_AMY/app.R#L10-L14
## Related to https://github.com/iSEE/iSEE/issues/568
colData(sce_pseudo) <- cbind(
    colData(sce_pseudo)[, !colnames(colData(sce_pseudo)) %in% c("subject", "BayesSpace")],
    colData(sce_pseudo)[, c("BayesSpace", "subject")]
)

sce_pseudo$subject <- as.factor(sce_pseudo$subject)

sce_pseudo <-
    registerAppOptions(sce_pseudo, color.maxlevels = length(colData(sce_pseudo)$BayesSpace_colors))

iSEE(
    sce_pseudo,
    appTitle = "spatialDLPFC, Visium, Sp09, pseudo-bulked",
    initial = initial,
    colormap = ExperimentColorMap(colData = list(
        BayesSpace = function(n) {
            return(colData(sce_pseudo)$BayesSpace_colors)
        }
    ))
)
