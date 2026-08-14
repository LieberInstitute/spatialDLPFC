library("here")
library("sessioninfo")
library("SpatialExperiment")
library("spatialLIBD")
library("ggplot2")
library("BayesSpace")

load(file = here::here("processed-data", "rdata", "spe", "01_build_spe", "spe_filtered_final.Rdata"), verbose = TRUE)

set.seed(030122)


## Set the array coords in a way that it'll work with 
## different versions of BayesSpace
spe$array_row <- spe$row
spe$array_col <- spe$col
spe$pxl_col_in_fullres <- spatialCoords(spe)[, "pxl_col_in_fullres"]
spe$pxl_row_in_fullres <- spatialCoords(spe)[, "pxl_row_in_fullres"]

## Check visually that there's no overlap between the samples
pdf(file = here::here("plots", "03_BayesSpace", paste0("BayesSpace_version_", packageVersion("BayesSpace"), "_offset_check.pdf")))
clusterPlot(spe, "subject", color = NA) + # make sure no overlap between samples
    labs(fill = "Subject", title = "Offset check")
dev.off()

### BayesSpace on Batch Corrected
spe <- spatialCluster(spe, use.dimred = "HARMONY", q = 9, nrep = 10000)

spe$bayesSpace_temp <- spe$spatial.cluster
bayesSpace_name <- paste0("BayesSpace_version_", packageVersion("BayesSpace"))
colnames(colData(spe))[ncol(colData(spe))] <- bayesSpace_name

cluster_export(
    spe,
    bayesSpace_name,
    cluster_dir = here::here("processed-data", "rdata", "spe", "clustering_results", "BayesSpace_versions")
)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
