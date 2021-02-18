## Automatically style the code in this script:
styler::style_file(
    here::here("analysis", "01a_check_scran_quick_cluster.R"),
    transformers = biocthis::bioc_style()
)

## This script requires R 4.1
# module load conda_R/devel

## utils
library("here")
library("sessioninfo")

## reading the data
library("SpatialExperiment")

## vis
library("spatialLIBD")

## Load SPE data
load(here::here("rdata", "spe", "spe.Rdata"), verbose = TRUE)

## Check cluster 89 from before
load(
    file = here::here("analysis", "clusters_nonzero.rda"),
    verbose = TRUE
)

table(spe$scran_quick_cluster == clusters)
## Clusters did change!
# FALSE  TRUE
# 41244  8755

## old cluster "89" from https://github.com/LTLA/scuttle/issues/7
## is not the same one now, despite the number of spots matching
identical(which(spe$scran_quick_cluster == "89"), which(clusters == "89")) # FALSE
stopifnot(identical(sum(spe$scran_quick_cluster == "89"), sum(clusters == "89")))

## Hm... maybe it's because the samples are in a different order!
## due to https://github.com/LieberInstitute/spatialDLPFC/blob/main/analysis/01_build_SPE.R#L72
stopifnot(sum(table(spe$scran_quick_cluster) - table(clusters)) == 0)
## Ok, that matches. Let's check more closely then.
load(
    file = here::here("analysis", "sce_nonzero.rda"),
    verbose = TRUE
)
m <-
    match(spe$key, with(colData(sce_nonzero), paste0(barcode_id, "_", sample_name)))
## Ok, the quick clusters are the same ones!
stopifnot(identical(spe$scran_quick_cluster, clusters[m]))

## Remove the older R objects since we don't need them
rm(clusters, sce_nonzero)

table(spe$scran_quick_cluster == "89")
# FALSE  TRUE
# 49727   272
spe$quick_cluster_89 <-
    factor(spe$scran_quick_cluster == "89", levels = c("TRUE", "FALSE"))

summary(as.data.frame(colData(spe)[spe$scran_quick_cluster == "89", c(
    "scran_low_lib_size",
    "scran_low_n_features",
    "scran_high_subsets_Mito_percent",
    "scran_discard"
)]))
# scran_low_lib_size scran_low_n_features scran_high_subsets_Mito_percent
# TRUE : 62          TRUE : 86            TRUE :  2
# FALSE:210          FALSE:186            FALSE:270
# scran_discard
# TRUE : 87
# FALSE:185
vis_grid_clus(
    spe = spe,
    clustervar = "quick_cluster_89",
    pdf = here::here("plots", paste0("scuttle_", "quick_cluster_89", ".pdf")),
    sort_clust = FALSE,
    colors = c("FALSE" = "grey90", "TRUE" = "orange")
)

## Now check cluster 64 since that's the one giving the error at
## https://github.com/LieberInstitute/spatialDLPFC/blob/main/analysis/01_build_SPE.R#L338
spe$quick_cluster_64 <-
    factor(spe$scran_quick_cluster == "64", levels = c("TRUE", "FALSE"))
summary(as.data.frame(colData(spe)[spe$scran_quick_cluster == "64", c(
    "scran_low_lib_size",
    "scran_low_n_features",
    "scran_high_subsets_Mito_percent",
    "scran_discard"
)]))
# scran_low_lib_size scran_low_n_features scran_high_subsets_Mito_percent
# TRUE :  0          TRUE :  0            TRUE : 56
# FALSE:542          FALSE:542            FALSE:486
# scran_discard
# TRUE : 56
# FALSE:486
vis_grid_clus(
    spe = spe,
    clustervar = "quick_cluster_64",
    pdf = here::here("plots", paste0("scuttle_", "quick_cluster_64", ".pdf")),
    sort_clust = FALSE,
    colors = c("FALSE" = "grey90", "TRUE" = "orange")
)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
