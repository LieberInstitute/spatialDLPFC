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

## Load SPE raw data
load(here::here("rdata", "spe", "spe_raw.Rdata"), verbose = TRUE)

## Filter down to spots in tissue
spe <- spe_raw[, which(inTissue(spe_raw))]

## Remove spots without counts
spe <- spe[, -which(colSums(counts(spe)) == 0)]


## Check cluster 89 from before
load(file = here::here("analysis", "clusters_nonzero.rda"), verbose = TRUE)
table(clusters == "89")
# FALSE  TRUE
# 49727   272
stopifnot(length(clusters) == ncol(spe))
spe$quick_cluster_89 <- factor(clusters == "89", levels = c("TRUE", "FALSE"))
colSums(as.matrix(qcfilter[clusters == "89", ]))
# low_lib_size            low_n_features high_subsets_Mito_percent                   discard
#            2                         3                         0                         3
vis_grid_clus(
    spe = spe,
    clustervar = "quick_cluster_89",
    pdf = here::here("plots", paste0("scuttle_", "quick_cluster_89", ".pdf")),
    sort_clust = FALSE,
    colors = c("FALSE" = "grey90", "TRUE" = "orange")
)


## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
