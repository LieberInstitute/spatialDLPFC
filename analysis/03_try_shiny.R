## Automatically style the code in this script:
styler::style_file(
    here::here("analysis", "03_try_shiny.R"),
    transformers = biocthis::bioc_style()
)

## This script requires local bioc-devel installation


## utils
library("here")
library("sessioninfo")

## reading the data
library("SpatialExperiment")

## vis
library("spatialLIBD")

## Make sure you manually download the data from JHPCE to your laptop
## Load SPE data
load(here::here("rdata", "spe", "spe.Rdata"), verbose = TRUE)

## Add random layers cluster
spe$layer_guess_reordered_short <-
    sample(c(paste0("L", 1:6), "WM"), size = ncol(spe), replace = TRUE)

## Try to run shiny
run_app(spe,
    spe_discrete_vars = "overlaps_tissue"
)
# Error in check_image_path(image_path, spe) :
#   all(file.exists(file.path(image_path, unique(spe$sample_id),  .... is not TRUE

## Current R info:
# SpatialExperiment      * 1.1.432    2021-02-17 [1] Github (drighelli/SpatialExperiment@8407ee8)
# spatialLIBD            * 1.3.5      2021-02-17 [1] Github (LieberInstitute/spatialLIBD@3fbc875)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
