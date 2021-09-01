## Automatically style the code in this script:
styler::style_file(here::here("analysis", "03_try_shiny.R"),
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

## Make it easy to copy over the data for the images in the format expected by
## spatialLIBD
dir.create(here::here("outputs", "images_spatialLIBD"), showWarnings = FALSE)
sapply(imgData(spe)$sample_id, function(x) {
    dir.create(here::here("outputs", "images_spatialLIBD", x), showWarnings = FALSE)
    new_image <-  here::here("outputs", "images_spatialLIBD", x, "tissue_lowres_image.png")
    file.copy(imgPath(spe, sample_id = x), new_image)
    new_image
})

## Add random layers cluster
set.seed(20210217)
spe$layer_guess_reordered_short <-
    sample(c(paste0("L", 1:6), "WM"), size = ncol(spe), replace = TRUE)

## Try to run shiny
run_app(
    spe,
    spe_discrete_vars = "overlaps_tissue",
    spe_continuous_vars = c(
        "sum_umi",
        "sum_gene",
        "expr_chrM",
        "expr_chrM_ratio"
    ),
    #image_path = "/Users/lcollado/Dropbox/Code/spatialDLPFC/outputs/images_spatialLIBD"
    image_path = here::here("outputs", "images_spatialLIBD")
)

## Current R info:
# SpatialExperiment      * 1.1.432    2021-02-17 [1] Github (drighelli/SpatialExperiment@8407ee8)
# spatialLIBD            * 1.3.5      2021-02-17 [1] Github (LieberInstitute/spatialLIBD@3fbc875)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
