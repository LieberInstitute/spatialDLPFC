## Automatically style the code in this script:
styler::style_file(
    here::here("analysis", "xx_template_script.R"),
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

## Or load the filtered SPE data (whichever is most appropriate)
load(here::here("rdata", "spe", "spe.Rdata"), verbose = TRUE)


## new code


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
