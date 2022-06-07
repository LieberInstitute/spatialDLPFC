library("SpatialExperiment")
library("lobstr")
library("here")
library("sessioninfo")

## Load the spe object
load("spe_filtered_final_with_clusters.Rdata", verbose = TRUE)

## Check the size in GB
lobstr::obj_size(spe) / 1024^3
# 6.442109 B

table(imgData(spe)$image_id)
# aligned detected    hires   lowres
#      30       30       30       30

imgData(spe) <- imgData(spe)[imgData(spe)$image_id == "lowres", ]
## Check the size in GB
lobstr::obj_size(spe) / 1024^3
# 4.231281 B

## Drop the counts for now
counts(spe) <- NULL
## Check the size in GB
lobstr::obj_size(spe) / 1024^3
# 2.238795 B

save(spe, file = here::here("code", "deploy_app", "spe_subset.Rdata"))

lobstr::obj_size(spe_pseudo) / 1024^3
# 0.05300035 B

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
