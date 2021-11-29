library("here")
library("lobstr")
library("spatialLIBD")

load(here("processed-data", "rdata", "spe", "spe_merged_final.Rdata"), verbose = TRUE)

## Check memory size
lobstr::obj_size(spe) / 1024^3
# 7.092369

## Inspect how big are some parts of the spe object
lobstr::obj_size(counts(spe)) / 1024^3
# 2.196542
lobstr::obj_size(imgData(spe)) / 1024^3
# 2.515935
lobstr::obj_size(subset(imgData(spe), image_id == "hires")) / 1024^3
# 0.8242114
lobstr::obj_size(subset(imgData(spe), image_id == "detected")) / 1024^3
# 0.8113759

## Drop the counts, otherwise the object is too big since it's too close to the
## 8 GB limit on shinyapps.io
counts(spe) <- NULL
lobstr::obj_size(spe) / 1024^3
# 4.897936

## Also drop the detected & aligned images
imgData(spe) <- subset(imgData(spe), image_id %in% c("lowres", "hires"))
lobstr::obj_size(spe) / 1024^3
# 3.278569
## Other options:
# 4.084106 (dropping only hires)
# 2.463861 (keeping only lowres)

save(spe, file = here("code", "deploy_app", "spe_merged_final_nocounts.Rdata"))
