library("SingleCellExperiment")
library("spatialLIBD")
library("tidyverse")
library("scMerge")
library("here")
library("sessioninfo")

## Load data
genes_of_intrest <- c("SNAP25","FYN")

pseudobulk_fn <- list.files(here("processed-data", "rdata", "spe", "14_spatial_registration_PEC"), pattern = "pseudobulk", full.names = TRUE)
names(pseudobulk_fn) <- gsub("pseudobulk_|.rds","",basename(pseudobulk_fn))

all_pseudobulk <- map(pseudobulk_fn, readRDS)
map(all_pseudobulk, dim)

sce_pseudo <- readRDS(file = here(
  "processed-data", "rdata", "spe", "14_spatial_registration_PEC",
  paste0("pseudobulk_", dataset, ".rds")
))

