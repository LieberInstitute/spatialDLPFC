
library("SpatialExperiment")
library("spatialLIBD")
library("nnSVG")
library("here")
library("sessioninfo")

sample = "Br2720_ant"
cluster_k = "bayesSpace_harmony_16"
domain_x = 1
domain_y = 2

message("Running sample: ", sample, ", ", cluster_k, ": ", domain_x, " vs. ", domain_y)

message(Sys.time(), " - Loading spe")
load(here("processed-data", "rdata","spe", "01_build_spe", "spe_filtered_final_with_clusters.Rdata"))

stopifnot(sample %in% spe$sample_id)
stopifnot(domain_x %in% spe[[cluster_k]])
stopifnot(domain_y %in% spe[[cluster_k]])

message(Sys.time(), " - Subset spe")

spe <- spe[, spe$sample_id == sample & spe[[cluster_k]] %in% c(domain_x, domain_y)]
message("ncol: ", ncol(spe))

## stop if < 60 spots 
stopifnot(ncol(spe) > 60)

## filter out mito gnees
spe <- filter_genes(spe)

# re-calculate logcounts after filtering (?) - check this
