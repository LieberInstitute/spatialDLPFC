library(SpatialExperiment)
library(here)
library(pheatmap)
library(spatialLIBD)

#load spe object
load(file = here::here("processed-data","rdata","spe","01_build_spe","spe_filtered_final.Rdata"),verbose = TRUE)

#load clusters
spe <- cluster_import(
  spe,
  cluster_dir = here::here("processed-data", "rdata", "spe", "clustering_results"), 
  prefix = ""
)

#make table of bayesSpace k = 9 vs k = 28
my.table <- with(colData(spe), table(bayesSpace_harmony_9, bayesSpace_harmony_28))