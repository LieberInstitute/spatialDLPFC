library("here")
library("sessioninfo")
library("spatialLIBD")
library("RColorBrewer")
library(ggplot2)

k <- as.numeric(Sys.getenv("SGE_TASK_ID"))

#load spe object and clusters
load(file = here::here("processed-data","rdata","spe","01_build_spe","spe_filtered_final.Rdata"),verbose = TRUE)
spe <- cluster_import(
  spe,
  cluster_dir = here::here("processed-data", "rdata", "spe", "clustering_results"), 
  prefix = ""
)


vis_grid_clus(
  spe = spe,
  clustervar = paste0("bayesSpace_harmony_",k),
  pdf = here("plots", "03_BayesSpace", paste0("vis_grid_clus_BayesSpace_k",k,".pdf")),
  sort_clust = FALSE,
  colors = setNames(Polychrome::palette36.colors(k), 1:k),
  spatial = FALSE,
  point_size = 2,
  height = 24,
  width = 90
)
