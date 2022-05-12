library("here")
library("sessioninfo")
library("spatialLIBD")
library("RColorBrewer")
library(ggplot2)

#load spe object and clusters
load(file = here::here("processed-data","rdata","spe","01_build_spe","spe_filtered_final.Rdata"),verbose = TRUE)
spe <- cluster_import(
  spe,
  cluster_dir = here::here("processed-data", "rdata", "spe", "clustering_results"), 
  prefix = ""
)

sample_order = unique(spe$sample_id)

p_k2 <-vis_grid_clus(
  spe = spe,
  clustervar = "bayesSpace_harmony_2",
  sort_clust = FALSE,
  colors = setNames(Polychrome::palette36.colors(2), seq_len(2)),
  spatial = FALSE,
  point_size = 2.5,
  height = 24,
  width = 90,
  return_plots = TRUE
)
p_k2 <- lapply(sample_order, function(sampleid){
  p <- p_k2[[sampleid]]
  p + theme(legend.position = "none")
})
names(p_k2) <- sample_order

pdf(file = here::here("plots", "03_BayesSpace", "vis_grid_clus_sfigu_BayesSpace_k2.pdf"), height = 24, width = 36)
print(cowplot::plot_grid(plotlist = p_k2))
dev.off()


