library("here")
library("sessioninfo")
library("spatialLIBD")
library("RColorBrewer")
library(ggplot2)

k <- as.numeric(Sys.getenv("SGE_TASK_ID"))

#load spe object and clusters
load(file = here::here("processed-data","rdata","spe","01_build_spe","spe_filtered_final_with_clusters.Rdata"),verbose = TRUE)

mycolors <- brewer.pal(8, "Dark2")
sample_order = unique(spe$sample_id)

pk <- vis_grid_clus(
  spe = spe,
  clustervar = paste0("bayesSpace_harmony_",k),
  sort_clust = FALSE,
  colors = mycolors,
  spatial = FALSE,
  point_size = 1.5,
  #height = 24,
  #width = 90,
  return_plots = TRUE
)

pk <- lapply(sample_order, function(sampleid){
  p <- pk[[sampleid]]
  p + theme(legend.position = "none")
})
names(pk) <- sample_order

pdf(file = here::here("plots", "03_BayesSpace", paste0("test_vis_grid_clus_sfigu_BayesSpace_k",k,".pdf")), height = 5*8, width = 6*8)
print(cowplot::plot_grid(plotlist = pk))
dev.off()

# bayesSpace_name <- paste0("bayesSpace_harmony_", k)
# sample_ids <- unique(colData(spe)$sample_id)
# 
# pdf(file = here::here("plots","03_BayesSpace",paste0("test_vis_clus_bayesSpace_harmony_",k,".pdf")))
# for (i in seq_along(sample_ids)){
#   my_plot <- vis_clus(
#     spe = spe,
#     clustervar = bayesSpace_name,
#     sampleid = sample_ids[i],
#     colors = mycolors
#   )
#   print(my_plot)
# }
# dev.off()
