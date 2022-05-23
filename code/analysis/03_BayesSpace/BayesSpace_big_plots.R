library("here")
library("sessioninfo")
library("spatialLIBD")
library("RColorBrewer")
library(ggplot2)

k <- as.numeric(Sys.getenv("SGE_TASK_ID"))

#load spe object and clusters
load(file = here::here("processed-data","rdata","spe","01_build_spe","spe_filtered_final_with_clusters.Rdata"),verbose = TRUE)

mycolors <- Polychrome::palette36.colors(k)
names(mycolors) <- sort(unique(colData(spe)[[paste0("bayesSpace_harmony_",k)]]))
bayesSpace_name <- paste0("bayesSpace_harmony_", k)
sample_ids <- unique(colData(spe)$sample_id)

pdf(file = here::here("plots","03_BayesSpace",paste0("polychrome_vis_clus_bayesSpace_harmony_",k,".pdf")))
for (i in seq_along(sample_ids)){
  my_plot <- vis_clus(
    spe = spe,
    clustervar = bayesSpace_name,
    sampleid = sample_ids[i],
    colors = mycolors
  )
  print(my_plot)
}
dev.off()