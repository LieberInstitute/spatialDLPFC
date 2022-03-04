#set up arrayjob to run k=2 to k = 15
#don't use spatial preprocess. in order to do this you have to reset metadata
# increase nrep for spatialCluster??

library("here")
library("sessioninfo")
library("SpatialExperiment")
library("spatialLIBD")
library(BayesSpace)
library("RColorBrewer")

load(file = here::here("processed-data","rdata","spe","01_build_spe","spe_filtered_final.Rdata"),verbose = TRUE)

k <- as.numeric(Sys.getenv("SGE_TASK_ID"))

set.seed(030122)

##do offset so we can run BayesSpace
auto_offset_row <- as.numeric(factor(unique(spe$sample_id))) * 100
names(auto_offset_row) <-unique(spe$sample_id)
spe$row <- colData(spe)$array_row + auto_offset_row[spe$sample_id]
spe$col <- colData(spe)$array_col

save(spe, file = here::here("processed-data","rdata","spe","01_build_spe","spe_filtered_final.Rdata"),verbose = TRUE)

pdf(file=here::here("plots","03_BayesSpace", "BayesSpace_offset_check.pdf"))
clusterPlot(spe, "subject", color = NA) + #make sure no overlap between samples
  labs(fill = "Subject", title = "Offset check")
dev.off()

###BayesSpace on Batch Corrected
spe = spatialCluster(spe, use.dimred = "HARMONY", q = k, nrep = 10000) 

spe$bayesSpace_temp<-spe$spatial.cluster
bayesSpace_name <- paste0("bayesSpace_harmony_", k)
colnames(colData(spe))[ncol(colData(spe))] <- bayesSpace_name

cluster_export(
  spe,
  bayesSpace_name,
  cluster_dir = here::here("processed-data", "rdata", "spe", "clustering_results" )
)

sample_ids <- unique(colData(spe)$sample_id)
mycolors <- brewer.pal(7, "Dark2")

pdf(file = here::here("plots","03_BayesSpace",paste0("vis_clus_bayesSpace_harmony_",k,".pdf")))
for (i in seq_along(sample_ids)){
  my_plot <- vis_clus(
    spe = spe,
    clustervar = bayesSpace_name,
    sampleid = sample_ids[i],
    colors =  mycolors
  )
  print(my_plot)
}
dev.off()

###BayesSpace on Non-Batch Corrected
spe = spatialCluster(spe, use.dimred = "PCA", q = k, nrep = 10000) 

spe$bayesSpace_temp<-spe$spatial.cluster
bayesSpace_name <- paste0("bayesSpace_pca_", k)
colnames(colData(spe))[ncol(colData(spe))] <- bayesSpace_name

cluster_export(
  spe,
  bayesSpace_name,
  cluster_dir = here::here("processed-data", "rdata", "spe", "clustering_results" )
)

sample_ids <- unique(colData(spe)$sample_id)
mycolors <- brewer.pal(7, "Dark2")

pdf(file = here::here("plots","03_BayesSpace",paste0("vis_clus_bayesSpace_pca_",k,".pdf")))
for (i in seq_along(sample_ids)){
  my_plot <- vis_clus(
    spe = spe,
    clustervar = bayesSpace_name,
    sampleid = sample_ids[i],
    colors =  mycolors
  )
  print(my_plot)
}
dev.off()


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
