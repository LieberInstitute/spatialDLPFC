###########################
#fasthplus, need to use hpb()
##########################

#install_github(repo="ntdyjack/fasthplus", ref = "main")
library(fasthplus)
library(SpatialExperiment)
library(here)
library("sessioninfo")
library(spatialLIBD)

k <- as.numeric(Sys.getenv("SGE_TASK_ID"))  

#load spe object
load(file = here::here("processed-data","rdata","spe","01_build_spe","spe_filtered_final.Rdata"),verbose = TRUE)

spe <- cluster_import(
  spe,
  cluster_dir = here::here("processed-data", "rdata", "spe", "clustering_results"), 
  prefix = ""
)

##remove white matter
dim(spe)
# [1]  28916 113927
spe <- spe[,-which(colData(spe)$bayesSpace_harmony_2 == 1)]
dim(spe)
# [1] 28916 99574

#hpb estimate. t = pre-bootstrap sample size, D = reduced dimensions matrix, L = cluster labels, r = number of bootstrap iterations
dim(reducedDims(spe)$HARMONY)
#[1] 99574    50

set.seed(20220216)
fasthplus <- hpb(D= reducedDims(spe)$HARMONY,L=colData(spe)[[paste0("bayesSpace_harmony_",k)]],t=200,r=30) # t= 99574*0.01
results <- data.frame (k=k, fasthplus=fasthplus)
write.table(results,file = here::here("processed-data","rdata","spe","06_fasthplus","fasthplus_results_no_WM.csv"), append = TRUE)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

find_t <- function(L, proportion = 0.05) {
  initial_t <- floor(length(L) * proportion)
  smallest_cluster_size <- min(table(L))
  n_labels <- length(unique(L))
  ifelse(smallest_cluster_size > (initial_t / n_labels), initial_t, smallest_cluster_size * n_labels)
}

initial_t <- find_t(L=colData(spe)[[paste0("bayesSpace_harmony_",k)]],proportion = 0.01)
# sort(table(colData(spe)[[paste0("bayesSpace_harmony_",k)]]))

cluster_prop <- table(colData(spe)[[paste0("bayesSpace_harmony_",k)]]) / ncol(spe)
bad_clusters <- which(cluster_prop < 0.01 / k)
if(length(bad_clusters) > 0) {
  message("For k: ", k, " we are dropping small clusters: ", paste(names(bad_clusters), collapse = ", "))
  spe <- spe[, !colData(spe)[[paste0("bayesSpace_harmony_",k)]] %in% as.integer(names(bad_clusters))]
  updated_t <- find_t(colData(spe)[[paste0("bayesSpace_harmony_", k)]], 0.01)
  message("initial t: ", initial_t, "; updated t: ", updated_t)
}
