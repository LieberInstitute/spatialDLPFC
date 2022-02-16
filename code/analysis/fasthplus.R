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
load(file = here::here("processed-data","rdata","spe","spe_final.Rdata"),verbose = TRUE)

spe <- cluster_import(
  spe,
  cluster_dir = here::here("processed-data", "rdata", "spe", "clustering_results","bayesSpace_k"),
  prefix = ""
)

#hpb estimate. t = pre-bootstrap sample size, D = reduced dimensions matrix, L = cluster labels, r = number of bootstrap iterations
set.seed(20220216)
fasthplus <- hpb(D= reducedDims(spe)$HARMONY,L=colData(spe)[[paste0("bayesSpace_harmony_",k)]],t=1446,r=30)
results <- data.frame (k=k, fasthplus=fasthplus)
write.csv(results,file = here::here("processed-data","rdata","spe","fasthplus_resuts.csv"), append = TRUE)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()





