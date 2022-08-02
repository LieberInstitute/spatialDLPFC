# Code adapted from Anthony Ramnauth, July 11 2022
################################################################

suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(spatialLIBD)
  library(here)
  library(scater)
  library(scran)
  library(nnSVG)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(sessioninfo)
})

#spe <- readRDS(here::here("processed-data", "harmony_processed_spe", "harmony_spe.rds"))
load("/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/01_build_spe/spe_filtered_final_with_clusters.Rdata")

#create new cluster labels to combine pairs
spe$merged_cluster <- spe$bayesSpace_harmony_16
spe$merged_cluster[spe$bayesSpace_harmony_16 %in% c(5,9)] <- "A" 
spe$merged_cluster[spe$bayesSpace_harmony_16 %in% c(4,16)] <- "B"
spe$merged_cluster[spe$bayesSpace_harmony_16 %in% c(7,13)] <- "C" 
spe$merged_cluster <- as.factor(spe$merged_cluster)

# Create vector of samples for nnSVG on whole tissue
sample_ids <- unique(spe$sample_id)

# Run nnSVG once per sample whole tissue and store lists of top SVGs
res_list <- as.list(rep(NA, length(sample_ids)))
names(res_list) <- sample_ids

for (s in seq_along(sample_ids)) {
  
  # select sample_id
  ix <- colData(spe)$sample_id == sample_ids[s]
  spe_sub <- spe[, ix]
  
  #if less than 60 spots in spe_sub, next 
  if(ncol(spe_sub)<65){
    next
  }
  
  # run nnSVG filtering for mitochondrial gene and low-expressed genes
  spe_sub <- filter_genes(spe_sub)
  
  # re-calculate logcounts after filtering
  spe_sub <- logNormCounts(spe_sub)
  
  # create model matrix for BayesSpace clusters covariates
  X <- model.matrix(~ colData(spe_sub)$merged_cluster)
  dim(X)
  head(X)
  stopifnot(nrow(X) == ncol(spe_sub))
  
  # run nnSVG
  set.seed(12345)
  message("running nnSVG")
  message(Sys.time())
  spe_sub <- nnSVG(spe_sub, X = X, n_threads = 8)
  message(Sys.time())
  # store whole tissue results
  res_list[[s]] <- rowData(spe_sub)
  message("finished adding results to results object")
}

saveRDS(res_list, file = paste0("/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/13_nnSVG/pairwise/res_list_merged_pairs.rds"))
save(res_list, file = paste0("/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/13_nnSVG/pairwise/res_list_merged_pairs.Rdata"))
message("results saved")


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()