################################################################
# spatial_DG_lifespan project
# nnSVG per Capture Area, Average ranks, & BayesSpace covariates
# Anthony Ramnauth, July 11 2022
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


# Create vector of samples for nnSVG on whole tissue
sample_ids <- unique(spe$sample_id)

# Run nnSVG once per sample whole tissue and store lists of top SVGs

res_list <- as.list(rep(NA, length(sample_ids)))
names(res_list) <- sample_ids

for (s in seq_along(sample_ids)) {
  
  # select sample_id
  ix <- colData(spe)$sample_id == sample_ids[s]
  spe_sub <- spe[, ix]
  
  # run nnSVG filtering for mitochondrial gene and low-expressed genes
  spe_sub <- filter_genes(spe_sub)
  
  # re-calculate logcounts after filtering
  spe_sub <- logNormCounts(spe_sub)
  
  # run nnSVG
  set.seed(12345)
  spe_sub <- nnSVG(spe_sub, n_threads = 8)
  
  # store whole tissue results
  res_list[[s]] <- rowData(spe_sub)
}

# directory to save whole tissue results
dir_outputs <- here("processed-data", "nnSVG", "whole_tissue")

# save whole tissue nnSVG results
fn_out <- file.path(dir_outputs, "DG_nnSVG_results")
saveRDS(res_list, paste0(fn_out, ".rds"))
save(res_list, file = paste0(fn_out, ".RData"))


# Run nnSVG once per sample whole tissue with BayesSpace covariates

bayes_res_list <- as.list(rep(NA, length(sample_ids)))
names(bayes_res_list) <- sample_ids

for (s in seq_along(sample_ids)) {
  
  # select sample_id
  ix <- colData(spe)$sample_id == sample_ids[s]
  spe_subS <- spe[, ix]
  
  # run nnSVG filtering for mitochondrial gene and low-expressed genes
  spe_subS <- filter_genes(spe_subS)
  
  # re-calculate logcounts after filtering
  spe_subS <- logNormCounts(spe_subS)
  
  # create model matrix for BayesSpace clusters covariates
  X <- model.matrix(~ colData(spe_subS)$bayesSpace_harmony_8)
  dim(X)
  head(X)
  stopifnot(nrow(X) == ncol(spe_subS))
  
  # run nnSVG
  set.seed(12345)
  spe_subS <- nnSVG(spe_subS, X = X, n_threads = 8)
  
  # store whole tissue results
  bayes_res_list[[s]] <- rowData(spe_subS)
}


# directory to save bayesspace informed results
dir_outputs <- here("processed-data", "nnSVG", "BayesSpace")

# save bayesspace nnSVG results
fn_out <- file.path(dir_outputs, "DG_BayesSpace_nnSVG_results")
saveRDS(bayes_res_list, paste0(fn_out, ".rds"))
save(bayes_res_list, file = paste0(fn_out, ".RData"))


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()