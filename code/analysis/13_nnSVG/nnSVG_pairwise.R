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

# task ID
job <- as.numeric(Sys.getenv("SGE_TASK_ID"))
# Create vector of samples for nnSVG on whole tissue
sample_ids <- unique(spe$sample_id)

#create list of pairwise combinations
x <- c(1:16)
y <-c(1:16)
pairs <- expand.grid(x = x, y = y)
#256 pairs

# Run nnSVG once per sample whole tissue and store lists of top SVGs
res_list <- as.list(rep(NA, length(sample_ids)))
names(res_list) <- sample_ids

for (s in seq_along(sample_ids)) {
  
  # select sample_id
  ix <- colData(spe)$sample_id == sample_ids[s]
  spe_sub <- spe[, ix]
  
  # subset for pair of BayesSpace clusters
  cx <- colData(spe)$bayesSpace_harmony_16 %in% c(pairs[job,1],pairs[job,2])
  spe_sub <- spe[, cx]
  
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
dir_outputs <- here("processed-data", "rdata","spe", "13_nnSVG", "pairwise")

# save whole tissue nnSVG results
fn_out <- file.path(dir_outputs, paste0(pairs[job,1],"_",pairs[job,2]))
saveRDS(res_list, paste0(fn_out, ".rds"))
save(res_list, file = paste0(fn_out, ".RData"))


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()