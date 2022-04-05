library(SpatialExperiment)
library(spatialLIBD)
library(here)
library(edgeR)
library(scuttle)

k <- as.numeric(Sys.getenv("SGE_TASK_ID"))

#load spe object and clusters
load(file = here::here("processed-data","rdata","spe","01_build_spe","spe_filtered_final.Rdata"),verbose = TRUE)
spe <- cluster_import(
  spe,
  cluster_dir = here::here("processed-data", "rdata", "spe", "clustering_results"), 
  prefix = ""
)

## Pseudo-bulk for our current BayesSpace cluster results
spe_pseudo <- aggregateAcrossCells(
  spe,
  DataFrame(
    BayesSpace = colData(spe)[[paste0("bayesSpace_harmony_",k)]],
    sample_id = spe$sample_id
  )
)

#log normalize the counts
spe_pseudo <- logNormCounts(spe_pseudo,size.factors = NULL)

#find a good expression cutoff using edgeR::filterByExpr https://rdrr.io/bioc/edgeR/man/filterByExpr.html
rowData(spe_pseudo)$low_expr <- filterByExpr(spe_pseudo)
summary(rowData(spe_pseudo)$low_expr)
# Mode   FALSE    TRUE 
# logical   21059    7857 

save(sce_pseudobulk_bayesSpace, file = here::here("processed-data","rdata","spe","09_region_differential_expression",paste0("sce_pseudobulk_bayesSpace_k",k,".Rdata")))

