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

dim(spe_pseudo)
#[1] 28916    90

#find a good expression cutoff using edgeR::filterByExpr https://rdrr.io/bioc/edgeR/man/filterByExpr.html
rowData(spe_pseudo)$low_expr <- filterByExpr(spe_pseudo)
summary(rowData(spe_pseudo)$low_expr)
# Mode   FALSE    TRUE 
# logical   21059    7857 

spe_pseudo <- spe_pseudo[-which(rowData(spe_pseudo)$low_expr == TRUE),]
dim(spe_pseudo)
#[1] 21059    90

save(spe_pseudo, file = here::here("processed-data","rdata","spe","09_region_differential_expression",paste0("sce_pseudobulk_bayesSpace_k",k,".Rdata")))

#run PCA
pca_pseudo<- prcomp(t(assays(spe_pseudo)$logcounts))$x[, 1:50]

reducedDims(spe_pseudo) <- list(PCA=pca_pseudo, test=pca_pseudo) #need to figure out to add just one reduced dim https://bioconductor.org/packages/devel/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html


###plot PCA###
pdf(file = here::here("plots","09_region_differential_expression",paste0("sce_pseudobulk_pca_k",k,".Rdata")))
plotPCA(spe_pseudo, colour_by = "region", ncomponents = 8)
dev.off()


########installations###########
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("spatialLIBD")
# 
# install.packages("here")
# 
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("edgeR")
# 
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("scuttle")
