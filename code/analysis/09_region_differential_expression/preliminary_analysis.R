library(SpatialExperiment)
library(spatialLIBD)
library(here)
library(edgeR)
library(scuttle)
#library(scRNAseq)
library(scater)

k <- as.numeric(Sys.getenv("SGE_TASK_ID"))

#load spe object and clusters
# load(file = here::here("processed-data","rdata","spe","01_build_spe","spe_filtered_final.Rdata"),verbose = TRUE)
# spe <- cluster_import(
#   spe,
#   cluster_dir = here::here("processed-data", "rdata", "spe", "clustering_results"), 
#   prefix = ""
# )

load("/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/01_build_spe/spe_filtered_final_with_clusters.Rdata")

## Pseudo-bulk for our current BayesSpace cluster results
spe_pseudo <- aggregateAcrossCells(
  spe,
  DataFrame(
    BayesSpace = colData(spe)[[paste0("bayesSpace_harmony_",k)]],
    sample_id = spe$sample_id
  )
)
spe_pseudo$BayesSpace <- factor(spe_pseudo$BayesSpace)

# add this line https://github.com/LieberInstitute/Visium_IF_AD/blob/5e3518a9d379e90f593f5826cc24ec958f81f4aa/code/11_grey_matter_only/01_create_pseudobulk_data.R#L125

#save here for differential expression analysis
save(spe_pseudo, file = here::here("processed-data","rdata","spe","pseudo_bulked_spe",paste0("spe_pseudobulk_bayesSpace_k",k,".Rdata")))

#adapt this code to keep colData variables I want https://github.com/LieberInstitute/Visium_IF_AD/blob/5e3518a9d379e90f593f5826cc24ec958f81f4aa/code/11_grey_matter_only/01_create_pseudobulk_data.R#L181-L196

#drop samples with low numbers of spots here https://github.com/LieberInstitute/Visium_IF_AD/blob/5e3518a9d379e90f593f5826cc24ec958f81f4aa/code/11_grey_matter_only/01_create_pseudobulk_data.R#L207
summary(spe_pseudo$ncells)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 80.0   394.8  1305.0  1898.8  3391.5  4629.0 
#spe_pseudo <- spe_pseudo[, spe_pseudo$ncells >=15]

# write csv of this https://github.com/LieberInstitute/Visium_IF_AD/blob/5e3518a9d379e90f593f5826cc24ec958f81f4aa/code/11_grey_matter_only/01_create_pseudobulk_data.R#L200

#find a good expression cutoff using edgeR::filterByExpr https://rdrr.io/bioc/edgeR/man/filterByExpr.html
rowData(spe_pseudo)$filter_expr <- filterByExpr(spe_pseudo) #add group by region keep <- filterByExpr(y, group=current$tomato)
summary(rowData(spe_pseudo)$filter_expr)
# Mode   FALSE    TRUE 
# logical   17914   11002 

rowData(spe_pseudo)$filter_expr_group_sample_id <- filterByExpr(spe_pseudo, group = spe_pseudo$sample_id)
summary(rowData(spe_pseudo)$filter_expr_group_sample_id )
# Mode   FALSE    TRUE 
# logical   13824   15092 

rowData(spe_pseudo)$filter_expr_group_cluster <- filterByExpr(spe_pseudo, group = colData(spe_pseudo)[[paste0("bayesSpace_harmony_",k)]])
summary(rowData(spe_pseudo)$filter_expr_group_cluster )
# Mode   FALSE    TRUE 
# logical   16762   12154 

rowData(spe_pseudo)$filter_expr_group_region <- filterByExpr(spe_pseudo, group = spe_pseudo$region)
summary(rowData(spe_pseudo)$filter_expr_group_region )

with(rowData(spe_pseudo), table(filter_expr,filter_expr_group_cluster))
# filter_expr_group_cluster
# filter_expr FALSE  TRUE
# FALSE 16788  1126
# TRUE      0 11002

with(rowData(spe_pseudo), table(filter_expr_group_sample_id, filter_expr_group_cluster))
# low_expr_group_cluster
# low_expr_group_sample_id FALSE  TRUE
# FALSE 15196     0
# TRUE   1566 12154

spe_pseudo_filter_cluster <- spe_pseudo[which(rowData(spe_pseudo)$filter_expr_group_cluster),]
dim(spe_pseudo_filter_cluster)

spe_pseudo_filter_region <- spe_pseudo[which(rowData(spe_pseudo)$filter_expr_group_region),]
dim(spe_pseudo_filter_region)

#log normalize the counts
# spe_pseudo <- logNormCounts(spe_pseudo,size.factors = NULL)

#replace calcNormFactors code with this https://github.com/LieberInstitute/Visium_IF_AD/blob/5e3518a9d379e90f593f5826cc24ec958f81f4aa/code/11_grey_matter_only/01_create_pseudobulk_data.R#L254-L258
x <- edgeR::cpm(edgeR::calcNormFactors(spe_pseudo_filter_cluster), 
                log = TRUE, 
                prior.count = 1)
## Verify that the gene order hasn't changed
stopifnot(identical(rownames(x), rownames(spe_pseudo_filter_cluster)))
## Fix the column names. DGEList will have samples names as Sample1 Sample2 etc
dimnames(x) <- dimnames(spe_pseudo_filter_cluster)
## Store the log normalized counts on the SingleCellExperiment object
logcounts(spe_pseudo_filter_cluster) <- x
## We don't need this 'x' object anymore
rm(x)

dim(spe_pseudo_filter_cluster)
#[1] 12128    60

#remove parts of pseudobulked objects to make them smaller https://github.com/LieberInstitute/Visium_IF_AD/blob/5e3518a9d379e90f593f5826cc24ec958f81f4aa/code/11_grey_matter_only/01_create_pseudobulk_data.R#L389-L399

# adpate this to use the polychrome scale I'm already using for my other figures https://github.com/LieberInstitute/Visium_IF_AD/blob/5e3518a9d379e90f593f5826cc24ec958f81f4aa/code/11_grey_matter_only/01_create_pseudobulk_data.R#L403-L404

#update this to indicate this version of the object is normalized and filtered
save(spe_pseudo_filter_cluster, file = here::here("processed-data","rdata","spe","pseudo_bulked_spe",paste0("spe_pseudobulk_bayesSpace_normalized_filtered_cluster_k",k,".Rdata")))

#log normalize the counts
# spe_pseudo <- logNormCounts(spe_pseudo,size.factors = NULL)
x <- edgeR::cpm(edgeR::calcNormFactors(spe_pseudo_filter_region), log = TRUE, prior.count = 1)
## Verify that the gene order hasn't changed
stopifnot(identical(rownames(x), rownames(spe_pseudo_filter_region)))
## Fix the column names. DGEList will have samples names as Sample1 Sample2 etc
dimnames(x) <- dimnames(spe_pseudo_filter_region)
## Store the log normalized counts on the SingleCellExperiment object
logcounts(spe_pseudo_filter_region) <- x
## We don't need this 'x' object anymore
rm(x)

dim(spe_pseudo_filter_region)
#[1] 28916    90

save(spe_pseudo_filter_region, file = here::here("processed-data","rdata","spe","pseudo_bulked_spe",paste0("spe_pseudobulk_bayesSpace_normalized_filtered_region_k",k,".Rdata")))



## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
