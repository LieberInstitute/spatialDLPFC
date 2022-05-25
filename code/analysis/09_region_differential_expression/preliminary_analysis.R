library(SpatialExperiment)
library(spatialLIBD)
library(here)
library(edgeR)
library(scuttle)
library(sessioninfo)
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

# https://github.com/LieberInstitute/Visium_IF_AD/blob/5e3518a9d379e90f593f5826cc24ec958f81f4aa/code/11_grey_matter_only/01_create_pseudobulk_data.R#L125
colnames(spe_pseudo) <- paste0(spe_pseudo$sample_id, "_", spe_pseudo$BayesSpace)


#adapt this code to keep colData variables I want 
#https://github.com/LieberInstitute/Visium_IF_AD/blob/5e3518a9d379e90f593f5826cc24ec958f81f4aa/code/11_grey_matter_only/01_create_pseudobulk_data.R#L181-L196
colData(spe_pseudo) <- colData(spe_pseudo)[, sort(c(
  "age",
  "sample_id",
  "BayesSpace",
  "subject",
  "sex",
  "region",
  "diagnosis",
  "ncells"
))]


# write csv of pseudobulked colData in order to look at number of spots in each pseudobulked sample
#https://github.com/LieberInstitute/Visium_IF_AD/blob/5e3518a9d379e90f593f5826cc24ec958f81f4aa/code/11_grey_matter_only/01_create_pseudobulk_data.R#L200
write.csv(as.data.frame(colData(spe_pseudo)),file = here::here("processed-data","rdata","spe","pseudo_bulked_spe",paste0("pseudobulk_colData_k",k,".csv")))

summary(spe_pseudo$ncells)

#drop samples with low numbers of spots here 
#https://github.com/LieberInstitute/Visium_IF_AD/blob/5e3518a9d379e90f593f5826cc24ec958f81f4aa/code/11_grey_matter_only/01_create_pseudobulk_data.R#L207
spe_pseudo <- spe_pseudo[, spe_pseudo$ncells >=10]


rowData(spe_pseudo)$filter_expr_group_sample_id <- filterByExpr(spe_pseudo, group = spe_pseudo$sample_id)
summary(rowData(spe_pseudo)$filter_expr_group_sample_id )


rowData(spe_pseudo)$filter_expr_group_cluster <- filterByExpr(spe_pseudo, group = spe_pseudo$BayesSpace)
summary(rowData(spe_pseudo)$filter_expr_group_cluster )


rowData(spe_pseudo)$filter_expr_group_region <- filterByExpr(spe_pseudo, group = spe_pseudo$region)
summary(rowData(spe_pseudo)$filter_expr_group_region )

with(rowData(spe_pseudo), table(filter_expr,filter_expr_group_cluster))
with(rowData(spe_pseudo), table(filter_expr_group_sample_id, filter_expr_group_cluster))

spe_pseudo_filter_cluster <- spe_pseudo[which(rowData(spe_pseudo)$filter_expr_group_cluster),]
dim(spe_pseudo_filter_cluster)

spe_pseudo_filter_region <- spe_pseudo[which(rowData(spe_pseudo)$filter_expr_group_region),]
dim(spe_pseudo_filter_region)

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

#remove parts of pseudobulked objects to make them smaller 
#https://github.com/LieberInstitute/Visium_IF_AD/blob/5e3518a9d379e90f593f5826cc24ec958f81f4aa/code/11_grey_matter_only/01_create_pseudobulk_data.R#L389-L399
rowData(spe_pseudo_filter_cluster)$gene_search <-
  paste0(
    rowData(spe_pseudo_filter_cluster)$gene_name,
    "; ",
    rowData(spe_pseudo_filter_cluster)$gene_id
  )

## Drop things we don't need
spatialCoords(spe_pseudo_filter_cluster) <- NULL
imgData(spe_pseudo_filter_cluster) <- NULL


# adpate this to use the polychrome scale I'm already using for my other figures 
#https://github.com/LieberInstitute/Visium_IF_AD/blob/5e3518a9d379e90f593f5826cc24ec958f81f4aa/code/11_grey_matter_only/01_create_pseudobulk_data.R#L403-L404
source(here("code", "analysis","colors_bayesSpace.R"), echo = TRUE, max.deparse.length = 500)
spe_pseudo$BayesSpace_colors <- colors_bayesSpace[as.character(spe_pseudo$BayesSpace)]

#update this to indicate this version of the object is normalized and filtered
saveRDS(spe_pseudo_filter_cluster, 
        file = here::here("processed-data","rdata","spe","pseudo_bulked_spe",paste0("spe_pseudobulk_bayesSpace_normalized_filtered_cluster_k",k,".RDS")))


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
