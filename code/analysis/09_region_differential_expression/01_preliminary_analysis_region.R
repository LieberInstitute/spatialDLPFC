library(SpatialExperiment)
library(spatialLIBD)
library(here)
library(edgeR)
library(scuttle)
library(sessioninfo)
library(scater)

k <- as.numeric(Sys.getenv("SGE_TASK_ID"))

load("/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/01_build_spe/spe_filtered_final_with_clusters.Rdata")

## Pseudo-bulk for our current BayesSpace cluster results
spe_pseudo <- aggregateAcrossCells(
    spe,
    DataFrame(
        BayesSpace = colData(spe)[[paste0("bayesSpace_harmony_", k)]],
        sample_id = spe$sample_id
    )
)
spe_pseudo$BayesSpace <- factor(spe_pseudo$BayesSpace)

# https://github.com/LieberInstitute/Visium_IF_AD/blob/5e3518a9d379e90f593f5826cc24ec958f81f4aa/code/11_grey_matter_only/01_create_pseudobulk_data.R#L125
colnames(spe_pseudo) <- paste0(spe_pseudo$sample_id, "_", spe_pseudo$BayesSpace)


# adapt this code to keep colData variables I want
# https://github.com/LieberInstitute/Visium_IF_AD/blob/5e3518a9d379e90f593f5826cc24ec958f81f4aa/code/11_grey_matter_only/01_create_pseudobulk_data.R#L181-L196
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


summary(spe_pseudo$ncells)

# drop samples with low numbers of spots here
# https://github.com/LieberInstitute/Visium_IF_AD/blob/5e3518a9d379e90f593f5826cc24ec958f81f4aa/code/11_grey_matter_only/01_create_pseudobulk_data.R#L207
spe_pseudo <- spe_pseudo[, spe_pseudo$ncells >= 10]

rowData(spe_pseudo)$filter_expr_group_region <- filterByExpr(spe_pseudo, group = spe_pseudo$region)
summary(rowData(spe_pseudo)$filter_expr_group_region)

spe_pseudo_filter_region <- spe_pseudo[which(rowData(spe_pseudo)$filter_expr_group_region), ]
dim(spe_pseudo_filter_region)

# replace calcNormFactors code with this https://github.com/LieberInstitute/Visium_IF_AD/blob/5e3518a9d379e90f593f5826cc24ec958f81f4aa/code/11_grey_matter_only/01_create_pseudobulk_data.R#L254-L258
x <- edgeR::cpm(edgeR::calcNormFactors(spe_pseudo_filter_region),
    log = TRUE,
    prior.count = 1
)
## Verify that the gene order hasn't changed
stopifnot(identical(rownames(x), rownames(spe_pseudo_filter_region)))
## Fix the column names. DGEList will have samples names as Sample1 Sample2 etc
dimnames(x) <- dimnames(spe_pseudo_filter_region)
## Store the log normalized counts on the SingleCellExperiment object
logcounts(spe_pseudo_filter_region) <- x
## We don't need this 'x' object anymore
rm(x)

dim(spe_pseudo_filter_region)

# run PCA
pca <- prcomp(t(assays(spe_pseudo_filter_region)$logcounts)) # this will be computed in script that creates pseudobulked object
message(Sys.time(), " % of variance explained for the top 20 PCs:")
metadata(spe_pseudo_filter_region)
metadata(spe_pseudo_filter_region) <- list("PCA_var_explained" = jaffelab::getPcaVars(pca)[seq_len(20)])
metadata(spe_pseudo_filter_region)
pca_pseudo <- pca$x[, seq_len(50)]
colnames(pca_pseudo) <- paste0("PC", sprintf("%02d", seq_len(ncol(pca_pseudo))))
reducedDims(spe_pseudo_filter_region) <- list(PCA = pca_pseudo)

## Compute some reduced dims
set.seed(20220423)
spe_pseudo_filter_region <- scater::runMDS(spe_pseudo_filter_region, ncomponents = 20)
spe_pseudo_filter_region <- scater::runPCA(spe_pseudo_filter_region, name = "runPCA")


# https://github.com/LieberInstitute/Visium_IF_AD/blob/5e3518a9d379e90f593f5826cc24ec958f81f4aa/code/11_grey_matter_only/01_create_pseudobulk_data.R#L389-L399
## For the spatialLIBD shiny app
rowData(spe_pseudo_filter_region)$gene_search <-
    paste0(
        rowData(spe_pseudo_filter_region)$gene_name,
        "; ",
        rowData(spe_pseudo_filter_region)$gene_id
    )

# remove parts of pseudobulked objects to make them smaller
spatialCoords(spe_pseudo_filter_region) <- NULL
imgData(spe_pseudo_filter_region) <- NULL


# adpate this to use the polychrome scale I'm already using for my other figures
# https://github.com/LieberInstitute/Visium_IF_AD/blob/5e3518a9d379e90f593f5826cc24ec958f81f4aa/code/11_grey_matter_only/01_create_pseudobulk_data.R#L403-L404
source(here("code", "analysis", "colors_bayesSpace.R"), echo = TRUE, max.deparse.length = 500)
spe_pseudo_filter_region$BayesSpace_colors <- colors_bayesSpace[as.character(spe_pseudo_filter_region$BayesSpace)]

# update this to indicate this version of the object is normalized and filtered
saveRDS(spe_pseudo_filter_region,
    file = here::here("processed-data", "rdata", "spe", "pseudo_bulked_spe", paste0("spe_pseudobulk_bayesSpace_normalized_filtered_region_k", k, ".RDS"))
)


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
