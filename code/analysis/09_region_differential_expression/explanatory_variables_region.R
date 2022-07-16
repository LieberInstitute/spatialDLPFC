library(SpatialExperiment)
library(spatialLIBD)
library(here)
library(edgeR)
library(scuttle)
# library(scRNAseq)
library(scater)

k <- as.numeric(Sys.getenv("SGE_TASK_ID"))

load(file = here::here("processed-data", "rdata", "spe", "pseudo_bulked_spe", paste0("spe_pseudobulk_bayesSpace_normalized_filtered_region_k", k, ".Rdata")))


# run PCA
pca <- prcomp(t(assays(spe_pseudo_filter_region)$logcounts)) # this will be computed in script that creates pseudobulked object
jaffelab::getPcaVars(pca)[seq_len(50)]
#  [1] 14.600  4.520  3.150  1.390  1.220  1.080  0.942  0.929  0.889  0.873
# [11]  0.858  0.855  0.845  0.831  0.818  0.810  0.800  0.796  0.787  0.785
# [21]  0.774  0.768  0.765  0.754  0.748  0.742  0.737  0.725  0.720  0.716
# [31]  0.706  0.704  0.695  0.691  0.682  0.671  0.670  0.666  0.657  0.647
# [41]  0.641  0.628  0.625  0.615  0.606  0.605  0.595  0.584  0.578  0.575
pca_pseudo <- pca$x[, seq_len(50)]


reducedDims(spe_pseudo_filter_region) <- list(PCA = pca_pseudo)
# define variables for which we are going to explain variance https://github.com/LieberInstitute/Visium_IF_AD/blob/5e3518a9d379e90f593f5826cc24ec958f81f4aa/code/11_grey_matter_only/02_explore_expr_variability.R#L61-L69

# code adapted from: http://bioconductor.org/packages/release/bioc/vignettes/scater/inst/doc/overview.html#2_Diagnostic_plots_for_quality_control
### plot PCA###
# make this a for loophttps://github.com/LieberInstitute/Visium_IF_AD/blob/5e3518a9d379e90f593f5826cc24ec958f81f4aa/code/11_grey_matter_only/02_explore_expr_variability.R#L73-L88
pdf(file = here::here("plots", "09_region_differential_expression", paste0("spe_pseudobulk_pca_k", k, ".pdf")), width = 14, height = 14)
plotPCA(spe_pseudo_filter_region, colour_by = "subject", ncomponents = 12, point_size = 1)
plotPCA(spe_pseudo_filter_region, colour_by = "region", ncomponents = 12, point_size = 1)
plotPCA(spe_pseudo_filter_region, colour_by = "sex", ncomponents = 12, point_size = 1)
plotPCA(spe_pseudo_filter_region, colour_by = "age", ncomponents = 12, point_size = 1)
plotPCA(spe_pseudo_filter_region, colour_by = "BayesSpace", ncomponents = 12, point_size = 1)
plotPCA(spe_pseudo_filter_region, colour_by = "sample_id", ncomponents = 12, point_size = 1)
dev.off()

# use these to plotPCA function https://github.com/LieberInstitute/Visium_IF_AD/blob/5e3518a9d379e90f593f5826cc24ec958f81f4aa/code/11_grey_matter_only/02_explore_expr_variability.R#L80-L81
#### plot explanatory variables ####

# uses linear regression model

vars <- getVarianceExplained(spe_pseudo_filter_region,
    variables = c("subject", "region", "sex", "age", "BayesSpace", "sample_id")
)
head(vars)

# subject    region       sex         age
# ENSG00000243485  6.149209 0.8049553 0.9891498 0.684494271
# ENSG00000238009 22.059419 0.3825081 6.0221979 0.758922808
# ENSG00000239945 10.112360 2.2471910 0.7490637 0.006556772
# ENSG00000241860 15.910843 1.7697885 2.7770465 3.317220956
# ENSG00000229905 10.113068 0.5911937 1.0110245 0.266643458
# ENSG00000237491 17.499031 9.9500577 0.2762791 0.570564318

pdf(file = here::here("plots", "09_region_differential_expression", paste0("plot_explanatory_vars_region_k", k, ".pdf")))
plotExplanatoryVariables(vars)
dev.off()

## For the spatialLIBD shiny app
rowData(spe_pseudo_filter_region)$gene_search <-
    paste0(
        rowData(spe_pseudo_filter_region)$gene_name,
        "; ",
        rowData(spe_pseudo_filter_region)$gene_id
    )

## Drop things we don't need
spatialCoords(spe_pseudo_filter_region) <- NULL
imgData(spe_pseudo_filter_region) <- NULL

## Simplify the colData()  for the pseudo-bulked data
colData(spe_pseudo_filter_region) <- colData(spe_pseudo_filter_region)[, sort(c(
    "age",
    "sample_id",
    "bayesSpace_harmony_9",
    "subject",
    "sex",
    "diagnosis"
))]

# ## Load pathology colors
# ## This info is used by spatialLIBD v1.7.18 or newer
# source(here("code", "colors_pathology.R"), echo = TRUE, max.deparse.length = 500)
# spe_pseudo$path_groups_colors <- colors_pathology[as.character(spe_pseudo$path_groups)]

## save RDS file
saveRDS(
    spe_pseudo_filter_region,
    file = here::here("processed-data", "rdata", "spe", "pseudo_bulked_spe", paste0("sce_pseudobulk_bayesSpace_normalized_filtered_region_k", k, "_shiny.Rdata"))
)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
