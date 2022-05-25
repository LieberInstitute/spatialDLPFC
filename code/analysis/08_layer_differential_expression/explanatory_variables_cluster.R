library(SpatialExperiment)
library(spatialLIBD)
library(here)
library(edgeR)
library(scuttle)
library(sessioninfo)
library(scater)

k <- as.numeric(Sys.getenv("SGE_TASK_ID"))

spe_pseudo_filter_cluster <- readRDS(
  file = here::here("processed-data",
                    "rdata",
                    "spe",
                    "pseudo_bulked_spe",
                    paste0("spe_pseudobulk_bayesSpace_normalized_filtered_cluster_k",k,".RDS")))


#define variables for which we are going to explain variance 
#https://github.com/LieberInstitute/Visium_IF_AD/blob/5e3518a9d379e90f593f5826cc24ec958f81f4aa/code/11_grey_matter_only/02_explore_expr_variability.R#L61-L69
vars <- c(
  "age",
  "sample_id",
  "BayesSpace",
  "subject",
  "sex",
  "region"
)

# code adapted from: http://bioconductor.org/packages/release/bioc/vignettes/scater/inst/doc/overview.html#2_Diagnostic_plots_for_quality_control
###plot PCA###
#make this a for loop
#https://github.com/LieberInstitute/Visium_IF_AD/blob/5e3518a9d379e90f593f5826cc24ec958f81f4aa/code/11_grey_matter_only/02_explore_expr_variability.R#L73-L88
source(here("code", "analysis","colors_bayesSpace.R"), echo = TRUE, max.deparse.length = 500)
pdf(file = here::here("plots",
                      "08_layer_differential_expression",
                      paste0("plot_explanatory_variables_k",k,".pdf")), 
    width = 14, height = 14)

for (var in vars) {
  p <- plotPCA(
    spe_pseudo_filter_cluster,
    colour_by = var,
    ncomponents = 12,
    point_size = 1,
    label_format = c("%s %02i", " (%i%%)"),
    percentVar = metadata(spe_pseudo_filter_cluster)$PCA_var_explained
  )
  if (var == "BayesSpace") {
    p <- p + scale_color_manual("BayesSpace", values = colors_bayesSpace)
  }
  print(p)
}
dev.off()

pdf(file = here::here("plots","09_region_differential_expression",paste0("spe_pseudobulk_pca_k",k,".pdf")), width = 14, height = 14)
plotPCA(spe_pseudo_filter_cluster, colour_by = "subject", ncomponents = 12, point_size = 1) 
plotPCA(spe_pseudo_filter_cluster, colour_by = "region", ncomponents = 12, point_size = 1) 
plotPCA(spe_pseudo_filter_cluster, colour_by = "sex", ncomponents = 12, point_size = 1) 
plotPCA(spe_pseudo_filter_cluster, colour_by = "age", ncomponents = 12, point_size = 1) 
plotPCA(spe_pseudo_filter_cluster, colour_by = "BayesSpace", ncomponents = 12, point_size = 1) 
plotPCA(spe_pseudo_filter_cluster, colour_by = "sample_id", ncomponents = 12, point_size = 1)
dev.off()

#use these to plotPCA function https://github.com/LieberInstitute/Visium_IF_AD/blob/5e3518a9d379e90f593f5826cc24ec958f81f4aa/code/11_grey_matter_only/02_explore_expr_variability.R#L80-L81
####plot explanatory variables ####

#uses linear regression model

vars <- getVarianceExplained(spe_pseudo_filter_cluster, 
                             variables=c("subject", "region", "sex", "age", "BayesSpace","sample_id")) 
head(vars)

# subject    region       sex         age
# ENSG00000243485  6.149209 0.8049553 0.9891498 0.684494271
# ENSG00000238009 22.059419 0.3825081 6.0221979 0.758922808
# ENSG00000239945 10.112360 2.2471910 0.7490637 0.006556772
# ENSG00000241860 15.910843 1.7697885 2.7770465 3.317220956
# ENSG00000229905 10.113068 0.5911937 1.0110245 0.266643458
# ENSG00000237491 17.499031 9.9500577 0.2762791 0.570564318

pdf(file = here::here("plots","09_region_differential_expression",paste0("plot_explanatory_vars_region_k",k,".pdf")))
plotExplanatoryVariables(vars)
dev.off()

## For the spatialLIBD shiny app
rowData(spe_pseudo_filter_cluster)$gene_search <-
  paste0(
    rowData(spe_pseudo_filter_cluster)$gene_name,
    "; ",
    rowData(spe_pseudo_filter_cluster)$gene_id
  )

## Drop things we don't need
spatialCoords(spe_pseudo_filter_cluster) <- NULL
imgData(spe_pseudo_filter_cluster) <- NULL

## Simplify the colData()  for the pseudo-bulked data
colData(spe_pseudo_filter_cluster) <- colData(spe_pseudo_filter_cluster)[, sort(c(
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
  spe_pseudo_filter_cluster,
  file = here::here("processed-data","rdata","spe","pseudo_bulked_spe",paste0("sce_pseudobulk_bayesSpace_normalized_filtered_region_k",k,"_shiny.Rdata"))
)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()