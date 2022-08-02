library(SpatialExperiment)
library(spatialLIBD)
library(here)
library(edgeR)
library(scuttle)
library(sessioninfo)
library(scater)

k <- as.numeric(Sys.getenv("SGE_TASK_ID"))

spe_pseudo_filter_region <- readRDS(
  file = here::here(
    "processed-data",
    "rdata",
    "spe",
    "pseudo_bulked_spe",
    paste0("spe_pseudobulk_bayesSpace_normalized_filtered_region_k", k, ".RDS")
  )
)


# define variables for which we are going to explain variance
# https://github.com/LieberInstitute/Visium_IF_AD/blob/5e3518a9d379e90f593f5826cc24ec958f81f4aa/code/11_grey_matter_only/02_explore_expr_variability.R#L61-L69
vars <- c(
  "age",
  "sample_id",
  "BayesSpace",
  "subject",
  "sex",
  "region"
)

# code adapted from: http://bioconductor.org/packages/release/bioc/vignettes/scater/inst/doc/overview.html#2_Diagnostic_plots_for_quality_control
### plot PCA###
# make this a for loop
# https://github.com/LieberInstitute/Visium_IF_AD/blob/5e3518a9d379e90f593f5826cc24ec958f81f4aa/code/11_grey_matter_only/02_explore_expr_variability.R#L73-L88
source(here("code", "analysis", "colors_bayesSpace.R"), echo = TRUE, max.deparse.length = 500)
pdf(
  file = here::here(
    "plots",
    "09_regionr_differential_expression",
    paste0("pca_explanatory_variables_k", k, ".pdf")
  ),
  width = 14, height = 14
)

for (var in vars) {
  p <- plotPCA(
    spe_pseudo_filter_region,
    colour_by = var,
    ncomponents = 12,
    point_size = 1,
    label_format = c("%s %02i", " (%i%%)"),
    percentVar = metadata(spe_pseudo_filter_region)$PCA_var_explained
  )
  if (var == "BayesSpace") {
    p <- p + scale_color_manual("BayesSpace", values = colors_bayesSpace)
  }
  print(p)
}
dev.off()


# uses linear regression model
vars <- getVarianceExplained(spe_pseudo_filter_region,
                             variables = c("subject", "region", "sex", "age", "BayesSpace", "sample_id")
)
head(vars)
pdf(file = here::here(
  "plots",
  "09_region_differential_expression",
  paste0("plot_explanatory_vars_region_k", k, ".pdf")
))
plotExplanatoryVariables(vars)
dev.off()


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
