library("spatialLIBD")
library("markdown") ## Hm... to avoid this error
# 2021-11-11T05:30:49.941401+00:00 shinyapps[5096402]: Listening on http://127.0.0.1:32863
# 2021-11-11T05:30:50.218127+00:00 shinyapps[5096402]: Warning: Error in loadNamespace: there is no package called ‘markdown’
# 2021-11-11T05:30:50.222437+00:00 shinyapps[5096402]:   111: <Anonymous>

## spatialLIBD uses golem
options("golem.app.prod" = TRUE)

## You need this to enable shinyapps to install Bioconductor packages
options(repos = BiocManager::repositories())

## Load the data
load("spe_merged_final_nocounts.Rdata", verbose = TRUE)
## If you want to see the unfiltered data locally (it's too big for shinyapps.io)
# load(here::here("processed-data", "rdata", "spe", "spe_merged_final.Rdata"), verbose = TRUE)

## Don't include scran_quick_cluster since it's too big (288 levels)
spe$scran_quick_cluster <- NULL

## Rename 'count' to make it more intuitive
spe$cell_count <- spe$count
spe$count <- NULL

vars <- colnames(colData(spe))

## Deploy the website
spatialLIBD::run_app(
    spe,
    sce_layer = NULL,
    modeling_results = NULL,
    sig_genes = NULL,
    title = "spatialDLPFC, Spangler et al, 2021",
    spe_discrete_vars = c(
        vars[grep("10x_|scran_", vars)],
        "ManualAnnotation",
        "BayesSpace",
        "BayesSpace_cluster.init"
    ),
    spe_continuous_vars = c(
        "sum_umi",
        "sum_gene",
        "expr_chrM",
        "expr_chrM_ratio",
        "cell_count"
    ),
    default_cluster = "BayesSpace"
)
