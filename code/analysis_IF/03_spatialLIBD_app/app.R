library("spatialLIBD")
library("markdown") ## Hm... to avoid this error
# 2021-11-11T05:30:49.941401+00:00 shinyapps[5096402]: Listening on http://127.0.0.1:32863
# 2021-11-11T05:30:50.218127+00:00 shinyapps[5096402]: Warning: Error in loadNamespace: there is no package called ‘markdown’
# 2021-11-11T05:30:50.222437+00:00 shinyapps[5096402]:   111: <Anonymous>
library("Polychrome")

## spatialLIBD uses golem
options("golem.app.prod" = TRUE)

## You need this to enable shinyapps to install Bioconductor packages
options(repos = BiocManager::repositories())

## Load the spe object
spe_IF <- readRDS("spe.rds")

# speB$BayesSpace <- speB$spatial.cluster
# speB$BayesSpace_initial <- speB$cluster.init
vars <- colnames(colData(spe_IF))

## Deploy the website
spatialLIBD::run_app(
    spe_IF,
    sce_layer = NULL,
    modeling_results = NULL,
    sig_genes = NULL,
    title = "spatialDLPFC, Visium SPG",
    spe_discrete_vars = c(
        vars[grep("^10x_", vars)],
        "ManualAnnotation"
    ),
    spe_continuous_vars = c(
        vars[grep("^(broad|layer|cart)_", vars)],
        "sum_umi",
        "sum_gene",
        "expr_chrM",
        "expr_chrM_ratio",
        "VistoSeg_count_deprecated",
        "cellpose_count"
    ),
    default_cluster = "10x_graphclust",
    docs_path = "www"
)
