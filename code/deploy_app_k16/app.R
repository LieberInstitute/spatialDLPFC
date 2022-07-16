library("spatialLIBD")
library("markdown") ## Hm... to avoid this error
# 2021-11-11T05:30:49.941401+00:00 shinyapps[5096402]: Listening on http://127.0.0.1:32863
# 2021-11-11T05:30:50.218127+00:00 shinyapps[5096402]: Warning: Error in loadNamespace: there is no package called ‘markdown’
# 2021-11-11T05:30:50.222437+00:00 shinyapps[5096402]:   111: <Anonymous>

## spatialLIBD uses golem
options("golem.app.prod" = TRUE)

## You need this to enable shinyapps to install Bioconductor packages
options(repos = BiocManager::repositories())

## Load the spe object
load("spe_subset.Rdata", verbose = TRUE)
# load the pseudobulked object spe_pseudo
spe_pseudo <- readRDS("spe_pseudobulk_bayesSpace_normalized_filtered_cluster_k16.RDS")
# load modeling results for k9 clustering/pseudobulking
load("parsed_modeling_results_k16.Rdata", verbose = TRUE)
load("sig_genes_subset_k16.Rdata", verbose = TRUE)

spe$BayesSpace <- spe$bayesSpace_harmony_16
vars <- colnames(colData(spe))
# https://github.com/LieberInstitute/Visium_IF_AD/blob/5e3518a9d379e90f593f5826cc24ec958f81f4aa/code/05_deploy_app_wholegenome/app.R#L61-L72

## Deploy the website
spatialLIBD::run_app(
    spe,
    sce_layer = spe_pseudo,
    modeling_results = modeling_results,
    sig_genes = sig_genes,
    title = "spatialDLPFC_k16, Spangler et al, 2022",
    spe_discrete_vars = c( # this is the variables for the spe object not the spe_pseudo object
        vars[grep("10x_|scran_", vars)],
        "ManualAnnotation",
        vars[grep("bayesSpace_harmony", vars)],
        vars[grep("bayesSpace_pca", vars)],
        "graph_based_PCA_within",
        "PCA_SNN_k10_k7",
        "Harmony_SNN_k10_k7",
        "BayesSpace"
    ),
    spe_continuous_vars = c(
        "sum_umi",
        "sum_gene",
        "expr_chrM",
        "expr_chrM_ratio",
        "count"
    ),
    default_cluster = "BayesSpace"
)
