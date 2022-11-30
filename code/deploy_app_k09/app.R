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
load("spe_subset.Rdata", verbose = TRUE)
# load the pseudobulked object sce_pseudo
sce_pseudo <- readRDS("sce_pseudo_BayesSpace_k09.rds")
# load modeling results for k9 clustering/pseudobulking
load("modeling_results_BayesSpace_k09.Rdata", verbose = TRUE)
load("sig_genes_subset_k09.Rdata", verbose = TRUE)

spe$BayesSpace <- spe$bayesSpace_harmony_9
vars <- colnames(colData(spe))
# https://github.com/LieberInstitute/Visium_IF_AD/blob/5e3518a9d379e90f593f5826cc24ec958f81f4aa/code/05_deploy_app_wholegenome/app.R#L61-L72

colors_bayesSpace <- Polychrome::palette36.colors(28)
names(colors_bayesSpace) <- c(1:28)
m <- match(as.character(spe$bayesSpace_harmony_9), names(colors_bayesSpace))
stopifnot(all(!is.na(m)))
spe$BayesSpace_colors <- spe$bayesSpace_harmony_9_colors <- colors_bayesSpace[m]

## Deploy the website
spatialLIBD::run_app(
    spe,
    sce_layer = sce_pseudo,
    modeling_results = modeling_results,
    sig_genes = sig_genes,
    title = "spatialDLPFC, Spangler et al, 2022",
    spe_discrete_vars = c( # this is the variables for the spe object not the sce_pseudo object
        vars[grep("10x_|scran_", vars)],
        "ManualAnnotation",
        vars[grep("bayesSpace_harmony", vars)],
        vars[grep("bayesSpace_pca", vars)],
        "graph_based_PCA_within",
        "PCA_SNN_k10_k7",
        "Harmony_SNN_k10_k7",
        "BayesSpace",
        "BayesSpace_colors"
    ),
    spe_continuous_vars = c(
        "sum_umi",
        "sum_gene",
        "expr_chrM",
        "expr_chrM_ratio",
        "count",
        vars[grep("^tangram_", vars)],
        vars[grep("^cell2location_", vars)],
        vars[grep("^spotlight_", vars)]
    ),
    default_cluster = "BayesSpace"
)
