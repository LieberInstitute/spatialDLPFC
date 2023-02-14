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
spe <- readRDS("spe_subset_for_spatialLIBD.rds")
# load the pseudobulked object sce_pseudo
sce_pseudo <- readRDS("sce_pseudo_BayesSpace_k09.rds")
# load modeling results for k9 clustering/pseudobulking
load("modeling_results_BayesSpace_k09.Rdata", verbose = TRUE)
load("sig_genes_subset_k09.Rdata", verbose = TRUE)

spe$BayesSpace <- spe$BayesSpace_harmony_09
vars <- colnames(colData(spe))
# https://github.com/LieberInstitute/Visium_IF_AD/blob/5e3518a9d379e90f593f5826cc24ec958f81f4aa/code/05_deploy_app_wholegenome/app.R#L61-L72

colors_BayesSpace <- Polychrome::palette36.colors(28)
names(colors_BayesSpace) <- c(1:28)
m <- match(as.character(spe$BayesSpace_harmony_09), names(colors_BayesSpace))
stopifnot(all(!is.na(m)))
spe$BayesSpace_colors <- spe$BayesSpace_harmony_09_colors <- colors_BayesSpace[m]

## Deploy the website
spatialLIBD::run_app(
    spe,
    sce_layer = sce_pseudo,
    modeling_results = modeling_results,
    sig_genes = sig_genes,
    title = "spatialDLPFC, Visium, Sp09",
    spe_discrete_vars = c( # this is the variables for the spe object not the sce_pseudo object
        "BayesSpace",
        "ManualAnnotation",
        vars[grep("^SpaceRanger_|^scran_", vars)],
        vars[grep("^BayesSpace_harmony", vars)],
        vars[grep("^BayesSpace_pca", vars)],
        "graph_based_PCA_within",
        "PCA_SNN_k10_k7",
        "Harmony_SNN_k10_k7",
        "manual_layer_label",
        "wrinkle_type",
        "BayesSpace_colors"
    ),
    spe_continuous_vars = c(
        "sum_umi",
        "sum_gene",
        "expr_chrM",
        "expr_chrM_ratio",
        vars[grep("^VistoSeg_", vars)],
        vars[grep("^layer_", vars)],
        vars[grep("^broad_", vars)]
    ),
    default_cluster = "BayesSpace",
    docs_path = "www"
)
