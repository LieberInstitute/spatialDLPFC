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
#load("spe_merged_final_nocounts.Rdata", verbose = TRUE)
load("spe_filtered_final_with_clusters.Rdata", verbose = TRUE)
## If you want to see the unfiltered data locally (it's too big for shinyapps.io)
# load(here::here("processed-data", "rdata", "spe", "spe_merged_final.Rdata"), verbose = TRUE)


## Don't include scran_quick_cluster since it's too big (288 levels)
spe$scran_quick_cluster <- NULL

## Rename 'count' to make it more intuitive
#spe$cell_count <- spe$count
#spe$count <- NULL

#load spe_pseudo
load("sce_pseudobulk_bayesSpace_normalized_filtered_k9.Rdata",verbose = TRUE)

#load modeling results
load("parsed_modeling_results_k9.Rdata",verbose = TRUE)

## For sig_genes_extract_all() to work https://github.com/LieberInstitute/Visium_IF_AD/blob/5e3518a9d379e90f593f5826cc24ec958f81f4aa/code/05_deploy_app_wholegenome/app.R#L37-L44
spe_pseudo$spatialLIBD <- spe_pseudo$bayesSpace_harmony_9
# sig_genes <- sig_genes_extract_all(
#   n = nrow(spe_pseudo),
#   modeling_results = modeling_results,
#   sce_layer = spe_pseudo
# )

vars <- colnames(colData(spe))

## Deploy the website
spatialLIBD::run_app(
    spe,
    sce_layer = spe_pseudo,
    modeling_results = modeling_results,
    sig_genes = NULL, #change this to use sig_genes object
    title = "spatialDLPFC, Spangler et al, 2021",
    spe_discrete_vars = c(
        vars[grep("10x_|scran_", vars)],
        "ManualAnnotation",
        "bayesSpace_harmony_9",
        "bayesSpace_harmony_16",
        "bayesSpace_harmony_28"
    ),
    spe_continuous_vars = c(
        "sum_umi",
        "sum_gene",
        "expr_chrM",
        "expr_chrM_ratio",
        "count"
    ),
    default_cluster = "bayesSpace_harmony_9"
)
