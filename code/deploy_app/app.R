library("spatialLIBD")
library("markdown") ## Hm... to avoid this error
# 2021-11-11T05:30:49.941401+00:00 shinyapps[5096402]: Listening on http://127.0.0.1:32863
# 2021-11-11T05:30:50.218127+00:00 shinyapps[5096402]: Warning: Error in loadNamespace: there is no package called ‘markdown’
# 2021-11-11T05:30:50.222437+00:00 shinyapps[5096402]:   111: <Anonymous>
library("scater") ## to compute some reduced dimensions

## spatialLIBD uses golem
options("golem.app.prod" = TRUE)

## You need this to enable shinyapps to install Bioconductor packages
options(repos = BiocManager::repositories())

## Load the data
load("spe_filtered_final_with_clusters.Rdata", verbose = TRUE)
#load spe_pseudo
spe_pseudo <- readRDS("spe_pseudobulk_bayesSpace_normalized_filtered_cluster_k9.RDS")
#load modeling results
load("parsed_modeling_results_k9.Rdata",verbose = TRUE)


## For sig_genes_extract_all() to work https://github.com/LieberInstitute/Visium_IF_AD/blob/5e3518a9d379e90f593f5826cc24ec958f81f4aa/code/05_deploy_app_wholegenome/app.R#L37-L44
spe_pseudo$spatialLIBD <- spe_pseudo$BayesSpace
sig_genes <- sig_genes_extract_all(
  n = nrow(spe_pseudo),
  modeling_results = modeling_results,
  sce_layer = spe_pseudo
)
sig_genes <- sig_genes_extract_all(
  n = nrow(spe_pseudo),
  modeling_results = test,
  sce_layer = spe_pseudo
)

## Extract FDR < 5%
## From
## https://github.com/LieberInstitute/brainseq_phase2/blob/be2b7f972bb2a0ede320633bf06abe1d4ef2c067/supp_tabs/create_supp_tables.R#L173-L181
# fix_csv <- function(df) {
#   for (i in seq_len(ncol(df))) {
#     if (any(grepl(",", df[, i]))) {
#       message(paste(Sys.time(), "fixing column", colnames(df)[i]))
#       df[, i] <- gsub(",", ";", df[, i])
#     }
#   }
#   return(df)
# }
# z <- fix_csv(as.data.frame(subset(sig_genes, fdr < 0.05)))
# write.csv(z, file = file.path(dir_rdata, "Visium_IF_AD_wholegenome_model_results_FDR5perc.csv"))

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
