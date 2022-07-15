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
#load the pseudobulked object spe_pseudo
spe_pseudo <- readRDS("spe_pseudobulk_bayesSpace_normalized_filtered_cluster_k16.RDS")
#load modeling results for k9 clustering/pseudobulking
load("parsed_modeling_results_k16.Rdata",verbose = TRUE)


## For sig_genes_extract_all() to work https://github.com/LieberInstitute/Visium_IF_AD/blob/5e3518a9d379e90f593f5826cc24ec958f81f4aa/code/05_deploy_app_wholegenome/app.R#L37-L44
spe_pseudo$spatialLIBD <- spe_pseudo$BayesSpace

## Check that we have the right number of tests
k <- 16
tests <- lapply(modeling_results, function(x) { colnames(x)[grep("stat", colnames(x))]})
stopifnot(length(tests$anova) == 1) ## assuming only noWM
stopifnot(length(tests$enrichment) == k)
#stopifnot(length(tests$pairwise) == choose(k, 2) * 2) this doesn't work. not sure why.

sig_genes <- sig_genes_extract_all(
  n = 4000, #only extract the top 5000 significant genes or else it's too big
  modeling_results = modeling_results,
  sce_layer = spe_pseudo
)

#subset genes after for fdr less than 0.05.  Create this smaller sig_genes in different script and just load it here.
#lobstr::obj_size(modeling_results) 
#0.03099124 B
#lobstr::obj_size(sig_genes) 
#0.7587798 B

dim(sig_genes)
#[1] 1285000      12


## Check that we have the right number of tests.
## the + 1 at the end assumes only noWM
# stopifnot(length(unique(sig_genes$test)) == choose(k, 2) * 2 + k + 1)

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

spe$BayesSpace <- spe$bayesSpace_harmony_16
vars <- colnames(colData(spe))
#https://github.com/LieberInstitute/Visium_IF_AD/blob/5e3518a9d379e90f593f5826cc24ec958f81f4aa/code/05_deploy_app_wholegenome/app.R#L61-L72

## Deploy the website
spatialLIBD::run_app(
    spe,
    sce_layer = spe_pseudo,
    modeling_results = modeling_results,
    sig_genes = sig_genes,
    title = "spatialDLPFC_k16, Spangler et al, 2022",
    spe_discrete_vars = c( #this is the variables for the spe object not the spe_pseudo object
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
