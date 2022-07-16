library("spatialLIBD")
library("lobstr")
library("here")
library("sessioninfo")

## Load the spe object
load(here::here("code", "deploy_app", "spe_subset.Rdata"),
    verbose = TRUE
)

# load the pseudobulked object spe_pseudo
spe_pseudo <-
    readRDS(
        here(
            "code",
            "deploy_app_k16",
            "spe_pseudobulk_bayesSpace_normalized_filtered_cluster_k16.RDS"
        )
    )

lobstr::obj_size(spe_pseudo)
# 75.50 MB

# load modeling results for k16 clustering/pseudobulking
load(here(
    "code",
    "deploy_app_k16",
    "parsed_modeling_results_k16.Rdata"
),
verbose = TRUE
)
lobstr::obj_size(modeling_results)
# 33.28 MB

## For sig_genes_extract_all() to work https://github.com/LieberInstitute/Visium_IF_AD/blob/5e3518a9d379e90f593f5826cc24ec958f81f4aa/code/05_deploy_app_wholegenome/app.R#L37-L44
spe_pseudo$spatialLIBD <- spe_pseudo$BayesSpace

## Check that we have the right number of tests
k <- 16
tests <-
    lapply(modeling_results, function(x) {
        colnames(x)[grep("stat", colnames(x))]
    })
stopifnot(length(tests$anova) == 1) ## assuming only noWM
stopifnot(length(tests$enrichment) == k)
stopifnot(length(tests$pairwise) == choose(k, 2))

sig_genes <- sig_genes_extract_all(
    n = nrow(spe_pseudo),
    modeling_results = modeling_results,
    sce_layer = spe_pseudo
)

## Check that we have the right number of tests.
## the + 1 at the end assumes only noWM
stopifnot(length(unique(sig_genes$test)) == choose(k, 2) * 2 + k + 1)

## Fix the pairwise colors
sig_genes$test <- gsub("BayesSpace", "", sig_genes$test)

lobstr::obj_size(sig_genes)
# 2.77 GB

dim(sig_genes)
# [1] 2463859      12

## Drop parts we don't need to reduce the memory
sig_genes$in_rows <- NULL
sig_genes$in_rows_top20 <- NULL
lobstr::obj_size(sig_genes)
# 199.40 MB

# ## Subset sig_genes
# sig_genes <- subset(sig_genes, fdr < 0.05)
# dim(sig_genes)
# # [1] 679801     10
# lobstr::obj_size(sig_genes)
# # 59.64 MB

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

save(sig_genes,
    file = here::here("code", "deploy_app_k16", "sig_genes_subset_k16.Rdata")
)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
