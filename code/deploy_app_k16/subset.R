library("spatialLIBD")
library("lobstr")
library("here")
library("sessioninfo")

## Set up soft links if needed
withr::with_dir(
    here("code", "deploy_app_k16"),
    system("ln -s ../../processed-data/rdata/spe/07_layer_differential_expression/modeling_results_BayesSpace_k16.Rdata modeling_results_BayesSpace_k16.Rdata")
)
withr::with_dir(
    here("code", "deploy_app_k16"),
    system("ln -s ../../processed-data/rdata/spe/07_layer_differential_expression/sce_pseudo_BayesSpace_k16.rds sce_pseudo_BayesSpace_k16.rds")
)
withr::with_dir(
    here("code", "deploy_app_k16"),
    system("ln -s ../../processed-data/rdata/spe/01_build_spe/spe_subset_for_spatialLIBD.rds spe_subset_for_spatialLIBD.rds")
)
withr::with_dir(
    here("code", "deploy_app_k16"),
    system("ln -s ../deploy_app_k09/www www")
)

# load the pseudobulked object sce_pseudo
sce_pseudo <-
    readRDS(
        here(
            "code",
            "deploy_app_k16",
            "sce_pseudo_BayesSpace_k16.rds"
        )
    )

lobstr::obj_size(sce_pseudo)
# 75.05 MB

# load modeling results for k16 clustering/pseudobulking
load(
    here(
        "code",
        "deploy_app_k16",
        "modeling_results_BayesSpace_k16.Rdata"
    ),
    verbose = TRUE
)
lobstr::obj_size(modeling_results)
# 33.43 MB

## For sig_genes_extract_all() to work https://github.com/LieberInstitute/Visium_IF_AD/blob/5e3518a9d379e90f593f5826cc24ec958f81f4aa/code/05_deploy_app_wholegenome/app.R#L37-L44
sce_pseudo$spatialLIBD <- sce_pseudo$BayesSpace

## Check that we have the right number of tests
k <- 16
tests <- lapply(modeling_results, function(x) {
    colnames(x)[grep("stat", colnames(x))]
})
stopifnot(length(tests$anova) == 1) ## assuming only "all"
stopifnot(length(tests$enrichment) == k)
stopifnot(length(tests$pairwise) == choose(k, 2))

sig_genes <- sig_genes_extract_all(
    n = nrow(sce_pseudo),
    modeling_results = modeling_results,
    sce_layer = sce_pseudo
)

## Check that we have the right number of tests.
## the + 1 at the end assumes only "all"
stopifnot(length(unique(sig_genes$test)) == choose(k, 2) * 2 + k + 1)

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
# # [1] 342449     10
# lobstr::obj_size(sig_genes)
# # 28.37 MB

## Extract FDR < 5%
## From
## https://github.com/LieberInstitute/brainseq_phase2/blob/be2b7f972bb2a0ede320633bf06abe1d4ef2c067/supp_tabs/create_supp_tables.R#L173-L181
fix_csv <- function(df) {
    for (i in seq_len(ncol(df))) {
        if (any(grepl(",", df[, i]))) {
            message(paste(Sys.time(), "fixing column", colnames(df)[i]))
            df[, i] <- gsub(",", ";", df[, i])
        }
    }
    return(df)
}
z <- fix_csv(as.data.frame(subset(sig_genes, fdr < 0.05)))
dim(z)
# [1] 684783     10
dim(subset(z, top <= 25))
# [1] 6127   10
write.csv(
    subset(z, top <= 25),
    file = here(
        "processed-data",
        "rdata",
        "spe",
        "07_layer_differential_expression",
        "spatialDLPFC_model_results_FDR5perc_top25_k16.csv"
    )
)

save(sig_genes,
    file = here::here("code", "deploy_app_k16", "sig_genes_subset_k16.Rdata")
)


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
