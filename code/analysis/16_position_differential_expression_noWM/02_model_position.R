# library(sgejobs)
# sgejobs::job_single(
#     name = "02_model_position",
#     create_shell = TRUE,
#     queue = "bluejay",
#     memory = "5G",
#     task_num = 28,
#     tc = 10
# )
# To execute the script builder, use: qsub 02_model_position.sh

k <- as.numeric(Sys.getenv("SGE_TASK_ID"))

## For testing
if (FALSE) {
    k <- 9
}


library("here")
library("sessioninfo")
library("spatialLIBD")

## output directory
dir_input <- here::here(
    "processed-data",
    "rdata",
    "spe",
    "16_position_differential_expression_noWM"
)
dir_rdata <- here::here(
    "processed-data",
    "rdata",
    "spe",
    "16_position_differential_expression_noWM"
)
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)
stopifnot(file.exists(dir_rdata)) ## Check that it was created successfully
dir_plots <- here::here(
    "plots",
    "16_position_differential_expression_noWM"
)
dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)
stopifnot(file.exists(dir_plots))

## load sce_pseudo data
sce_pseudo <-
    readRDS(
        file.path(
            dir_input,
            paste0("sce_pseudo_BayesSpace_k", sprintf("%02d", k), ".rds")
        )
    )

## To avoid having to change parameters later on
sce_pseudo$registration_variable <- sce_pseudo$position
sce_pseudo$registration_sample_id <- sce_pseudo$sample_id

## Set arguments used in spatialLIBD::registration_wrapper()
covars <- c("BayesSpace", "age", "sex")
gene_ensembl <- "gene_id"
gene_name <- "gene_name"
suffix <- "all"

## Taken from spatialLIBD::registration_wrapper()
## https://github.com/LieberInstitute/spatialLIBD/blob/master/R/registration_wrapper.R
registration_mod <-
    registration_model(sce_pseudo, covars = covars)

block_cor <-
    registration_block_cor(sce_pseudo, registration_model = registration_mod)

results_enrichment <-
    registration_stats_enrichment(
        sce_pseudo,
        block_cor = block_cor,
        covars = covars,
        gene_ensembl = gene_ensembl,
        gene_name = gene_name
    )
results_pairwise <-
    registration_stats_pairwise(
        sce_pseudo,
        registration_model = registration_mod,
        block_cor = block_cor,
        gene_ensembl = gene_ensembl,
        gene_name = gene_name
    )
results_anova <-
    registration_stats_anova(
        sce_pseudo,
        block_cor = block_cor,
        covars = covars,
        gene_ensembl = gene_ensembl,
        gene_name = gene_name,
        suffix = suffix
    )

modeling_results <- list(
    "anova" = results_anova,
    "enrichment" = results_enrichment,
    "pairwise" = results_pairwise
)

## Save the final results
save(
    modeling_results,
    file = file.path(
        dir_rdata,
        paste0("modeling_results_position_k", sprintf("%02d", k), ".Rdata")
    )
)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
