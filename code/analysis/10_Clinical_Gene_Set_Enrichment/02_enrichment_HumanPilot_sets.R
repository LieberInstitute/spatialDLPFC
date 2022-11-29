# library(sgejobs)
# sgejobs::job_single(
#     name = "02_enrichment_HumanPilot_sets",
#     create_shell = TRUE,
#     queue = "bluejay",
#     memory = "5G"
# )
# To execute the script builder, use: qsub 02_enrichment_HumanPilot_sets.sh

library("here")
library("purrr")
library("sessioninfo")

## Input dir
dir_input <- here::here(
    "processed-data",
    "rdata",
    "spe",
    "07_layer_differential_expression"
)

## Output directories
dir_rdata <- here::here(
    "processed-data",
    "rdata",
    "spe",
    "10_Clinical_Gene_Set_Enrichment"
)
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)
stopifnot(file.exists(dir_rdata)) ## Check that it was created successfully
dir_plots <- here::here(
    "plots",
    "10_Clinical_Gene_Set_Enrichment"
)
dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)
stopifnot(file.exists(dir_plots))

## Load the gene sets from the HumanPilot project
load(here(dir_rdata, "gene_sets_HumanPilot.Rdata"), verbose = TRUE)

## Load the modeling results from the BayesSpace models
bayesSpace_registration_fn <-
    map(k_list, ~ here(
        dir_input,
        paste0(
            "modeling_results_BayesSpace_k",
            sprintf("%02d", .x),
            ".Rdata"
        )
    ))
bayesSpace_registration <-
    lapply(bayesSpace_registration_fn, function(x) {
        get(load(x))
    })



## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
