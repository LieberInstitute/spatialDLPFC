library("SingleCellExperiment")
library("spatialLIBD")
library("jaffelab")
library("here")
library("sessioninfo")

#### load dataset  ####
args <- commandArgs(trailingOnly = TRUE)
dataset <- args[1]
message("Running - ", dataset)

## for testing
# dataset <- "MultiomeBrain-DLPFC"
# dataset <- "UCLA-ASD"

sce_pseudo <- readRDS(file = here(
    "processed-data", "rdata", "spe", "14_spatial_registration_PEC",
    paste0("pseudobulk_", dataset, ".rds")
))

message("\nSCE Dimesions:")
dim(sce_pseudo)

#### Add Dx col ####
dx_data <- read.csv(here(
    "processed-data", "rdata", "spe", "14_spatial_registration_PEC",
    paste0("primaryDiagnosis_", dataset, ".csv")
))

dx_data |> dplyr::count(primaryDiagnosis)

sce_pseudo$primaryDiagnosis <- dx_data$primaryDiagnosis[match(sce_pseudo$individualID, dx_data$individualID)]

table(sce_pseudo$registration_variable, sce_pseudo$primaryDiagnosis)

## exclude any Dx NAs
sce_pseudo <- sce_pseudo[, !is.na(sce_pseudo$primaryDiagnosis)]

# Drop pb samples with < 10 nuc input
min_ncells <- 10
message(
    Sys.time(),
    " dropping ",
    sum(sce_pseudo$ncells < min_ncells),
    " pseudo-bulked samples that are below 'min_ncells'."
)
sce_pseudo <- sce_pseudo[, sce_pseudo$ncells >= min_ncells]


#### Run models ####
dx_split <- rafalib::splitit(sce_pseudo$primaryDiagnosis)
results_enrichment <- purrr::map2(dx_split, names(dx_split), function(i, dx_name) {
    sce_temp <- sce_pseudo[, i]
    message("Dx: ", dx_name, " - ", ncol(sce_temp), " cols")
    var_tab <- table(sce_temp$registration_variable)

    if (any(var_tab == 0)) message("Dropping Empty Levels: ", paste0(names(var_tab)[var_tab == 0], collpase = " "))
    ## Drop Levels
    sce_temp$registration_variable <- droplevels(sce_temp$registration_variable)


    registration_mod <- registration_model(sce_temp)
    block_cor <- registration_block_cor(sce_temp, registration_model = registration_mod)

    gene_name <- "gene_name"
    gene_ensembl <- "featureid"

    registration_stats_enrichment(
        sce_temp,
        block_cor = block_cor,
        gene_ensembl = gene_ensembl,
        gene_name = gene_name
    )
})


## Save results
saveRDS(results_enrichment,
    file = here(
        "processed-data", "rdata", "spe", "14_spatial_registration_PEC",
        paste0("registration_stats_Dx_", dataset, ".rds")
    )
)

# # Run with job_loop
# datasets <- c("CMC",
#               "DevBrain",
#               "IsoHuB",
#               "SZBDMulti",
#               "UCLA-ASD",
#               "Urban-DLPFC")

# sgejobs::job_single('05_compute_registration_Dx', create_shell = TRUE, queue= 'bluejay', memory = '10G', command = "Rscript 05_compute_registration_Dx.R")

#
# job_loop(
#   loops = list(dataset = datasets),
#   name = '05_compute_registration_Dx',
#   create_shell = TRUE
# )

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
