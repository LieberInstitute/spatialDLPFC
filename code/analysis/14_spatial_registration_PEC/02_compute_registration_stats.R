# library('zellkonverter')
library("SingleCellExperiment")
library("spatialLIBD")
library("here")
library("sessioninfo")

#### load dataset  ####
args <- commandArgs(trailingOnly = TRUE)
dataset <- args[1]
message("Running - ", dataset)

sce_pseudo <- readRDS(file = here(
    "processed-data", "rdata", "spe", "14_spatial_registration_PEC",
    paste0("pseudobulk_", dataset, ".rds")
))

message("\nSCE Dimesions:")
dim(sce_pseudo)

min_ncells <- 10
message(
    Sys.time(),
    " dropping ",
    sum(sce_pseudo$ncells < min_ncells),
    " pseudo-bulked samples that are below 'min_ncells'."
)
sce_pseudo <- sce_pseudo[, sce_pseudo$ncells >= min_ncells]

message("Cell Types:")
## must be syntactically valid
(var_tab <- table(sce_pseudo$registration_variable))

if (any(var_tab == 0)) message("Dropping Empty Levels: ", paste0(names(var_tab)[var_tab == 0], collpase = " "))
## Drop Levels
sce_pseudo$registration_variable <- droplevels(sce_pseudo$registration_variable)


#### Run models ####
registration_mod <- registration_model(sce_pseudo)
block_cor <- registration_block_cor(sce_pseudo, registration_model = registration_mod)

gene_name <- "gene_name"
gene_ensembl <- "featureid"

results_enrichment <-
    registration_stats_enrichment(
        sce_pseudo,
        block_cor = block_cor,
        gene_ensembl = gene_ensembl,
        gene_name = gene_name
    )

## Save results
saveRDS(results_enrichment,
    file = here(
        "processed-data", "rdata", "spe", "14_spatial_registration_PEC",
        paste0("registration_stats_", dataset, ".rds")
    )
)

# sgejobs::job_single('02_compute_registration_stats_DevBrain', create_shell = TRUE, memory = '25G', command = "Rscript 02_compute_registration_stats.R DevBrain")
# sgejobs::job_single('02_compute_registration_stats_CMC', create_shell = TRUE, memory = '25G', command = "Rscript 02_compute_registration_stats.R CMC")
# sgejobs::job_single('02_compute_registration_stats_IsoHuB', create_shell = TRUE, memory = '25G', command = "Rscript 02_compute_registration_stats.R IsoHuB")
# sgejobs::job_single('02_compute_registration_stats_UCLA', create_shell = TRUE, memory = '25G', command = "Rscript 02_compute_registration_stats.R UCLA-ASD")

# sgejobs::job_single('02_compute_registration_stats_SZBD', create_shell = TRUE, memory = '25G', command = "Rscript 02_compute_registration_stats.R SZBDMulti")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
