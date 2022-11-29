# library(sgejobs)
# sgejobs::job_single(
#     name = "02_enrichment_HumanPilot_sets",
#     create_shell = TRUE,
#     queue = "bluejay",
#     memory = "5G"
# )
# To execute the script builder, use: qsub 02_enrichment_HumanPilot_sets.sh

library("here")
library("sessioninfo")

## output directory
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


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
