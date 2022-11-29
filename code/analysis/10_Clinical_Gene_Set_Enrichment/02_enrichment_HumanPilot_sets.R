# library(sgejobs)
# sgejobs::job_single(
#     name = "03_model_BayesSpace",
#     create_shell = TRUE,
#     queue = "bluejay",
#     memory = "5G",
#     task_num = 28,
#     tc = 10
# )
# To execute the script builder, use: qsub 03_model_BayesSpace.sh

k <- as.numeric(Sys.getenv("SGE_TASK_ID"))

## For testing
if (FALSE) {
    k <- 2
}


library("here")
library("sessioninfo")
library("spatialLIBD")

## output directory
dir_rdata <- here::here(
    "processed-data",
    "rdata",
    "spe",
    "07_layer_differential_expression"
)
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)
stopifnot(file.exists(dir_rdata)) ## Check that it was created successfully
dir_plots <- here::here(
    "plots",
    "07_layer_differential_expression"
)
dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)
stopifnot(file.exists(dir_plots))


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
