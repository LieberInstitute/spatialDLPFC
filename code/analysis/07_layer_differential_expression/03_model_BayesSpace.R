# library(sgejobs)
# sgejobs::job_loop(
#    loops = list(spetype = c(
#        "wholegenome", "targeted"
#     )),
#     name = "03_model_pathology",
#     create_shell = TRUE,
#     queue = "bluejay",
#     memory = "15G")
# To execute the script builder, use: sh 03_model_pathology.sh

# Required libraries
library("getopt")

## Specify parameters
spec <- matrix(c(
    "spetype", "s", 2, "character", "SPE spetype: wholegenome or targeted",
    "help", "h", 0, "logical", "Display help"
), byrow = TRUE, ncol = 5)
opt <- getopt(spec = spec)

## if help was asked for print a friendly message
## and exit with a non-zero error code
if (!is.null(opt$help)) {
    cat(getopt(spec, usage = TRUE))
    q(status = 1)
}

## For testing
if (FALSE) {
    opt <- list(spetype = "wholegenome")
}


library("here")
library("sessioninfo")
library("spatialLIBD")
stopifnot(packageVersion("spatialLIBD") >= "1.9.19")


## output directory
dir_rdata <- here::here("processed-data", "11_grey_matter_only", opt$spetype)
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)
stopifnot(file.exists(dir_rdata)) ## Check that it was created successfully
dir_plots <- here::here("plots", "11_grey_matter_only", opt$spetype)
dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)
stopifnot(file.exists(dir_plots))

## load spe data
sce_pseudo <-
    readRDS(
        file.path(
            dir_rdata,
            paste0("sce_pseudo_pathology_", opt$spetype, ".rds")
        )
    )

## To avoid having to change parameters later on
sce_pseudo$registration_variable <- sce_pseudo$path_groups
sce_pseudo$registration_sample_id <- sce_pseudo$sample_id

## Set arguments used in spatialLIBD::registration_wrapper()
covars <- NULL
gene_ensembl <- "gene_id"
gene_name <- "gene_name"
suffix <- "noWM"

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
        "Visium_IF_AD_modeling_results.Rdata"
    )
)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
