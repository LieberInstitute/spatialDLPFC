# library(sgejobs)
# sgejobs::job_single(
#     name = "02_explore_expr_variability",
#     create_shell = TRUE,
#     queue = "bluejay",
#     memory = "5G",
#     task_num = 28,
#     tc = 10
# )
# To execute the script builder, use: qsub 02_explore_expr_variability.sh

k <- as.numeric(Sys.getenv("SGE_TASK_ID"))

## For testing
if (FALSE) {
    k <- 2
}

library("here")
library("sessioninfo")
library("SingleCellExperiment")
library("scater")

## Load BayesSpace colors
source(here("code", "analysis", "colors_bayesSpace.R"), echo = TRUE, max.deparse.length = 500)
names(colors_bayesSpace) <-
    paste0("Sp", sprintf("%02d", k), "D", sprintf("%02d", as.integer(names(colors_bayesSpace))))

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

## load sce_pseudo data
sce_pseudo <-
    readRDS(
        file.path(
            dir_rdata,
            paste0("sce_pseudo_BayesSpace_k", sprintf("%02d", k), ".rds")
        )
    )

## Define variables to use
vars <- c(
    "age",
    "sample_id",
    "BayesSpace",
    "subject",
    "sex",
    "position"
)

## Plot PCs with different colors
## Each point here is a sample
pdf(file = file.path(dir_plots, paste0("sce_pseudo_PCs_k", sprintf("%02d", k), ".pdf")), width = 14, height = 14)
for (var in vars) {
    p <- plotPCA(
        sce_pseudo,
        colour_by = var,
        ncomponents = 12,
        point_size = 1,
        label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(sce_pseudo)$PCA_var_explained
    )
    if (var == "BayesSpace") {
        p <- p + scale_color_manual("BayesSpace", values = colors_bayesSpace)
    }
    print(p)
}
dev.off()


## Obtain percent of variance explained at the gene level
## using scater::getVarianceExplained()
vars <- getVarianceExplained(sce_pseudo,
    variables = vars
)

## Now visualize the percent of variance explained across all genes
pdf(file = file.path(dir_plots, paste0("sce_pseudo_gene_explanatory_vars_k", sprintf("%02d", k), ".pdf")))
plotExplanatoryVariables(vars)
dev.off()

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
