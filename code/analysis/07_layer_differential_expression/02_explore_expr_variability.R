# library(sgejobs)
# sgejobs::job_loop(
#    loops = list(spetype = c(
#        "wholegenome", "targeted"
#     )),
#     name = "02_explore_expr_variability",
#     create_shell = TRUE,
#     queue = "bluejay",
#     memory = "15G")
# To execute the script builder, use: sh 02_explore_expr_variability.sh

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
library("SingleCellExperiment")
library("scater")

## Load pathology colors
source(here("code", "colors_pathology.R"), echo = TRUE, max.deparse.length = 500)

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

## Define variables to use
vars <- c(
    "age",
    "sample_id",
    "path_groups",
    "subject",
    "sex",
    "pmi",
    "APOe"
)

## Plot PCs with different colors
## Each point here is a sample
pdf(file = file.path(dir_plots, paste0("sce_pseudo_pca.pdf")), width = 14, height = 14)
for (var in vars) {
    p <- plotPCA(
        sce_pseudo,
        colour_by = var,
        ncomponents = 12,
        point_size = 1,
        label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(sce_pseudo)$PCA_var_explained
    )
    if (var == "path_groups") {
        p <- p + scale_color_manual("path_groups", values = colors_pathology)
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
pdf(file = file.path(dir_plots, paste0("sce_pseudo_gene_explanatory_vars.pdf")))
plotExplanatoryVariables(vars)
dev.off()

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
