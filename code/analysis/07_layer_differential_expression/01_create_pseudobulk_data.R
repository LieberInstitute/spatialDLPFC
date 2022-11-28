# library(sgejobs)
# sgejobs::job_single(
#     name = "01_create_pseudobulk_data",
#     create_shell = TRUE,
#     queue = "bluejay",
#     memory = "20G",
#     task_num = 28,
#     tc = 10
# )
# To execute the script builder, use: qsub 01_create_pseudobulk_data.sh

k <- as.numeric(Sys.getenv("SGE_TASK_ID"))

## For testing
if (FALSE) {
    k <- 2
}

library("here")
library("spatialLIBD")
library("sessioninfo")
library("scater")

## output directory
dir_rdata <- here::here(
    "processed-data",
    "rdata",
    "spe",
    "07_layer_differential_expression"
)
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)
stopifnot(file.exists(dir_rdata)) ## Check that it was created successfully

## load spe data
load(
    here(
        "processed-data",
        "rdata",
        "spe",
        "01_build_spe",
        "spe_filtered_final_with_clusters.Rdata"
    ),
    verbose = TRUE
)

## Convert from character to a factor
spe$BayesSpace <-
    factor(
        paste0("Sp", sprintf("%02d", k), "D", sprintf("%02d", colData(spe)[[paste0("bayesSpace_harmony_", k)]]))
    )

## pseudobulk across a given BayesSpace k
sce_pseudo <-
    registration_pseudobulk(spe,
        var_registration = "BayesSpace",
        var_sample_id = "sample_id",
        min_ncells = 10
    )
dim(sce_pseudo)

## Rename "region" into "position" for consistency with
## https://github.com/LieberInstitute/DLPFC_snRNAseq
sce_pseudo$position <- sce_pseudo$region

## Simplify the colData()  for the pseudo-bulked data
colData(sce_pseudo) <- colData(sce_pseudo)[, sort(c(
    "age",
    "sample_id",
    "BayesSpace",
    "subject",
    "sex",
    "position",
    "diagnosis",
    "ncells"
))]

## Explore the resulting data
options(width = 400)
as.data.frame(colData(sce_pseudo))

## Compute PCs
## Adapted from https://github.com/LieberInstitute/spatialDLPFC/blob/f47daafa19b02e6208c7e0a9bc068367f806206c/code/analysis/09_region_differential_expression/preliminary_analysis.R#L60-L68
pca <- prcomp(t(assays(sce_pseudo)$logcounts))
message(Sys.time(), " % of variance explained for the top 20 PCs:")
metadata(sce_pseudo)
metadata(sce_pseudo) <- list("PCA_var_explained" = jaffelab::getPcaVars(pca)[seq_len(20)])
metadata(sce_pseudo)
pca_pseudo <- pca$x[, seq_len(20)]
colnames(pca_pseudo) <- paste0("PC", sprintf("%02d", seq_len(ncol(pca_pseudo))))
reducedDims(sce_pseudo) <- list(PCA = pca_pseudo)

## Compute some reduced dims
set.seed(20220423)
sce_pseudo <- scater::runMDS(sce_pseudo, ncomponents = 20)
sce_pseudo <- scater::runPCA(sce_pseudo, name = "runPCA")

## We don't want to model the pathology groups as integers / numeric
## so let's double check this
stopifnot(is.factor(sce_pseudo$BayesSpace))

## For the spatialLIBD shiny app
rowData(sce_pseudo)$gene_search <-
    paste0(
        rowData(sce_pseudo)$gene_name,
        "; ",
        rowData(sce_pseudo)$gene_id
    )

## Load pathology colors
## This info is used by spatialLIBD v1.7.18 or newer
source(here("code", "analysis", "colors_bayesSpace.R"), echo = TRUE, max.deparse.length = 500)
names(colors_bayesSpace) <-
    paste0("Sp", sprintf("%02d", k), "D", sprintf("%02d", as.integer(names(colors_bayesSpace))))
sce_pseudo$BayesSpace_colors <- colors_bayesSpace[as.character(sce_pseudo$BayesSpace)]

## save RDS file
saveRDS(
    sce_pseudo,
    file = file.path(
        dir_rdata,
        paste0("sce_pseudo_BayesSpace_k", sprintf("%02d", k), ".rds")
    )
)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
