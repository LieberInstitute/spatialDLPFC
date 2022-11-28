library("SingleCellExperiment")
library("spatialLIBD")
library("zellkonverter")
library("here")
library("sessioninfo")

#### load dataset  ####
dataset <- "UCLA-ASD"

## what assays are in the sce files?
list.files(here("raw-data", "psychENCODE", "version2", "UCLA-ASD"))
# [1] "Counts"                                  "SYNAPSE_METADATA_MANIFEST.tsv"
# [3] "UCLA-ASD_annotated.h5ad"                 "UCLA-ASD_cell_metadata.tsv"
# [5] "UCLA-ASD_sample_provenance.tsv"          "UCLA-ASD-snRNAseq_annotated.h5ad"
# [7] "UCLA-ASD-snRNAseq_cell_metadata.tsv"     "UCLA-ASD-snRNAseq_sample_provenance.tsv"

## file syntax matches other datasets "DATASET-snRNAseq_annotated.h5ad
filepath_1 <- here("raw-data", "psychENCODE", "version2", "UCLA-ASD", "UCLA-ASD-snRNAseq_annotated.h5ad")
message(Sys.time(), " - Reading data from: ", filepath_1)
sce1 <- readH5AD(file = filepath_1)
names(assays(sce1))
# [1] "X"

## Not empty??
assays(sce1)$X[1:5, 1:5]
# 5 x 5 sparse Matrix of class "dgCMatrix"
#                        18BW-AGCTTCCTCTCAACGA 18BW-AGGAAATCATTCACCC 18BW-TATTGCTAGTAAATGC 18BW-CGAAGGAGTCCTTGTC
# AL627309.1             .                             .              .                     .
# AL627309.3             .                             .              .                     .
# AL627309.5             0.8294494                     .              1.034114              2.111415
# AP006222.2             .                             .              .                     .
# AL669831.2             .                             .              .                     .
#                         18BW-TTCTAACAGACTGTTC
# AL627309.1              .
# AL627309.3              .
# AL627309.5              1.053666
# AP006222.2              .
# AL669831.2              .

all(assays(sce1)$X[, 1:100] == 0)
# [1] FALSE ## cool!

## Other h5ad file (not sure relation to the first, pattern in IsoHub, CMC too)
filepath_2 <- here("raw-data", "psychENCODE", "version2", "UCLA-ASD", "UCLA-ASD_annotated.h5ad")
message(Sys.time(), " - Reading data from: ", filepath_2)
sce2 <- readH5AD(file = filepath_2)
names(assays(sce2))



datasets <- c("DevBrain", "IsoHuB", "CMC", "UCLA-ASD", "SZBDMulti")
names(datasets) <- datasets
purrr::map(datasets, ~ list.files(here("raw-data", "psychENCODE", "version2", .x), pattern = ".h5ad"))

#### Load pseudobulk data ####
sce_pseudo <- readRDS(file = here(
    "processed-data", "rdata", "spe", "14_spatial_registration_PEC",
    paste0("pseudobulk_", dataset, ".rds")
))

sce_pseudo_dev <- readRDS(file = here(
    "processed-data", "rdata", "spe", "14_spatial_registration_PEC",
    paste0("pseudobulk_", "DevBrain", ".rds")
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

(var_tab2 <- table(sce_pseudo$registration_variable))
# Astro Chandelier       Endo    L2.3.IT      L4.IT    L5.6.NP      L5.ET      L5.IT      L6.CT      L6.IT
# 41         40         37         41         41         35         29         41         40         41
# L6.IT.Car3        L6b      Lamp5 Lamp5.Lhx6  Micro.PVM      Oligo        OPC       Pax6      Pvalb       Sncg
# 39         38         40         38         40         41         41         30         41         40
# Sst  Sst.Chodl        Vip       VLMC
# 41          2         41         41

#### Run models ####
registration_mod <- registration_model(sce_pseudo)
# block_cor <- registration_block_cor(sce_pseudo, registration_model = registration_mod)
block_cor <- NaN

gene_name <- "gene_name"
gene_ensembl <- "featureid"

results_enrichment <-
    registration_stats_enrichment(
        sce_pseudo,
        block_cor = block_cor,
        gene_ensembl = gene_ensembl,
        gene_name = gene_name
    )

# 2022-09-20 12:31:24 computing enrichment statistics
# 2022-09-20 12:32:23 extract and reformat enrichment results
# Error in h(simpleError(msg, call)) :
#   error in evaluating the argument 'i' in selecting a method for function '[[': object 'gene_ensembl' not found
# In addition: There were 24 warnings (use warnings() to see them)

## 24 identical messages
# Warning messages:
#   1: In fitFDist(var, df1 = df, covariate = covariate) :
#   More than half of residual variances are exactly zero: eBayes unreliable

## Debug registration_stats_enrichment
covars <- NULL
var_registration <- "registration_variable"
var_sample_id <- "registration_sample_id"
gene_ensembl <- NULL
gene_name <- NULL

cluster_idx <- split(seq(along = sce_pseudo[[var_registration]]), sce_pseudo[[var_registration]])
sapply(cluster_idx, length)

table(sce_pseudo[[var_registration]])
# Astro Chandelier       Endo    L2.3.IT      L4.IT    L5.6.NP      L5.ET      L5.IT      L6.CT      L6.IT
# 41         41         41         41         41         40         41         41         41         41
# L6.IT.Car3        L6b      Lamp5 Lamp5.Lhx6  Micro.PVM      Oligo        OPC       Pax6      Pvalb       Sncg
# 41         41         41         41         41         41         41         41         41         41
# Sst  Sst.Chodl        Vip       VLMC
# 41         34         41         41

message(Sys.time(), " computing enrichment statistics")
eb0_list_cluster <- lapply(cluster_idx, function(x) {
    res <- rep(0, ncol(sce_pseudo))
    res[x] <- 1
    if (!is.null(covars)) {
        res_formula <-
            eval(str2expression(paste(
                "~", "res", "+", paste(covars, collapse = " + ")
            )))
    } else {
        res_formula <- eval(str2expression(paste("~", "res")))
    }
    m <- model.matrix(res_formula, data = colData(sce_pseudo))

    if (is.finite(block_cor)) {
        res <- limma::eBayes(limma::lmFit(
            logcounts(sce_pseudo),
            design = m,
            block = sce_pseudo[[var_sample_id]],
            correlation = block_cor
        ))
    } else {
        res <- limma::eBayes(limma::lmFit(logcounts(sce_pseudo),
            design = m
        ))
    }
    return(res)
})

## Just Astro
res <- rep(0, ncol(sce_pseudo))
res[cluster_idx$Astro] <- 1
res_formula <- eval(str2expression(paste("~", "res")))
m <- model.matrix(res_formula, data = colData(sce_pseudo))
res <- limma::eBayes(limma::lmFit(logcounts(sce_pseudo),
    design = m
))

## work with block cor?
res2 <- limma::eBayes(limma::lmFit(
    logcounts(sce_pseudo),
    design = m,
    block = sce_pseudo[[var_sample_id]],
    correlation = 0.09
))

# Warning message:
#   In fitFDist(var, df1 = df, covariate = covariate) :
#   More than half of residual variances are exactly zero: eBayes unreliable

## Save results
saveRDS(results_enrichment,
    file = here(
        "processed-data", "rdata", "spe", "14_spatial_registration_PEC",
        paste0("registration_stats_", dataset, ".rds")
    )
)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
