library("zellkonverter")
library("SingleCellExperiment")
library("spatialLIBD")
library("jaffelab")
library("here")
library("sessioninfo")

#### load dataset  ####
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
message("input = ", input_file)

dataset <- dirname(input_file)
message("\n#### Running: ", dataset, " ####")
# filepath <- here("raw-data", "psychENCODE", "version2", dataset, paste0(dataset, "-snRNAseq_annotated.h5ad"))

## for v3 data
filepath <- here("raw-data", "psychENCODE", "version3", "scRNAseq_AllenCTHarmonized", input_file)
stopifnot(file.exists(filepath))

message(Sys.time(), " - Reading data from: ", filepath)
sce <- readH5AD(file = filepath)

message("\nSCE Dimesions:")
dim(sce)

print(colnames(colData(sce)))
message("Cell Types:")
## must be syntactically valid
colData(sce)$cellType <- as.factor(make.names(colData(sce)$subclass))
table(sce$cellType)

# identify annotation/cluster labels
rowData(sce)$gene_name <- rownames(sce) # save gene name as column of rowData
rownames(sce) <- rowData(sce)$featureid # have to make row names of object the ensembl id instead of gene names

## Logcounts
# default “X” contain the log-normalized counts
message(Sys.time(), " revert to counts")

## check for all 0s (just first 100 cols for mem)
stopifnot(any(assays(sce)$X[, 1:100] != 0))

counts(sce) <- assays(sce)$X # change to normalized counts
# counts(sce)[counts(sce) != 0] <- (2^counts(sce)[counts(sce) != 0])-1 # Replace just non-zero values
counts(sce)@x <- 2^(counts(sce)@x) - 1 ## remove log2(counts + 1)

#### Pseudobulk ####
message(Sys.time(), " Pseudobulk")
sce_pseudo <- registration_pseudobulk(sce, var_registration = "cellType", var_sample_id = "sampleID", covars = NULL)

message("\nSCE Pseudobulk Dimesions:")
dim(sce_pseudo)

## Save results
saveRDS(sce_pseudo,
    file = here(
        "processed-data", "rdata", "spe", "14_spatial_registration_PEC",
        paste0("pseudobulk_", dataset, ".rds")
    )
)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
