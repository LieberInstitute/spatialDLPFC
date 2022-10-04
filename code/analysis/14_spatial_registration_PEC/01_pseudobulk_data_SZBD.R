library("zellkonverter")
library("SingleCellExperiment")
library("spatialLIBD")
library("here")
library("sessioninfo")

#### load dataset  ####

dataset <- "SZBDMulti"
message("\n#### Running: ", dataset, " ####")
## Read data from multiple files
filepath_1 <- here("raw-data", "psychENCODE", "version2", dataset, paste0(dataset, "-Seq_annotated_part1.h5ad"))
filepath_2 <- here("raw-data", "psychENCODE", "version2", dataset, paste0(dataset, "-Seq_annotated_part2.h5ad"))
stopifnot(all(file.exists(c(filepath_1, filepath_2))))

message(Sys.time(), " - Reading data from: ", filepath_1)
sce <- readH5AD(file = filepath_1)

rowData(sce)$varm <- NULL

message("\nSCE Dimesions:")
dim(sce)

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

counts(sce) <- assays(sce)$X # change to normalized counts
# counts(sce)[counts(sce) != 0] <- (2^counts(sce)[counts(sce) != 0])-1 # Replace just non-zero values
counts(sce)@x <- 2^(counts(sce)@x) - 1 ## remove log2(counts + 1)

sce_subset <- sce[,1:1000]
#### Pseudobulk ####
sce_pseudo <- registration_pseudobulk(sce, var_registration = "cellType", var_sample_id = "sampleID", covars = NULL)

message("\nSCE Pseudobulk Dimesions:")
dim(sce_pseudo)

#### Load sce part 2 ####
message(Sys.time(), " - Reading data from: ", filepath_2)
sce <- readH5AD(file = filepath_2)

rowData(sce)$varm <- NULL

message("\nSCE Dimesions:")
dim(sce)

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

counts(sce) <- assays(sce)$X # change to normalized counts
# counts(sce)[counts(sce) != 0] <- (2^counts(sce)[counts(sce) != 0])-1 # Replace just non-zero values
counts(sce)@x <- 2^(counts(sce)@x) - 1 ## remove log2(counts + 1)

sce_subset2 <- sce[,1:1000]
#### Pseudobulk ####
sce_pseudo2 <- registration_pseudobulk(sce, var_registration = "cellType", var_sample_id = "sampleID", covars = NULL)


save(sce_subset, sce_pseudo,sce_subset2, sce_pseudo2,
     file = here(
       "processed-data", "rdata", "spe", "14_spatial_registration_PEC",
       "sce_subsets_SZBD.Rdata")
)

message(Sys.time(), " Combine two parts")

sce_pseudo2 <- sce_pseudo2[rownames(sce_pseudo),]
sce_pseudo <- cbind(sce_pseudo, sce_pseudo2)

message("\n Full SCE Pseudobulk Dimesions:")
dim(sce_pseudo)

## Save results
saveRDS(sce_pseudo,
    file = here(
        "processed-data", "rdata", "spe", "14_spatial_registration_PEC",
        paste0("pseudobulk_", dataset, ".rds")
    )
)

# sgejobs::job_single('01_pseudobulk_data_SZBD', create_shell = TRUE, memory = '200G', command = "Rscript 01_pseudobulk_data_SZBD.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
