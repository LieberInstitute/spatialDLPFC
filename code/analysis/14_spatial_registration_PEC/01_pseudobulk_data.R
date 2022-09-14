library('zellkonverter')
library("SingleCellExperiment")
library("spatialLIBD")
library("here")
library("sessioninfo")

#### load dataset  ####
args <- commandArgs(trailingOnly = TRUE)
dataset <- args[1]
message("\n#### Running: ", dataset, " ####")
filepath <- here("raw-data","psychENCODE","version2",dataset, paste0(dataset, "-snRNAseq_annotated.h5ad"))
stopifnot(file.exists(filepath))

message(Sys.time(), " - Reading data from: ", filepath)
sce <- readH5AD(file = filepath)

message("\nSCE Dimesions:")
dim(sce)

message("Cell Types:")
## must be syntactically valid
colData(sce)$cellType <- as.factor(make.names(colData(sce)$subclass))
table(sce$cellType)

# identify annotation/cluster labels
rowData(sce)$gene_name <- rownames(sce) #save gene name as column of rowData
rownames(sce) <- rowData(sce)$featureid # have to make row names of object the ensembl id instead of gene names

## Logcounts 
# default “X” contain the log-normalized counts
message(Sys.time(), " revert to counts")

# message("Sparsity: ", sum(assays(sce)$X == 0)/product(dim(assays(sce)$X)))

counts(sce) <- assays(sce)$X # change to normalized counts
counts(sce)[counts(sce) != 0] <- (2^counts(sce)[counts(sce) != 0])-1 # Replace just non-zero values

message(Sys.time(), " make pseudobulk object")
sce_pseudo <- scuttle::aggregateAcrossCells(
  sce,
  DataFrame(
    registration_variable = sce$cellType,
    registration_sample_id = sce$sampleID),
  use.assay.type = "counts"
)

colnames(sce_pseudo) <-
  paste0(
    sce_pseudo$registration_sample_id,
    "_",
    sce_pseudo$registration_variable
  )

message(Sys.time(), " normalize expression")
logcounts(sce_pseudo) <-
  edgeR::cpm(edgeR::calcNormFactors(sce_pseudo),
             log = TRUE,
             prior.count = 1
  )


## Save results
saveRDS(sce_pseudo,
        file = here("processed-data","rdata","spe","14_spatial_registration_PEC",
                    paste0("pseudobulk_",dataset,".rds")))

# sgejobs::job_single('01_pseudobulk_data_DevBrain', create_shell = TRUE, memory = '25G', command = "Rscript 01_pseudobulk_data.R DevBrain")
# sgejobs::job_single('01_pseudobulk_data_SZBD', create_shell = TRUE, memory = '25G', command = "Rscript 01_pseudobulk_data.R SZBD")
# sgejobs::job_single('01_pseudobulk_data_CMC', create_shell = TRUE, memory = '25G', command = "Rscript 01_pseudobulk_data.R CMC")
# sgejobs::job_single('01_pseudobulk_data_IsoHUB', create_shell = TRUE, memory = '25G', command = "Rscript 01_pseudobulk_data.R IsoHUB")
# sgejobs::job_single('01_pseudobulk_data_UCLA', create_shell = TRUE, memory = '25G', command = "Rscript 01_pseudobulk_data.R UCLA")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
