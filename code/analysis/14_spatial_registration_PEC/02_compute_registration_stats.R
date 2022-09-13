library('zellkonverter')
library("SingleCellExperiment")
library("spatialLIBD")
library("here")
library("sessioninfo")


#### load dataset  ####
args <- commandArgs(trailingOnly = TRUE)
dataset <- args[1]
message("Running -", dataset)
filepath = here("raw-data","psychENCODE","version2",dataset, paste0(dataset, "-snRNAseq_annotated.h5ad"))
stopifnot(file.exists(filepath))

message(Sys.time(), " - Reading data from: ", filepath)
sce <- readH5AD(file = here("raw-data/psychENCODE/version2/DevBrain/DevBrain-snRNAseq_annotated.h5ad"))

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
names(assays(sce)) <- "counts" ## this is wrong --- try for now

registration_stats <- registration_wrapper(sce,
                                     var_registration = "cellType",
                                     var_sample_id = "individualID",
                                     covars = NULL,
                                     gene_ensembl = "featureid",
                                     gene_name = "gene_name",
                                     prefix = "",
                                     min_ncells = 10)

## Save results
saveRDS(registration_stats,
        file = here("processed-data","rdata","spe","14_spatial_registration_PEC",
                    paste0("registration_stats_",dataset,".rds")))

# sgejobs::job_single('01_compute_registration_stats_DevBrain', create_shell = TRUE, memory = '25G', command = "Rscript 01_compute_registration_stats.R DevBrain")
# sgejobs::job_single('01_compute_registration_stats_SZBD', create_shell = TRUE, memory = '25G', command = "Rscript 01_compute_registration_stats.R SZBD")
# sgejobs::job_single('01_compute_registration_stats_CMC', create_shell = TRUE, memory = '25G', command = "Rscript 01_compute_registration_stats.R CMC")
# sgejobs::job_single('01_compute_registration_stats_IsoHUB', create_shell = TRUE, memory = '25G', command = "Rscript 01_compute_registration_stats.R IsoHUB")
# sgejobs::job_single('01_compute_registration_stats_UCLA', create_shell = TRUE, memory = '25G', command = "Rscript 01_compute_registration_stats.R UCLA")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
