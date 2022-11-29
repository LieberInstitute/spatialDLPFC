library("SpatialExperiment")
library("spatialLIBD")
library("here")
library("sessioninfo")

data_dir <- here("processed-data", "rdata", "spe", "12_spatial_registration_sn")

#### Load sn data & exclude drop cells ####
load(file = here(data_dir, "sce_DLPFC.Rdata"), verbose = TRUE)
sce <- sce[, sce$cellType_hc != "drop"]
sce$cellType_hc <- droplevels(sce$cellType_hc)

table(sce$cellType_hc)

## Factor categorical variables used as covariates
colData(sce)$Position <- as.factor(colData(sce)$Position)
colData(sce)$sex <- as.factor(colData(sce)$sex)

## Use all unique ensembl IDs as rownames
rownames(sce) <- rowData(sce)$gene_id

## Run spatial registration
sn_hc_registration <- registration_wrapper(
    sce = sce,
    var_registration = "cellType_hc",
    var_sample_id = "Sample",
    covars = c("Position", "age", "sex"),
    gene_ensembl = "gene_id",
    gene_name = "gene_name"
)

save(sn_hc_registration, file = here(data_dir, "sn_hc_registration.Rdata"))

# sgejobs::job_single('01_spatial_registration_sn', create_shell = TRUE, memory = '25G', command = "Rscript 01_spatial_registration_sn.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
