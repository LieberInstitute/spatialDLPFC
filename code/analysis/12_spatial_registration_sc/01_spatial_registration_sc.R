library("SpatialExperiment")
library("spatialLIBD")
library("here")
library("sessioninfo")

# load sn data
load(file = "/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/DLPFC_snRNAseq/processed-data/sce/sce_DLPFC.Rdata", verbose = TRUE)

## Factor categorical variables used as covariates
colData(sce)$Position <- as.factor(colData(sce)$Position)
colData(sce)$sex <- as.factor(colData(sce)$sex)

## Run spatial registration
sn_hc_registration <- registration_wrapper(sce = sce,
                                           var_registration = "cellType_hc", 
                                           var_sample_id = "Sample",
                                           covars = c("Position", "age", "sex"),
                                           gene_ensembl = "gene_id",
                                           gene_name = "gene_name",
                                           prefix = "")

save(sn_hc_registration, file = here("processed-data","rdata","spe","12_spatial_registration_sc","sn_hc_registration.RDS"))

# sgejobs::job_single('01_spatial_registration_sn', create_shell = TRUE, memory = '25G', command = "Rscript 01_spatial_registration_sn.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

