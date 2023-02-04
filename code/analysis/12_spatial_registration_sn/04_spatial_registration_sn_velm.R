library("SpatialExperiment")
library("spatialLIBD")
library("here")
library("sessioninfo")

data_dir <- here("processed-data", "rdata", "spe", "12_spatial_registration_sn")

load("/dcl01/ajaffe/data/lab/singleCell/velmeshev2019/analysis_MNT/SCE_asd-velmeshev-etal_MNT.rda",
    verbose = T
) # sce.asd
sce.asd

table(sce.asd$diagnosis, sce.asd$region)
#         ACC   PFC
# ASD     19984 32019
# Control 22409 30147

## filter to just PFC
sce.asd <- sce.asd[, sce.asd$region == "PFC"]

table(sce.asd$cluster, sce.asd$diagnosis)
#                   ASD Control
# AST-FB           1226    1106
# AST-PP           3625    1836
# Endothelial       557     680
# IN-PV            1113    1123
# IN-SST           1315    1414
# IN-SV2C           612     442
# IN-VIP           2000    1357
# L2/3             4925    3701
# L4               2864    2245
# L5/6              976     996
# L5/6-CC          1093    1253
# Microglia        1222    1320
# Neu-mat          1091     963
# Neu-NRGN-I        652     799
# Neu-NRGN-II      1795    2066
# Oligodendrocytes 3656    5892
# OPC              3297    2954


## Factor categorical variables used as covariates
colData(sce.asd)$sex <- as.factor(colData(sce.asd)$sex)

## remove slashes from cell types
sce.asd$cellType <- make.names(sce.asd$cluster)
table(sce.asd$cellType)

# AST.FB           AST.PP      Endothelial            IN.PV           IN.SST          IN.SV2C           IN.VIP
# 2332             5461             1237             2236             2729             1054             3357
# L2.3               L4             L5.6          L5.6.CC        Microglia          Neu.mat       Neu.NRGN.I
# 8626             5109             1972             2346             2542             2054             1451
# Neu.NRGN.II Oligodendrocytes              OPC
# 3861             9548             6251

## Add rowData
rowData(sce.asd) <- DataFrame(gene_id = rownames(sce.asd))

## revert to counts - might not be the best...
counts(sce.asd) <- 2^(logcounts(sce.asd)) - 1 ## remove log2(counts + 1)

## Run spatial registration
message(Sys.time(), "- Run full PFC data ####")
sn_velm_registration <- registration_wrapper(
    sce = sce.asd,
    var_registration = "cellType",
    var_sample_id = "sample",
    covars = c("age", "sex"),
    gene_ensembl = "gene_id",
    gene_name = "gene_name"
)


sn_velm_registration_dx <- purrr::map(c(asd = "ASD", control = "Control"), function(dx) {
    message(Sys.time(), "- Run ", dx, " data ####")
    sce.asd.temp <- sce.asd[, sce.asd$diagnosis == dx]
    sn_velm_registration <- registration_wrapper(
        sce = sce.asd.temp,
        var_registration = "cellType",
        var_sample_id = "sample",
        covars = c("age", "sex"),
        gene_ensembl = "gene_id",
        gene_name = "gene_name"
    )
})

## try combo Dx cell type?
sce.asd$combo <- paste0(sce.asd$cellType, "_", sce.asd$diagnosis)

message(Sys.time(), "- Run full PFC data with combo ####")
sn_velm_registration_combo <- registration_wrapper(
    sce = sce.asd,
    var_registration = "combo",
    var_sample_id = "sample",
    covars = c("age", "sex"),
    gene_ensembl = "gene_id",
    gene_name = "gene_name"
)

message("#### SAVING ####")
velm_registration <- list(all = sn_velm_registration, combo = sn_velm_registration_combo)
velm_registration <- c(velm_registration, sn_velm_registration_dx)

names(velm_registration)
# [1] "all"     "combo"   "asd"     "control"

save(velm_registration, file = here(data_dir, "sn_velm_registration_velm.Rdata"))

# sgejobs::job_single('04_spatial_registration_sn_velm', create_shell = TRUE, memory = '75G', command = "Rscript 04_spatial_registration_sn_velm.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
