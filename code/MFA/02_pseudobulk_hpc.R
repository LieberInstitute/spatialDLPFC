#   Pseudobulk the spatialHPC Visium dataset by donor and spatial domain;
#   pseudobulk the spatialHPC snRNA-seq dataset by donor and n=60 cell type

library(spatialLIBD)
library(HDF5Array)
library(here)
library(sessioninfo)

in_sn_path = '/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/snRNAseq_hpc/processed-data/sce/sce_class_final.rda'
in_visium_path = '/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/06_clustering/PRECAST/spe_precast_HE_domain.rda'
pseudo_sn_path = here(
    'processed-data', 'MFA', 'pseudobulk_rds', 'hpc_sn_pb.rds'
)
pseudo_visium_path = here(
    'processed-data', 'MFA', 'pseudobulk_rds', 'hpc_visium_pb.rds'
)

dir.create(dirname(pseudo_sn_path), showWarnings = FALSE)

################################################################################
#   Prep snRNA-seq data
################################################################################

message(Sys.time(), ' | Loading and prepping datasets')
sce = get(load(in_sn_path))
reducedDims(sce) = list()
rownames(sce) = rowData(sce)$gene_id
sce$cell.type = factor(paste0('HPC_sn_', as.character(sce$cell.type)))

################################################################################
#   Prep Visium data
################################################################################

spe = get(load(in_visium_path))
spe = as(spe, "SingleCellExperiment") # throws cryptic error otherwise
reducedDims(spe) = list()
spe$domain = factor(paste0('HPC_visium_', as.character(spe$domain)))
stopifnot(setequal(spe$brnum, sce$brnum))

################################################################################
#   Pseudobulk
################################################################################

#   Pseudobulk by donor and spatial domain (note this groups multiple capture
#   areas per donor)!
message(Sys.time(), ' | Pseudobulking Visium data')
spe_pb = registration_pseudobulk(
    spe,
    var_registration = 'domain',
    var_sample_id = 'brnum',
    pseudobulk_rds_file = pseudo_visium_path
)

message(Sys.time(), ' | Pseudobulking snRNA-seq data')
sce_pb = registration_pseudobulk(
    sce,
    var_registration = 'cell.type',
    var_sample_id = 'brnum',
    pseudobulk_rds_file = pseudo_sn_path
)

session_info()
