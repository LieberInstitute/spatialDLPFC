library(spatialLIBD)
library(HDF5Array)
library(here)
library(sessioninfo)

in_visium_dir = '/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/code/06_deploy_app/spe_shiny'
in_visium_path = '/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/05_harmony_BayesSpace/04-preprocess_and_harmony/spe_harmony.rds'
in_sn_path = '/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/12_snRNA/sce_CellType_noresiduals.Rds'
pseudo_sn_path = here(
    'processed-data', 'MFA', 'pseudobulk_rds', 'nac_sn_pb.rds'
)
pseudo_visium_path = here(
    'processed-data', 'MFA', 'pseudobulk_rds', 'nac_visium_pb.rds'
)

dir.create(dirname(pseudo_sn_path), showWarnings = FALSE)

################################################################################
#   Prep snRNA-seq data
################################################################################

message(Sys.time(), ' | Loading and prepping datasets')
sce = readRDS(in_sn_path)
reducedDims(sce) = list()
sce = sce[, sce$CellType.Final != 'Neuron_Ambig']
sce$CellType.Final = factor(paste0('NAc_sn_', sce$CellType.Final))

################################################################################
#   Prep Visium data
################################################################################

spe = loadHDF5SummarizedExperiment(in_visium_dir)
spe = as(spe, "SingleCellExperiment") # throws cryptic error otherwise
rownames(spe) = rowData(spe)$gene_id
reducedDims(spe) = list()
spe$spatial_domains = paste0(
    'NAc_visium_', gsub('[ /]', '_', spe$spatial_domains)
)
spe$spatial_domains[spe$spatial_domains == 'NAc_visium_Endothelial_Ependymal'] = 'NAc_visium_Endo_Ependymal'
spe$spatial_domains = factor(spe$spatial_domains)
stopifnot(setequal(spe$donor, sce$Brain_ID))

#   counts assay was dropped for the shiny app; read it back in from a different
#   object
spe_raw = readRDS(in_visium_path)
stopifnot(identical(rownames(spe), rownames(spe_raw)))
stopifnot(all(colnames(spe) %in% colnames(spe_raw)))
spe_raw = spe_raw[, colnames(spe)]
assays(spe) = list(counts = assays(spe_raw)$counts)

################################################################################
#   Pseudobulk
################################################################################

#   Pseudobulk by donor and spatial domain (note this groups multiple capture
#   areas per donor)!
message(Sys.time(), ' | Pseudobulking Visium data')
spe_pb = registration_pseudobulk(
    spe,
    var_registration = 'spatial_domains',
    var_sample_id = 'donor',
    pseudobulk_rds_file = pseudo_visium_path
)

message(Sys.time(), ' | Pseudobulking snRNA-seq data')
sce_pb = registration_pseudobulk(
    sce,
    var_registration = 'CellType.Final',
    var_sample_id = 'Brain_ID',
    pseudobulk_rds_file = pseudo_sn_path
)

session_info()
