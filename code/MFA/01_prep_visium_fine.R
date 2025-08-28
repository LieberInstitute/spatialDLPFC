#   Pseudobulk the spatialDLPFC Visium dataset by donor and spatial domain;
#   pseudobulk the spatialDLPFC snRNA-seq dataset by donor and fine cell type

library(spatialLIBD)
library(HDF5Array)
library(here)
library(sessioninfo)

pseudo_sn_path = here(
    'processed-data', 'MFA', 'pseudobulk_rds', 'sn_fine_pb.rds'
)
pseudo_visium_path = here(
    'processed-data', 'MFA', 'pseudobulk_rds', 'visium_pb.rds'
)

dir.create(dirname(pseudo_sn_path), showWarnings = FALSE)

message(Sys.time(), ' | Loading and prepping datasets')
sce_path = unzip(fetch_data("spatialDLPFC_snRNAseq"), exdir = tempdir())
sce = loadHDF5SummarizedExperiment(file.path(tempdir(), "sce_DLPFC_annotated"))

#   Drop ambiguous cells
sce = sce[, sce$cellType_broad_hc != "Ambiguous"]
stopifnot(!any(is.na(sce$cellType_layer)))

spe = fetch_data('spatialDLPFC_Visium')
spe = as(spe, "SingleCellExperiment") # throws cryptic error otherwise
stopifnot(setequal(spe$subject, sce$BrNum))

#   Use same notation as paper
spe$BayesSpace_harmony_09 = sprintf('Sp09D%02d', spe$BayesSpace_harmony_09)

#   For speed, bring counts into memory
assays(sce)$counts = as(assays(sce)$counts, "dgCMatrix")

#   Pseudobulk by donor and spatial domain (note this groups 3 capture areas per
#   donor)!
message(Sys.time(), ' | Pseudobulking Visium data')
spe_pb = registration_pseudobulk(
    spe,
    var_registration = 'BayesSpace_harmony_09',
    var_sample_id = 'subject',
    pseudobulk_rds_file = pseudo_visium_path
)

message(Sys.time(), ' | Pseudobulking snRNA-seq data')
sce_pb = registration_pseudobulk(
    sce,
    var_registration = 'cellType_layer',
    var_sample_id = 'BrNum',
    pseudobulk_rds_file = pseudo_sn_path
)

session_info()
