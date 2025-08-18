#   Pseudobulk the spatialDLPFC Visium dataset by donor and spatial domain;
#   pseudobulk the spatialDLPFC snRNA-seq dataset by donor and fine cell type

library(spatialLIBD)
library(HDF5Array)
library(here)
library(sessioninfo)

pseudo_sn_path = here('processed-data', 'MFA', 'sn_fine_pb.rds')
pseudo_visium_path = here('processed-data', 'MFA', 'visium_pb.rds')

message(Sys.time(), ' | Loading and prepping datasets')
sce_path = unzip(fetch_data("spatialDLPFC_snRNAseq"), exdir = tempdir())
sce = loadHDF5SummarizedExperiment(file.path(tempdir(), "sce_DLPFC_annotated"))

spe = fetch_data('spatialDLPFC_Visium')
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

session_info()
