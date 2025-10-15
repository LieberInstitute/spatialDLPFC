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
layer_anno = c(
    'Sp09D01_L1', 'Sp09D02_L1', 'Sp09D03_L2', 'Sp09D04_L5', 'Sp09D05_L3',
    'Sp09D06_WM', 'Sp09D07_L6', 'Sp09D08_L4', 'Sp09D09_WM'
)

dir.create(dirname(pseudo_sn_path), showWarnings = FALSE)

################################################################################
#   Prep snRNA-seq data
################################################################################

message(Sys.time(), ' | Loading and prepping datasets')
sce_path = unzip(fetch_data("spatialDLPFC_snRNAseq"), exdir = tempdir())
sce = loadHDF5SummarizedExperiment(file.path(tempdir(), "sce_DLPFC_annotated"))
rownames(sce) = rowData(sce)$gene_id

#   Drop ambiguous cells
sce = sce[
    ,
    (sce$cellType_broad_hc != "Ambiguous") &
    (sce$cellType_layer != "Excit_ambig")
]
stopifnot(!any(is.na(sce$cellType_layer)))

sce$cellType_layer = paste('DLPFC_sn', sce$cellType_layer, sep = '_')

#   For speed, bring counts into memory
assays(sce)$counts = as(assays(sce)$counts, "dgCMatrix")

################################################################################
#   Prep Visium data
################################################################################

spe = fetch_data('spatialDLPFC_Visium')
spe = as(spe, "SingleCellExperiment") # throws cryptic error otherwise
stopifnot(setequal(spe$subject, sce$BrNum))

#   Annotate spatial domains
spe$BayesSpace_harmony_09 = paste(
    'DLPFC_visium', layer_anno[spe$BayesSpace_harmony_09], sep = '_'
)

################################################################################
#   Pseudobulk
################################################################################

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
