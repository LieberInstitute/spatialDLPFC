suppressPackageStartupMessages(library('basilisk'))
suppressPackageStartupMessages(library('scRNAseq'))
suppressPackageStartupMessages(library('SingleCellExperiment'))
suppressPackageStartupMessages(library('zellkonverter'))
suppressPackageStartupMessages(library('data.table'))
suppressPackageStartupMessages(library('tidyverse'))
suppressPackageStartupMessages(library("DeconvoBuddies"))
suppressPackageStartupMessages(library("jaffelab"))
suppressPackageStartupMessages(library('sessioninfo'))
library("here")

#  Path to write the python AnnData object
dir.create(file.path(here::here(), "tangram_libd/processed-data/03_nn_run"), showWarnings = FALSE)
visium_out = file.path(here::here(), "tangram_libd/processed-data/03_nn_run/visium_dlpfc.h5ad")
sc_out = file.path(here::here(), "tangram_libd/processed-data/03_nn_run/sce_pan.h5ad")

print('Loading objects...')

#  snRNAseq and spatial objects, respectively
load(file.path(here::here(), "tangram_libd/raw-data/01_prepare_tangram/sce_pan.Rdata"))
load(file.path(here::here(), "tangram_libd/raw-data/01_prepare_tangram/Human_DLPFC_Visium_processedData_sce_scran_spatialLIBD.Rdata"))

#  Load Louise's marker stats for pan-brain
load(file.path(here::here(), "tangram_libd/raw-data/01_prepare_tangram/marker_stats_pan.Rdata"))

gc()

###############################################################################
#  Convert SCE R objects to AnnData python objects, as a preprocessing step to
#  running tangram
###############################################################################

write_anndata = function(sce, out_path) {
    invisible(
        basiliskRun(
            fun = function(sce, filename) {
                library('zellkonverter')
                library('reticulate')
                
                # Convert SCE to AnnData:
                adata <- SCE2AnnData(sce)
                
                #  Write AnnData object to disk
                adata$write(filename = filename)
                
                return()
            },
            env = zellkonverter:::anndata_env,
            sce = sce,
            filename = out_path
        )
    )
}

print('Writing AnnDatas...')
write_anndata(sce_pan, sc_out)
write_anndata(sce, visium_out)
gc()

###############################################################################
#  Find marker genes, starting with Louise's marker stats
###############################################################################

print('Determining and writing markers...')
marker_stats_filter <- marker_stats %>%
    filter(gene %in% rownames(sce)) %>%
    arrange(rank_ratio) %>%
    group_by(cellType.target) %>%
    slice(1:25)
    
marker_genes <- marker_stats_filter$Symbol

writeLines(
    marker_genes,
    con = file.path(here::here(), "tangram_libd/processed-data/03_nn_run/pan_markers.txt")
)

###############################################################################
#  Write the sample names to a file for use in the python script
###############################################################################

print('Writing sample names for use in python...')
ids = as.character(unique(colData(sce)$sample_name))
writeLines(
    ids,
    con = here("tangram_libd", "processed-data", "03_nn_run", "brain_samples.txt")
)

print('Done.')

session_info()
