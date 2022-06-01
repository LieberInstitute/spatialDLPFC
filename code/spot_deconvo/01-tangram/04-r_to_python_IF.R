suppressPackageStartupMessages(library("basilisk"))
suppressPackageStartupMessages(library("SingleCellExperiment"))
suppressPackageStartupMessages(library("SpatialExperiment"))
suppressPackageStartupMessages(library("zellkonverter"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("sessioninfo"))
suppressPackageStartupMessages(library("here"))

#   Number of marker genes to use per cell type
n_genes <- 25

#  Paths
marker_path_in = "/dcl01/lieber/ajaffe/lab/deconvolution_bsp2/data/marker_stats_pan.v2.Rdata"
sce_in = here("") # TODO!
spe_in = here(
    "processed-data", "rdata", "spe_IF", "01_build_spe_IF", "spe.rds"
)

sce_out = here(
    "processed-data", "spot_deconvo", "01-tangram", "IF", "sce.h5ad"
)
spe_out = here(
    "processed-data", "spot_deconvo", "01-tangram", "IF", "spe.h5ad"
)
marker_out = here(
    "processed-data", "spot_deconvo", "01-tangram", "IF", "markers.txt"
)
sample_out = here(
    "processed-data", "spot_deconvo", "01-tangram", "IF", "sample_ids.txt"
)

#  Make sure output directories exist
dir.create(dirname(sce_out), recursive = TRUE, showWarnings = FALSE)

###############################################################################
#  Functions
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
            env = zellkonverterAnnDataEnv,
            sce = sce,
            filename = out_path
        )
    )
}

###############################################################################
#   Convert snRNA-seq and spatial R objects to AnnData python objects, as a
#   preprocessing step to running tangram
###############################################################################

load(sce_in, verbose=TRUE)
spe = readRDS(spe_in)

#   Save imgData without actual images to the metadata slot, which is converted
#   successfully (unlike imgData as-is). We'll only need the hires-related
#   imgData
imgData(spe) = imgData(spe)[imgData(spe)$image_id == 'hires',]
metadata(spe) = list(
    sample_id = imgData(spe)$sample_id,
    image_id = imgData(spe)$image_id,
    scaleFactor = imgData(spe)$scaleFactor
)

#   Add spatialData to colData, since only colData is converted
colData(spe) = cbind(
    colData(spe), spatialCoords(spe)
)

print('Converting both objects to AnnDatas...')
write_anndata(sce, sce_out)
write_anndata(spe, visium_out)
gc()

###############################################################################
#  Find marker genes
###############################################################################

load(marker_path_in, verbose = TRUE)

#   TODO: Have code to decide what cell types we're keeping in the single-cell
#   data and make sure markers found from marker_stats are only for those cell
#   types. Only keep marker genes that are present in both spatial and single-
#   cell data. Then verify a reasonable number are left
print('Determining and writing markers...')
marker_stats = marker_stats %>%
    filter(rank_ratio <= n_genes & !(Symbol %in% c("MRC1","LINC00278")))

#  Save just marker-gene names for use in python
writeLines(
    marker_stats$gene,
    con = marker_out
)

#   TODO: write sample names to text file

session_info()
