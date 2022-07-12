#   Convert IF SpatialExperiment, non-IF SpatialExperiment, and snRNA-seq
#   SingleCellExperiment to AnnDatas. Also, find marker genes that will be
#   shared for both the IF and non-IF analyses.

suppressPackageStartupMessages(library("basilisk"))
suppressPackageStartupMessages(library("SingleCellExperiment"))
suppressPackageStartupMessages(library("SpatialExperiment"))
suppressPackageStartupMessages(library("DeconvoBuddies"))
suppressPackageStartupMessages(library("zellkonverter"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("sessioninfo"))
suppressPackageStartupMessages(library("here"))

#   Number of marker genes to use per cell type
n_markers_per_type <- 25

#  Paths
sce_in = "/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/DLPFC_snRNAseq/processed-data/sce/sce_DLPFC.Rdata"
spe_IF_in = here(
    "processed-data", "rdata", "spe_IF", "01_build_spe_IF", "spe.rds"
)
spe_nonIF_in = here(
    "processed-data", "rdata", "spe", "01_build_spe",
    "spe_filtered_final_with_clusters.Rdata"
)

sce_out = here(
    "processed-data", "spot_deconvo", "01-tangram", "sce.h5ad"
)
spe_IF_out = here(
    "processed-data", "spot_deconvo", "01-tangram", "IF", "spe.h5ad"
)
spe_nonIF_out = here(
    "processed-data", "spot_deconvo", "01-tangram", "nonIF", "spe.h5ad"
)
marker_out = here(
    "processed-data", "spot_deconvo", "01-tangram", "markers.txt"
)
sample_IF_out = here(
    "processed-data", "spot_deconvo", "01-tangram", "IF", "sample_ids.txt"
)
sample_nonIF_out = here(
    "processed-data", "spot_deconvo", "01-tangram", "nonIF", "sample_ids.txt"
)

#  Make sure output directories exist
dir.create(dirname(spe_IF_out), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(spe_nonIF_out), recursive = TRUE, showWarnings = FALSE)

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
            env = zellkonverterAnnDataEnv(),
            sce = sce,
            filename = out_path
        )
    )
}

###############################################################################
#   Convert snRNA-seq and spatial R objects to AnnData python objects, as a
#   preprocessing step to running tangram
###############################################################################

#   Load objects
load(sce_in, verbose=TRUE)
spe_IF = readRDS(spe_IF_in)
load(spe_nonIF_in, verbose = TRUE)
spe_nonIF = spe
rm(spe)
gc()

#   zellkonverter doesn't know how to convert the 'spatialCoords' slot. We'd
#   ultimately like the spatialCoords in the .obsm['spatial'] slot of the
#   resulting AnnDatas, which corresponds to reducedDims(spe)$spatial in R
reducedDims(spe_IF)$spatial = spatialCoords(spe_IF)
reducedDims(spe_nonIF)$spatial = spatialCoords(spe_nonIF)

print('Converting all 3 objects to AnnDatas...')
write_anndata(sce, sce_out)
write_anndata(spe_IF, spe_IF_out)
write_anndata(spe_nonIF, spe_nonIF_out)
gc()

###############################################################################
#  Find marker genes
###############################################################################

#   Drop rare cell types for single-cell data
print('Distribution of cells to keep (FALSE) vs. drop (TRUE):')
table(sce$cellType_broad_hc == "Endo.Mural")
sce = sce[, ! (sce$cellType_broad_hc == "Endo.Mural")]

#   Use Ensembl gene IDs for rownames (not gene symbol)
rownames(sce) = rowData(sce)$gene_id

print('Determining and writing markers...')
marker_stats = get_mean_ratio2(
    sce, cellType_col = 'cellType_broad_hc', assay_name = 'logcounts',
    add_symbol = TRUE
)

#   Take top N marker genes for each (non-rare) cell type
marker_stats = marker_stats %>% 
    filter(
        (rank_ratio <= n_markers_per_type) & 
        !(Symbol %in% c("MRC1","LINC00278"))
    )

markers_scratch = marker_stats$gene

#   It's technically possible to find a single gene that is used as a "marker"
#   for two different cell types via this method. Verify this is not the case,
#   because that would indicate the use of bad markers
stopifnot(
    length(unique(markers_scratch)) == n_markers_per_type * length(unique(marker_stats$cellType.target))
)

#   All the marker genes are present in the single-cell data (sanity check), but
#   one marker is in neither the IF nor non-IF spatial data
stopifnot(all(markers_scratch %in% rowData(sce)$gene_id))

start_len = length(markers_scratch)
markers_scratch = markers_scratch[
    (markers_scratch %in% rowData(spe_IF)$gene_id) &
    (markers_scratch %in% rowData(spe_nonIF)$gene_id)
]
print(
    paste(
        'Dropped', start_len - length(markers_scratch),
        'markers, which were not in the spatial data.'
    )
)

#   Write list of markers
writeLines(markers_scratch, con = marker_out)

#   Write sample names to text files
writeLines(unique(spe_IF$sample_id), con = sample_IF_out)
writeLines(unique(spe_nonIF$sample_id), con = sample_nonIF_out)

session_info()
