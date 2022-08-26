#   Convert IF SpatialExperiment, non-IF SpatialExperiment, and snRNA-seq
#   SingleCellExperiment to AnnDatas. Save a copy of the modified SCE as an R
#   object as well. This is a processing step shared by, and prior to, the
#   different deconvolution softwares (tangram, cell2location, SPOTlight).
#   Finally, run getMeanRatio2 to rank genes as markers, and save the resulting
#   object.

suppressPackageStartupMessages(library("basilisk"))
suppressPackageStartupMessages(library("SingleCellExperiment"))
suppressPackageStartupMessages(library("SpatialExperiment"))
suppressPackageStartupMessages(library("DeconvoBuddies"))
suppressPackageStartupMessages(library("zellkonverter"))
suppressPackageStartupMessages(library("sessioninfo"))
suppressPackageStartupMessages(library("here"))

#  Paths
sce_in <- "/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/DLPFC_snRNAseq/processed-data/sce/sce_DLPFC.Rdata"
spe_IF_in <- here(
    "processed-data", "rdata", "spe_IF", "01_build_spe_IF", "spe.rds"
)
spe_nonIF_in <- here(
    "processed-data", "rdata", "spe", "01_build_spe",
    "spe_filtered_final_with_clusters.Rdata"
)

sce_out <- here(
    "processed-data", "spot_deconvo", "01-tangram", "sce.h5ad"
)
sce_r_out <- here(
    "processed-data", "spot_deconvo", "sce.rds"
)
spe_IF_out <- here(
    "processed-data", "spot_deconvo", "01-tangram", "IF", "spe.h5ad"
)
spe_nonIF_out <- here(
    "processed-data", "spot_deconvo", "01-tangram", "nonIF", "spe.h5ad"
)
sample_IF_out <- here(
    "processed-data", "spot_deconvo", "01-tangram", "IF", "sample_ids.txt"
)
sample_nonIF_out <- here(
    "processed-data", "spot_deconvo", "01-tangram", "nonIF", "sample_ids.txt"
)
marker_object_out <- here(
    "processed-data", "spot_deconvo", "marker_stats.rds"
)

#  Make sure output directories exist
dir.create(dirname(spe_IF_out), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(spe_nonIF_out), recursive = TRUE, showWarnings = FALSE)

###############################################################################
#  Functions
###############################################################################

write_anndata <- function(sce, out_path) {
    invisible(
        basiliskRun(
            fun = function(sce, filename) {
                library("zellkonverter")
                library("reticulate")
                
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
#   Convert snRNA-seq and spatial R objects to AnnData python objects
###############################################################################

#   Load objects
load(sce_in, verbose = TRUE)
spe_IF <- readRDS(spe_IF_in)
load(spe_nonIF_in, verbose = TRUE)
spe_nonIF <- spe
rm(spe)
gc()

#   zellkonverter doesn't know how to convert the 'spatialCoords' slot. We'd
#   ultimately like the spatialCoords in the .obsm['spatial'] slot of the
#   resulting AnnDatas, which corresponds to reducedDims(spe)$spatial in R
reducedDims(spe_IF)$spatial <- spatialCoords(spe_IF)
reducedDims(spe_nonIF)$spatial <- spatialCoords(spe_nonIF)

#   Drop rare cell types for single-cell data
print("Distribution of cells to keep (FALSE) vs. drop (TRUE):")
table(sce$cellType_broad_hc == "EndoMural")
sce <- sce[, !(sce$cellType_broad_hc == "EndoMural")]

#   Use Ensembl gene IDs for rownames (not gene symbol)
rownames(sce) <- rowData(sce)$gene_id

#   Save a copy of the filtered + slightly modified sce as an R object, and
#   convert all objects to Anndatas
saveRDS(sce, sce_r_out)

# print('Converting all 3 objects to AnnDatas...')
# write_anndata(sce, sce_out)
# write_anndata(spe_IF, spe_IF_out)
# write_anndata(spe_nonIF, spe_nonIF_out)
gc()

#   Write sample names to text files
writeLines(unique(spe_IF$sample_id), con = sample_IF_out)
writeLines(unique(spe_nonIF$sample_id), con = sample_nonIF_out)

print('Running getMeanRatio2 to rank genes as markers...')
marker_stats <- get_mean_ratio2(
    sce,
    cellType_col = "cellType_broad_hc", assay_name = "logcounts"
)
saveRDS(marker_stats, marker_object_out)

session_info()
