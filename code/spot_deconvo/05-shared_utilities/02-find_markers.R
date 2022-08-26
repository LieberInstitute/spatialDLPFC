#   Find marker genes that will be shared for both the IF and non-IF analyses.

suppressPackageStartupMessages(library("SingleCellExperiment"))
suppressPackageStartupMessages(library("SpatialExperiment"))
suppressPackageStartupMessages(library("DeconvoBuddies"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("sessioninfo"))
suppressPackageStartupMessages(library("here"))

#   Number of marker genes to use per cell type. Note that cell2location seems
#   to need more markers than the other tools, motivating the exception below
n_markers_per_type <- 25
n_markers_per_type_c2l <- 100

#  Paths
sce_in <- here(
    "processed-data", "spot_deconvo", "sce.rds"
)
spe_IF_in <- here(
    "processed-data", "rdata", "spe_IF", "01_build_spe_IF", "spe.rds"
)
spe_nonIF_in <- here(
    "processed-data", "rdata", "spe", "01_build_spe",
    "spe_filtered_final_with_clusters.Rdata"
)
marker_object_in <- here(
    "processed-data", "spot_deconvo", "marker_stats.rds"
)

marker_out <- here(
    "processed-data", "spot_deconvo", "01-tangram", "markers.txt"
)
marker_c2l_out <- here(
    "processed-data", "spot_deconvo", "01-tangram", "markers.txt"
)

###############################################################################
#  Find marker genes
###############################################################################

#   Load objects
sce = readRDS(sce_in)
spe_IF <- readRDS(spe_IF_in)
load(spe_nonIF_in, verbose = TRUE)
spe_nonIF <- spe
rm(spe)
marker_stats <- readRDS(marker_object_in)
gc()

print("Writing markers...")

#   Take top N marker genes for each cell type
marker_stats <- marker_stats %>%
    filter(rank_ratio <= n_markers_per_type)

markers_scratch <- marker_stats$gene

#   It's technically possible to find a single gene that is used as a "marker"
#   for two different cell types via this method. Verify this is not the case,
#   because that would indicate the use of bad markers
stopifnot(
    length(unique(markers_scratch)) == n_markers_per_type * length(unique(marker_stats$cellType.target))
)

#   All the marker genes are present in the single-cell data (sanity check), but
#   one marker is in neither the IF nor non-IF spatial data
stopifnot(all(markers_scratch %in% rowData(sce)$gene_id))

start_len <- length(markers_scratch)
markers_scratch <- markers_scratch[
    (markers_scratch %in% rowData(spe_IF)$gene_id) &
        (markers_scratch %in% rowData(spe_nonIF)$gene_id)
]
print(
    paste(
        "Dropped", start_len - length(markers_scratch),
        "markers, which were not in the spatial data."
    )
)

#   Write list of markers
writeLines(markers_scratch, con = marker_out)

session_info()
