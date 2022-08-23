#   Previously, I had used DLPFC-only single-cell data, which included finer
#   and rare cell types. I also used training genes selected by some
#   cell2location functions. With these choices, cell2location performed very
#   poorly. Here, I use pan-brain single-cell data, exclude rare cell types,
#   use broad cell types only, and select marker genes for training using
#   Louise's method.

suppressPackageStartupMessages(library("basilisk"))
suppressPackageStartupMessages(library("SingleCellExperiment"))
suppressPackageStartupMessages(library("zellkonverter"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("DeconvoBuddies"))
suppressPackageStartupMessages(library("jaffelab"))
suppressPackageStartupMessages(library("sessioninfo"))
suppressPackageStartupMessages(library("spatialLIBD"))
suppressPackageStartupMessages(library("here"))

sce_in <- here(
    "tangram_libd", "raw-data", "03_nn_run", "SCE_DLPFC-n3_tran-etal.rda"
)
spe_out <- here("cell2location", "processed-data", "01-nn_run", "spe.h5ad")
sce_out <- here("cell2location", "processed-data", "01-nn_run", "sce_dlpfc.h5ad")

marker_path_out <- here(
    "cell2location", "processed-data", "01-nn_run", "dlpfc_markers.txt"
)

cell_types_to_drop <- c("Endo", "Macrophage", "Mural", "Tcell")

n_markers_per_type <- 100

###############################################################################
#   Functions
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
#  Convert SCE R objects to AnnData python objects, as a preprocessing step to
#  running cell2location
###############################################################################

dir.create(dirname(spe_out), recursive = TRUE, showWarnings = FALSE)

#  snRNAseq and spatial objects, respectively
print("Loading objects...")
load(sce_in, verbose = TRUE)
spe <- spatialLIBD::fetch_data("spe")

print("Cell types in single-cell originally:")
levels(sce.dlpfc.tran$cellType)

#   Manually merge finer cell types into the standard broad categories
sce.dlpfc.tran$cellType <- as.character(sce.dlpfc.tran$cellType)
sce.dlpfc.tran$cellType[substr(sce.dlpfc.tran$cellType, 1, 5) == "Excit"] <- "Excit"
sce.dlpfc.tran$cellType[substr(sce.dlpfc.tran$cellType, 1, 5) == "Inhib"] <- "Inhib"

#   Drop rare cell types for single-cell data
print("Distribution of cells to keep (FALSE) vs. drop (TRUE):")
table(sce.dlpfc.tran$cellType %in% cell_types_to_drop)
sce.dlpfc.tran <- sce.dlpfc.tran[
    , !(sce.dlpfc.tran$cellType %in% cell_types_to_drop)
]
sce.dlpfc.tran$cellType <- as.factor(sce.dlpfc.tran$cellType)

#   Use Ensembl gene IDs for rownames (not gene symbol)
rownames(sce.dlpfc.tran) <- rowData(sce.dlpfc.tran)$gene_id

#  Append 'spatialCoords' slot to 'colData', since in
#  conversion we're treating the spatialExperiment object as if it is a
#  singleCellExperiment, which doesn't have that additional slot
colData(spe) <- cbind(colData(spe), spatialCoords(spe))

print("Writing AnnDatas...")
write_anndata(sce.dlpfc.tran, sce_out)
# write_anndata(spe, spe_out)
gc()

###############################################################################
#  Find marker genes
###############################################################################

print("Determining and writing markers...")

marker_stats <- get_mean_ratio2(
    sce.dlpfc.tran,
    cellType_col = "cellType", assay_name = "logcounts"
)

#   Take top N marker genes for each (non-rare) cell type
marker_stats <- marker_stats %>%
    filter(rank_ratio <= n_markers_per_type)

markers_scratch <- marker_stats$gene

#   It's technically possible to find a single gene that is used as a "marker"
#   for two different cell types via this method. Verify this is not the case,
#   because that would indicate the use of bad markers
stopifnot(
    length(unique(markers_scratch)) == n_markers_per_type * length(unique(marker_stats$cellType.target))
)

#   All the marker genes are present in the single-cell data (sanity check) and
#   spatial data, as required
all(markers_scratch %in% rowData(sce.dlpfc.tran)$gene_id)
all(markers_scratch %in% rowData(spe)$gene_id)

#   Write list of markers
writeLines(markers_scratch, con = marker_path_out)

session_info()
