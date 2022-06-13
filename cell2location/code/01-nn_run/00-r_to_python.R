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

sce_in = here(
    "tangram_libd", "raw-data", "03_nn_run", "SCE_DLPFC-n3_tran-etal.rda"
)
spe_out = here("cell2location", "processed-data", "01-nn_run", "spe.h5ad")
sce_out = here("cell2location", "processed-data", "01-nn_run", "sce_dlpfc.h5ad")

marker_path_in = here(
    "tangram_libd", "raw-data", "03_nn_run", "marker_stats_pan.v2.Rdata"
)
marker_path_out = here(
    "cell2location", "processed-data", "01-nn_run", "dlpfc_markers.txt"
)

cell_types_to_drop = c('Endo', 'Macrophage', 'Mural', 'Tcell')

###############################################################################
#   Functions
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
#  Convert SCE R objects to AnnData python objects, as a preprocessing step to
#  running cell2location
###############################################################################

dir.create(dirname(spe_out), recursive = TRUE, showWarnings = FALSE)

#  snRNAseq and spatial objects, respectively
print('Loading objects...')
load(sce_in, verbose = TRUE)
spe = spatialLIBD::fetch_data("spe")

#   Drop rare cell types for single-cell data
print('Cell types in single-cell originally:')
levels(sce.dlpfc.tran$cellType)

#   Manually merge finer cell types into the standard broad categories
sce.dlpfc.tran$cellType = as.character(sce.dlpfc.tran$cellType)
sce.dlpfc.tran$cellType[substr(sce.dlpfc.tran$cellType, 1, 5) == "Excit"] = "Excit"
sce.dlpfc.tran$cellType[substr(sce.dlpfc.tran$cellType, 1, 5) == "Inhib"] = "Inhib"

print('Distribution of cells to keep (FALSE) vs. drop (TRUE):')
table(sce.dlpfc.tran$cellType %in% cell_types_to_drop)
sce.dlpfc.tran = sce.dlpfc.tran[, ! (sce.dlpfc.tran$cellType %in% cell_types_to_drop)]
sce.dlpfc.tran$cellType = as.factor(sce.dlpfc.tran$cellType)

#  Append 'spatialCoords' slot to 'colData', since in
#  conversion we're treating the spatialExperiment object as if it is a
#  singleCellExperiment, which doesn't have that additional slot
colData(spe) = cbind(colData(spe), spatialCoords(spe))

print('Writing AnnDatas...')
write_anndata(sce.dlpfc.tran, sce_out)
#write_anndata(spe, spe_out)
gc()

###############################################################################
#  Find marker genes
###############################################################################

print('Determining and writing markers...')

load(marker_path_in, verbose = TRUE)

#   Take top N marker genes for each (non-rare) cell type
n_genes <- 100
marker_stats = marker_stats %>% 
    filter(! cellType.target %in% cell_types_to_drop) %>%
    filter(rank_ratio <= n_genes)


#   All the marker genes are present in the single-cell data (sanity check) and
#   spatial data, as required
all(marker_stats$gene %in% rowData(sce.dlpfc.tran)$gene_id)
all(marker_stats$gene %in% rowData(spe)$gene_id)

#   Write list of markers
writeLines(marker_stats$gene, con = marker_path_out)

session_info()
