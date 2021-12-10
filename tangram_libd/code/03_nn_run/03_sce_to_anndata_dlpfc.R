#  In 01_sce_to_anndata_pan.R, we converted the spatial/Visium
#  spatialExperiment to a python AnnData, so here we only need to
#  convert the snRNA-seq singleCellExperiment to an AnnData.
#
#  As before, we use marker genes from Louise, without doing additional
#  filtering based on an expression_cutoff determined for either the spatial
#  or sn object. Therefore, we simply use the marker genes as determined with
#  the pan-brain sn object

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

#  Paths
out_dir_processed = here("tangram_libd", "processed-data", "03_nn_run")
out_dir_plots = here("tangram_libd", "plots", "03_nn_run")

marker_path_in = file.path(out_dir_processed, "marker_genes.Rdata")
sce_in = here("tangram_libd", "raw-data", "03_nn_run", "SCE_DLPFC-n3_tran-etal.rda")

sce_out = file.path(out_dir_processed, "sce_DLPFC.h5ad")
expr_plot_out = file.path(out_dir_plots, "expression_cutoffs_DLPFC.pdf")

#  Make sure output directories exist
dir.create(out_dir_processed, showWarnings = FALSE)
dir.create(out_dir_plots, showWarnings = FALSE)

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
#  Convert snRNA-seq R objects to AnnData python object, as a preprocessing
#  step to running tangram
###############################################################################

load(sce_in, verbose=TRUE)

print('Converting snRNAseq object to AnnData...')
write_anndata(sce.dlpfc.tran, sce_out)
gc()

###############################################################################
#  Perform sanity check on marker genes
###############################################################################

load(marker_path_in, verbose = TRUE)
stopifnot(all(marker_stats$gene %in% rowData(sce.dlpfc.tran)$gene_id))

###############################################################################
#  Examine expression cutoff for the snRNA-seq data (for visualization
#  purposes! We aren't filtering out any marker genes)
###############################################################################

seed <- 20211210

counts_sce = as.matrix(assays(sce.dlpfc.tran)$counts)

#  Plot the determined expression cutoffs
pdf(expr_plot_out)
cutoffs <- expression_cutoff(counts_sce, seed = seed)
dev.off()

#  Look at how many marker genes would be kept based on using the higher or
#  lower cutoff
counts_sce_markers = counts_sce[rowData(sce.dlpfc.tran)$gene_id %in% marker_stats$gene,]

print('Marker genes kept when filtering larger expression_cutoff for snRNAseq object:')
table(rowMeans(counts_sce_markers) > max(cutoffs))

print('Marker genes kept when filtering smaller expression_cutoff for snRNAseq object:')
table(rowMeans(counts_sce_markers) > min(cutoffs))

print('Done all tasks.')
session_info()
