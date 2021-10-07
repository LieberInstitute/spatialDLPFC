setwd("/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spython/")

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

#  Path to write the python AnnData object
dir.create(here("tangram_libd", "processed-data", "03_nn_run"), showWarnings = FALSE)
visium_out = here("tangram_libd", "processed-data", "03_nn_run", "visium_DLPFC.h5ad")
sce_out = here("tangram_libd", "processed-data", "03_nn_run", "sce_pan.h5ad")
expr_plot_out = here("tangram_libd", "plots", "03_nn_run", "expression_cutoffs.pdf")
overlaps_out = here("tangram_libd", "processed-data", "03_nn_run", "overlaps.Rda")
  
print('Loading objects...')

#  snRNAseq and spatial objects, respectively
load(here("tangram_libd", "raw-data", "03_nn_run", "sce_pan.v2.Rdata"))
visium_DLPFC = spatialLIBD::fetch_data("spe")

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
      env = zellkonverterAnnDataEnv,
      sce = sce,
      filename = out_path
    )
  )
}

print('Writing AnnDatas...')
write_anndata(sce_pan, sce_out)
write_anndata(visium_DLPFC, visium_out)
gc()

###############################################################################
#  Find marker genes
###############################################################################

marker_stats <- get_mean_ratio2(
  sce_pan,
  cellType_col = "cellType.Broad",
  assay_name = "logcounts",
  add_symbol = TRUE
)


print('Determining and writing markers...')

n_genes <- 25
marker_stats1 <- marker_stats %>% 
  mutate(Marker = case_when(rank_ratio <= n_genes & !(Symbol %in% c("MRC1","LINC00278"))  ~ "Marker QC",
                            rank_ratio <= n_genes ~ 'Marker top25',
                            TRUE ~ 'Non-Marker'))

marker_genes <- data.frame(marker_stats1[marker_stats1$Marker == 'Marker QC', 'Symbol'])

writeLines(
  marker_genes$Symbol,
  con = here("tangram_libd", "processed-data", "03_nn_run", "pan_markers.txt")
)

###############################################################################
#  Filter using expression cutoffs
###############################################################################

seed <- 20210324

# Find expression cutoff for each dataset and filter out genes below that cutoff
pdf(expr_plot_out)
counts_sce = as.matrix(assays(sce_pan)$counts)
cutoff_sce <- max(expression_cutoff(counts_sce, seed = seed))
sce_pan1 = sce_pan[rowMeans(counts_sce) > cutoff_sce, ]

counts_visium = as.matrix(assays(visium_DLPFC)$counts)
cutoff_visium <- max(expression_cutoff(counts_visium, seed = seed))
visium_DLPFC1 = visium_DLPFC[rowMeans(counts_visium) > cutoff_visium, ]
dev.off()

# Getting overlaps between spatial and scRNAseq data
overlaps <- sce_pan1[rowData(sce_pan1)$Symbol %in% rowData(visium_DLPFC1)$gene_name,]
save(overlaps, file = overlaps_out)

overlaps_genes = rowData(overlaps)$Symbol
overlaps_genes = rowData(visium_DLPFC1)$gene_name
common_genes = intersect(marker_genes$Symbol, overlaps_genes)
###############################################################################
#  Write the sample names to a file for use in the python script
###############################################################################

print('Writing sample names for use in python...')
ids = as.character(unique(colData(visium_DLPFC1)$sample_id))
writeLines(
  ids,
  con = here("tangram_libd", "processed-data", "03_nn_run", "brain_samples.txt")
)

print('Done.')

session_info()