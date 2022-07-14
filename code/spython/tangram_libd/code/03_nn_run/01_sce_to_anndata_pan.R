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

marker_path_in = here("tangram_libd", "raw-data", "03_nn_run", "marker_stats_pan.v2.Rdata")
marker_path_out = here("tangram_libd", "processed-data", "03_nn_run", "marker_genes.Rdata")
  
print('Loading objects...')

#  snRNAseq and spatial objects, respectively
load(here("tangram_libd", "raw-data", "03_nn_run", "sce_pan.v2.Rdata"), verbose = TRUE)
visium_DLPFC = spatialLIBD::fetch_data("spe")

#  Append 'spatialCoords' slot to 'colData', since in
#  conversion we're treating the spatialExperiment object as if it is a
#  singleCellExperiment, which doesn't have that additional slot. Note that
#  in this case, the 'spatialData' slot only contains duplicate information,
#  which is why we don't append it
colData(visium_DLPFC) = cbind(
    colData(visium_DLPFC), spatialCoords(visium_DLPFC)
)

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

load(marker_path_in, verbose = TRUE)

print('Determining and writing markers...')

n_genes <- 25
marker_stats <- marker_stats %>% 
  mutate(Marker = case_when(rank_ratio <= n_genes & !(Symbol %in% c("MRC1","LINC00278"))  ~ "Marker QC",
                            rank_ratio <= n_genes ~ 'Marker top25',
                            TRUE ~ 'Non-Marker'))

marker_stats <- data.frame(marker_stats[marker_stats$Marker == 'Marker QC', ])

#  All the marker genes are present in the spatial data, as required
all(marker_stats$gene %in% rowData(visium_DLPFC)$gene_id)

#  Save entire marker object
save(marker_stats, file=marker_path_out)

#  Save just gene names for use in python
writeLines(
    marker_stats$gene,
    con = here("tangram_libd", "processed-data", "03_nn_run", "pan_markers.txt")
)

###############################################################################
#  Filter using expression cutoffs
###############################################################################

seed <- 20210324

# Find expression cutoff for each dataset and filter out genes below that cutoff
#pdf(expr_plot_out)
#counts_sce = as.matrix(assays(sce_pan)$counts)
#cutoff_sce <- max(expression_cutoff(counts_sce, seed = seed))
#sce_pan1 = sce_pan[rowMeans(counts_sce) > cutoff_sce, ]
#dev_off()

counts_visium = as.matrix(assays(visium_DLPFC)$counts)
cutoff_visium <- expression_cutoff(counts_visium, seed = seed)
visium_DLPFC_H = visium_DLPFC[rowMeans(counts_visium) > max(cutoff_visium), ]
visium_DLPFC_L = visium_DLPFC[rowMeans(counts_visium) > 0.02, ]
visium_DLPFC_L1 = visium_DLPFC[rowMeans(counts_visium) > 0.03, ]

VisHigh_overlaps = intersect(marker_stats$gene,rowData(visium_DLPFC_H)$gene_id)
VisLow = intersect(marker_stats$gene,rowData(visium_DLPFC_L)$gene_id)
VisLow1 = intersect(marker_stats$gene,rowData(visium_DLPFC_L1)$gene_id)
VisLow_overlaps = setdiff(VisLow,VisLow1)

writeLines(
  VisLow_overlaps,
  con = here("tangram_libd", "processed-data", "03_nn_run", "VisLow_overlaps.txt")
)
writeLines(
  VisHigh_overlaps,
  con = here("tangram_libd", "processed-data", "03_nn_run", "VisHigh_overlaps.txt")
)

###############################################################################
#  Write the sample names to a file for use in the python script
###############################################################################

print('Writing sample names for use in python...')
ids = as.character(unique(colData(visium_DLPFC)$sample_id))
writeLines(
    ids,
    con = here("tangram_libd", "processed-data", "03_nn_run", "brain_samples.txt")
)

print('Done.')

session_info()
