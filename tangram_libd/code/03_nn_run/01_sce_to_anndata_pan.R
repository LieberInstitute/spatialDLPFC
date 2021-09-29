suppressPackageStartupMessages(library('basilisk'))
suppressPackageStartupMessages(library('SingleCellExperiment'))
suppressPackageStartupMessages(library('zellkonverter'))
suppressPackageStartupMessages(library('data.table'))
suppressPackageStartupMessages(library('tidyverse'))
suppressPackageStartupMessages(library("DeconvoBuddies"))
suppressPackageStartupMessages(library("jaffelab"))
suppressPackageStartupMessages(library('sessioninfo'))
suppressPackageStartupMessages(library('spatialLIBD'))
library("here")

#  Path to write the python AnnData object
dir.create(here("tangram_libd", "processed-data", "03_nn_run"), showWarnings = FALSE))
visium_out = here("tangram_libd", "processed-data", "03_nn_run", "visium_dlpfc.h5ad")
sc_out = here("tangram_libd", "processed-data", "03_nn_run", "sce_pan.h5ad")
expr_plot_out = here("tangram_libd", "plots", "03_nn_run", "expression_cutoffs.pdf")

print('Loading objects...')

#  snRNAseq and spatial objects, respectively
load(here("tangram_libd", "raw-data", "03_nn_run", "sce_pan.v2.Rdata"))
sce = spatialLIBD::fetch_data("spe")

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
#  Filter using expression cutoffs
###############################################################################

rna.sce <- sce_pan
spatial.seq <- sce

# prepairing for expression_cutoff
seed <- 20210324

# Find expression cutoff for each dataset and filter out genes below that
# cutoff
pdf(expr_plot_out)
counts_rna = as.matrix(assays(rna.sce)$counts)
cutoff_rna <- max(expression_cutoff(counts_rna, seed = seed))
rna.sce = rna.sce[rowMeans(counts_rna) > cutoff_rna, ]

counts_spat = as.matrix(assays(spatial.seq)$counts)
cutoff_spat <- max(expression_cutoff(counts_spat, seed = seed))
spatial.seq = spatial.seq[rowMeans(counts_spat) > cutoff_spat, ]
dev.off()

# Getting overlaps between spatial and scRNAseq data
overlaps <-
    rna.sce[rowData(rna.sce)$Symbol %in% rowData(spatial.seq)$gene_name,]

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
marker_stats <- marker_stats %>% 
  mutate(Marker = case_when(rank_ratio <= n_genes & !(Symbol %in% c("MRC1","LINC00278"))  ~ "Marker QC",
                            rank_ratio <= n_genes ~ 'Marker top25',
                            TRUE ~ 'Non-Marker')) %>%
  left_join(markers_mathys) %>%
  replace_na(list(marker_anno = FALSE))
    
marker_genes <- marker_stats[marker_stats$Marker == 'Marker QC', 'Symbol']

writeLines(
    marker_genes,
    con = here("tangram_libd", "processed-data", "03_nn_run", "pan_markers.txt")
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
