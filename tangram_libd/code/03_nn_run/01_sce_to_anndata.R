library('basilisk')
library('scRNAseq')
library('SingleCellExperiment')
library('zellkonverter')
library('data.table')
library('tidyverse')
library("DeconvoBuddies")
library("here")
library("jaffelab")

#  Path to write the python AnnData object
dir.create(file.path(here::here(), "tangram_libd/processed-data/03_nn_run"), showWarnings = FALSE)
visium_out = file.path(here::here(), "tangram_libd/processed-data/03_nn_run/visium_dlpfc.h5ad")
sc_out = file.path(here::here(), "tangram_libd/processed-data/03_nn_run/sce_dlpfc.h5ad")
expr_plot_out = file.path(here::here(), "tangram_libd/plots/03_nn_run/expression_cutoffs.pdf")

#  The name of the variable in the snRNAseq colData describing cell type 
cell_type_var = 'cellType'

#  snRNAseq and spatial objects, respectively
load(file.path(here::here(), "tangram_libd/raw-data/01_prepare_tangram/regionSpecific_DLPFC-n3_cleaned-combined_SCE_LAH2021.rda"))
load(file.path(here::here(), "tangram_libd/raw-data/01_prepare_tangram/Human_DLPFC_Visium_processedData_sce_scran_spatialLIBD.Rdata"))

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

write_anndata(sce.dlpfc, sc_out)
write_anndata(sce, visium_out)

rna.sce <- sce.dlpfc
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
    rna.sce[rowData(rna.sce)$gene_name %in% rowData(spatial.seq)$gene_name,]

# sce:  SingleCellExperiment object
# n:    number of marker genes to determine (approximately)
# combine: (logical) merge all inhibitory and all excitatory cell types
#       together?
#
# prints number of unique marker genes found and the median mean-ratio of
# those unique marker genes; returns the unique markers
get_markers = function(sce, n, combine=TRUE) {
    sce_copy = sce
    if (combine) {
        #  Merge together all inhibitory and all excitatory "cell types"
        colData(sce_copy)[,cell_type_var] = as.character(colData(sce_copy)[,cell_type_var])
        colData(sce_copy)[substr(colData(sce_copy)[,cell_type_var], 1, 5) == 'Excit', cell_type_var] = 'Excit'
        colData(sce_copy)[substr(colData(sce_copy)[,cell_type_var], 1, 5) == 'Inhib', cell_type_var] = 'Inhib'
        colData(sce_copy)[,cell_type_var] = as.factor(colData(sce_copy)[,cell_type_var])
    }
    
    mean_ratio = 
        get_mean_ratio2(
            sce_copy,
            cellType_col = cell_type_var,
            assay_name = "counts",
            add_symbol = TRUE
        )
    
    #  We are taking the top x markers per cell type, and totaling as close to
    #  n as possible (rounding down if n isn't divisible by x)
    n_per_group = n %/% nlevels(mean_ratio$cellType.target)
    n_actual = n_per_group * nlevels(mean_ratio$cellType.target)
    
    markers = mean_ratio %>% 
        group_by(cellType.target) %>% 
        filter(rank_ratio <= n_per_group) %>% 
        ungroup()
    n_uniq = length(unique(markers$Symbol))
    print(paste0('Unique markers: ', n_uniq, ' of ', n_actual, ' (', round(100 * n_uniq / n_actual, 1), '%).'))
    
    markers = markers[match(unique(markers$Symbol), markers$Symbol),]
    print(paste0('Median ratio: ', median(markers$ratio)))
    
    return(unique(markers$Symbol))
}

tangram_markers <- get_markers(overlaps, 200)
tangram_markers

writeLines(
    tangram_markers,
    con = file.path(here::here(), "tangram_libd/processed-data/03_nn_run/markers.txt")
)

pdf(file.path(here::here(), "tangram_libd/plots/03_nn_run/marker_stats.pdf"))
#### Plot ####
ratio_plot <- map(
    marker_stats,
    ~ .x %>%
        ggplot(aes(ratio, std.logFC)) +
        geom_point(size = 0.5) +
        facet_wrap( ~ cellType.target, scales = "free_x") +
        labs(x = "mean(target logcount)/mean(highest non-target logcount)") +
        theme_bw() +
        NULL
)

ratio_plot

dev.off()
