library('basilisk')
library('scRNAseq')
library('SingleCellExperiment')
library('zellkonverter')
library('data.table')
library('entropy')
library('tidyverse')
library("DeconvoBuddies")
library("here")
library("jaffelab")

args = commandArgs(trailingOnly=TRUE)

# args[1] = "JHPCE"

if (length(args) == 0) {
    print("Running on local PC")
    setwd("/home/arta/Documents/GitHub/spython/local_tangram_expdata")
} else if (args[1] == "JHPCE"){
    print("Running on JHPCE")
    setwd("/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spython/local_tangram_expdata")
}
#  Path to write the python AnnData object
# out_path = '/dcl01/lieber/ajaffe/Nick/spatial/tangram/sce_dlpfc.h5ad'
dir.create("out/", showWarnings = FALSE)
visium_out = 'out/visium_dlpfc.h5ad'
sc_out = 'out/sce_dlpfc.h5ad'

#  An example SingleCellExperiment object
load("data/SCE_DLPFC_tran-etal.rda")
load('data/sce_combined.rda')

rna.sce <- sce.dlpfc
spatial.seq <- sce

# Getting overlaps between spatial and scRNAseq data
overlaps <- rna.sce[rowData(rna.sce)$Symbol %in% rowData(spatial.seq)$gene_name, ]

# prepairing for expression_cutoff
seed <- 20210324
olaps_rpkm <- as.matrix(rowData(overlaps)[, 3:20])

# generating cutoff value
cutoff <- max(expression_cutoff(olaps_rpkm, seed = seed))

# generating rowMeans
## based on https://github.com/LieberInstitute/goesHyde_mdd_rnaseq/blob/3ee0ba2a77f2d3a111dba3e81c28594cdd4aa46f/exprs_cutoff/get_expression_cutoffs.R
rowData(overlaps)$rowMeans <- rowMeans(olaps_rpkm)

# subset to only rows with means greater than cutoff
overlaps <- overlaps[rowData(overlaps)$rowMeans > cutoff,]

# sort table in descending rowMeans order (for fun)
rowData(overlaps) <- rowData(overlaps)[order(rowData(overlaps)$rowMeans, decreasing = TRUE), ]

mean_ratio <- map("cell_type", ~get_mean_ratio2(overlaps, cellType_col = .x, assay_name = "counts", add_symbol = TRUE))
markers_1vAll <- map("cell_type", ~findMarkers_1vAll(overlaps, cellType_col = .x, assay_name = "counts", add_symbol = TRUE))

marker_stats <- map2(mean_ratio, markers_1vAll, 
                     ~left_join(.x, .y, by = c("gene", "cellType.target", "Symbol"))) 

tangram_markers <- marker_stats %>% as.data.table() %>% group_by("cellType") %>% slice_head(prop = 0.01) %>% ungroup() %>% select("Symbol")
tangram_markers

write.csv(tangram_markers, file =  "data/marker_stats.csv")

# save(marker_stats, file =  "data/marker_stats.Rdata")
# load(here("data", "marker_stats.Rdata"), verbose = TRUE)

dir.create("analysis/")

pdf("analysis/marker_stats.pdf")
#### Plot ####
ratio_plot <- map(marker_stats,
                  ~.x %>%
                      ggplot(aes(ratio, std.logFC)) +
                      geom_point(size = 0.5) +
                      facet_wrap(~cellType.target, scales = "free_x") +
                      labs(x = "mean(target logcount)/mean(highest non-target logcount)") +
                      theme_bw()+
                      NULL)

ratio_plot

dev.off()

rm(sce, sce.dlpfc)

###############################################################################
#  The main code we'll use in general to convert SCE R objects to AnnData
#  python objects, as a preprocessing step to running tangram on our own data
###############################################################################

write_anndata = function(sce, out_path) {
    invisible(basiliskRun(fun = function(sce, filename) {
        library('zellkonverter')
        library('reticulate')
        
        # Convert SCE to AnnData:
        adata <- SCE2AnnData(sce)
        
        #  Write AnnData object to disk
        adata$write(filename=filename)
        
        return()
    }, env = zellkonverter:::anndata_env, sce = sce, filename = out_path))
}

write_anndata(rna.sce, sc_out)
write_anndata(spatial.seq, visium_out)
