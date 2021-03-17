library('basilisk')
library('scRNAseq')
library('SingleCellExperiment')
library('zellkonverter')
library('data.table')
library('entropy')
library('tidyverse')
library("DeconvoBuddies")
library("here")

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

# overlaps <- rowData(spatial.seq)[rowData(spatial.seq)$gene_name %in% rowData(rna.sce)$Symbol,]
overlaps <- rna.sce[rowData(rna.sce)$Symbol %in% rowData(spatial.seq)$gene_name,]

# # Getting row sums
rowData(overlaps)$rowSums <- rowSums(as.data.table(rowData(overlaps))[, 3:20])
# 
# # taking only highly-expressed genes
rowData(overlaps) <- rowData(overlaps)[order(rowData(overlaps)$rowSums, decreasing = TRUE),]
# 
# # setting threshold to 5
overlaps <- overlaps[rowData(overlaps)$rowSums > 5,]

mean_ratio <- map("cell_type", ~get_mean_ratio2(overlaps, cellType_col = .x, assay_name = "counts", add_symbol = TRUE))
markers_1vAll <- map("cell_type", ~findMarkers_1vAll(overlaps, cellType_col = .x, assay_name = "counts", add_symbol = TRUE))

marker_stats <- map2(mean_ratio, markers_1vAll, 
                     ~left_join(.x, .y, by = c("gene", "cellType.target", "Symbol"))) 

# > dim(marker_stats %>% as.data.table())
# [1] 20239    16

tangram_markers <- marker_stats %>% as.data.table() %>% group_by("cellType") %>% slice_head(prop = 0.05) %>% ungroup() %>% select("Symbol")
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
