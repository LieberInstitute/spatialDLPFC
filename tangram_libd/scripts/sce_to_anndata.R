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
dir.create(file.path(here::here(), "tangram_libd/out/"), showWarnings = FALSE)
visium_out = file.path(here::here(), "tangram_libd/out/visium_dlpfc.h5ad")
sc_out = file.path(here::here(), "tangram_libd/out/sce_dlpfc.h5ad")

#  An example SingleCellExperiment object
load(file.path(here::here(), "tangram_libd/data/SCE_DLPFC_tran-etal.rda"))
load(file.path(here::here(), "tangram_libd/data/sce_combined.rda"))

rna.sce <- sce.dlpfc
spatial.seq <- sce

# prepairing for expression_cutoff
seed <- 20210324

# Find expression cutoff for each dataset and filter out genes below that
# cutoff
counts_rna = as.matrix(assays(rna.sce)$counts)
cutoff_rna <- max(expression_cutoff(counts_rna, seed = seed))
rna.sce = rna.sce[rowMeans(counts_rna) > cutoff_rna,]

counts_spat = as.matrix(assays(spatial.seq)$counts)
cutoff_spat <- max(expression_cutoff(counts_spat, seed = seed))
spatial.seq = spatial.seq[rowMeans(counts_spat) > cutoff_spat,]

# Getting overlaps between spatial and scRNAseq data
overlaps <- rna.sce[rowData(rna.sce)$Symbol %in% rowData(spatial.seq)$gene_name, ]

mean_ratio <- map("cell_type", ~get_mean_ratio2(overlaps, cellType_col = .x, assay_name = "counts", add_symbol = TRUE))
markers_1vAll <- map("cell_type", ~findMarkers_1vAll(overlaps, cellType_col = .x, assay_name = "counts", add_symbol = TRUE))

marker_stats <- map2(mean_ratio, markers_1vAll, 
                     ~left_join(.x, .y, by = c("gene", "cellType.target", "Symbol"))) 

tangram_markers <- marker_stats %>% as.data.table() %>% group_by("cellType") %>% slice_head(prop = 0.01) %>% ungroup() %>% select("Symbol")
tangram_markers

write.csv(tangram_markers, 
          file = file.path(here::here(), "tangram_libd/data/marker_stats.csv"),
          quote = FALSE,
          row.names = FALSE,
          col.names = FALSE)

pdf(file.path(here::here(), "tangram_libd/out/marker_stats.pdf"))
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
