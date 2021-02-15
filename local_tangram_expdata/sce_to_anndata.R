library('basilisk')
library('scRNAseq')
library('SingleCellExperiment')
library('zellkonverter')
library('data.table')
library('entropy')
library('tidyverse')

setwd("/home/arta/Documents/GitHub/spython/local_tangram_expdata")
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

rowData(rna.sce)$rowSums <- rowSums(as.data.table(rowData(rna.sce))[, 3:20])

# entropy(as.numeric(as.data.frame(rowData(rna.sce)[,3:19])))

markers_entropy <- as.data.table(rowData(rna.sce)) %>% 
    mutate(entropy = as.numeric(pmap(as.data.frame(rowData(rna.sce)[,3:19]), lift_vd(entropy)))) %>% 
    drop_na() %>% 
    arrange(desc(entropy)) %>% 
    filter(entropy != 0) %>% 
    select(Symbol, entropy) %>% 
    slice_head(., prop = .10)

markers <- intersect(rowData(spatial.seq)$gene_name, as.character(markers_entropy$Symbol)) 
    
# rowData(rna.sce[rowSums(as.data.table(rowData(rna.sce))[, 3:20]) != 0,])

# > hist(subset(rowData(rna.sce)$rowSums, rowData(rna.sce)$rowSums < 3))
# > hist(subset(rowData(rna.sce)$rowSums, rowData(rna.sce)$rowSums < 1))
# > hist(subset(rowData(rna.sce)$rowSums, rowData(rna.sce)$rowSums < .5))
# > hist(subset(rowData(rna.sce)$rowSums, rowData(rna.sce)$rowSums < .1))

rm(sce, sce.dlpfc)

# data.table(n = 1:length(markers), markers)

write.csv(as.data.table(markers), file = "out/markers.csv")

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

write_anndata(spatial.seq, visium_out)
write_anndata(rna.sce, sc_out)
