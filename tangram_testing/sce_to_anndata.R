library('basilisk')
library('scRNAseq')
library('SingleCellExperiment')
library('zellkonverter')

#  Path to write the python AnnData object
# out_path = '/dcl01/lieber/ajaffe/Nick/spatial/tangram/sce_dlpfc.h5ad'
out_path = '/dcl01/lieber/ajaffe/Nick/spatial/tangram/visium_dlpfc.h5ad'

#  An example SingleCellExperiment object
# load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/forAmazonS3/SCE_DLPFC_tran-etal.rda")
load('/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/analysis/sce_combined.rda')

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

write_anndata(sce, out_path)

###############################################################################
#  Example from zellkonverter docs
###############################################################################

library(basilisk)
library(scRNAseq)
seger <- SegerstolpePancreasData()

# These functions are designed to be run inside
# a specified Python environment
roundtrip <- basiliskRun(fun = function(sce) {
    # Convert SCE to AnnData:
    adata <- SCE2AnnData(sce)
    
    # Maybe do some work in Python on'adata':
    # BLAH BLAH BLAH
    
    # Convert back to an SCE:
    AnnData2SCE(adata)
}, env = zellkonverter:::anndata_env, sce = seger)
