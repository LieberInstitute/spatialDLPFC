#!/bin/sh

#  Download and unzip data required for the newer tutorial:
#  https://github.com/broadinstitute/Tangram/blob/master/example/1_tutorial_tangram.ipynb

mkdir -p ../../raw-data/01_prepare_tangram/

cp /dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/analysis/sce_combined.rda ../../raw-data/01_prepare_tangram/sce_combined.rda
cp /dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/forAmazonS3/SCE_DLPFC_tran-etal.rda ../../raw-data/01_prepare_tangram/SCE_DLPFC_tran-etal.rda

exit 0
