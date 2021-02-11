#!/bin/sh

#  Download and unzip data required for the newer tutorial:
#  https://github.com/broadinstitute/Tangram/blob/master/example/1_tutorial_tangram.ipynb
mkdir data/

echo "Downloading scRNAseq data..."
sftp aseyedia@jhpce-transfer01.jhsph.edu
get /dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/analysis/sce_combined.rda data/sce_combined.rda

echo "Downloading spatial data..."
get /dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/forAmazonS3/SCE_DLPFC_tran-etal.rda

exit 0
