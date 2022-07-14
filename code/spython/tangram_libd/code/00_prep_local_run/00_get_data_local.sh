#!/bin/sh

#  Download and unzip data required for the newer tutorial:
#  https://github.com/broadinstitute/Tangram/blob/master/example/1_tutorial_tangram.ipynb

mkdir -p ../../raw-data/00_get_data_local/

echo "Downloading scRNAseq and spatial data..."
sftp aseyedia@jhpce-transfer01.jhsph.edu << EOF
get /dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/analysis/sce_combined.rda ../../raw-data/00_get_data_local/sce_combined.rda
get /dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/forAmazonS3/SCE_DLPFC_tran-etal.rda ../../raw-data/00_get_data_local/SCE_DLPFC_tran-etal.rda
bye
EOF

echo "../../raw-data/00_get_data_local/*" >> ../../raw-data/00_get_data_local/.gitignore

exit 0
