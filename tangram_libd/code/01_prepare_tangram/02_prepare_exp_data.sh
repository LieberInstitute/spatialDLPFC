#!/bin/sh

#  Download and unzip data required for the newer tutorial:
#  https://github.com/broadinstitute/Tangram/blob/master/example/1_tutorial_tangram.ipynb

mkdir -p ../../raw-data/01_prepare_tangram/

#  Spatial data from the SpatialDLPFC repo and early snRNA-seq data from Matt
cp /dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/analysis/sce_combined.rda ../../raw-data/01_prepare_tangram/sce_combined.rda
cp /dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/forAmazonS3/SCE_DLPFC_tran-etal.rda ../../raw-data/01_prepare_tangram/SCE_DLPFC_tran-etal.rda

#  12 spatial samples (pilot DLPFC samples published in NN) and new 3-donor snRNAseq from Matt
cp /dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/Analysis/spatialLIBD_files/Human_DLPFC_Visium_processedData_sce_scran_spatialLIBD.Rdata ../../raw-data/01_prepare_tangram/
cp /dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/revision/regionSpecific_DLPFC-n3_cleaned-combined_SCE_LAH2021.rda ../../raw-data/01_prepare_tangram/

exit 0
