#!/bin/bash
#$ -pe local 8
#$ -l bluejay,mem_free=10G,h_vmem=20G,h_fsize=100G
#$ -V
#$ -cwd

# run on JHPCE cluster
# qsub scripts/spaceranger/spaceranger_round3_ant.sh

############################
# Script to run Space Ranger
############################

# locations of files:
# -------------------
# spaceranger reference: /dcl02/lieber/ajaffe/SpatialTranscriptomics/refdata-gex-GRCh38-2020-A
# fastq (NextSeq): /dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-07-09_ASpa061421_01/ and /dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-07-06_ASpa061421/
# images (raw): /dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/round3

# summary spreadsheet: /dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC
# - contains sample ID, sample name, slide serial number, capture area ID


# load spaceranger module
module use /jhpce/shared/jhpce/modulefiles/libd
module load spaceranger

# run in outputs directory (spaceranger can only save outputs in current working directory)
cwd=$(pwd)
# mkdir -p outputs/NextSeq
cd outputs/NextSeq/Round2


# run spaceranger count for each sample

spaceranger count \
--id=DLPFC_Br6471_mid_manual_alignment_all \
--transcriptome=/dcl02/lieber/ajaffe/SpatialTranscriptomics/refdata-gex-GRCh38-2020-A \
--fastqs=/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-07-09_ASpa061421_01/5_L002_ds.35e41e3c45c746f1943a08c3103e66cc,/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-07-09_ASpa061421_01/5_L001_ds.6d5df2c002f5449382759148ffca5ac4,/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-07-06_ASpa061421/Br6471_mid_L001_ds.7627a93f0f674d7dbd5ffc3b2c25842d,/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-07-06_ASpa061421/Br6471_mid_L002_ds.4ba78f3a7c044ca6a215e81c6ed493fe \
--image=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/round3/1000145_dlpfc_mid_round3_A1.tif \
--slide=V10B01-052 \
--area=A1 \
--loupe-alignment=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/images_raw_align_json/round3/V10B01-052-A1_Br6471_mid.json \
--jobmode=local \
--localcores=8 \
--localmem=64

spaceranger count \
--id=DLPFC_Br6522_mid_manual_alignment_all \
--transcriptome=/dcl02/lieber/ajaffe/SpatialTranscriptomics/refdata-gex-GRCh38-2020-A \
--fastqs=/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-07-09_ASpa061421_01/6_L002_ds.53c25680150c49ff8a3a40f0cf041d42,/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-07-09_ASpa061421_01/6_L001_ds.71618150320f455cad41e6a245f6ac3c,/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-07-06_ASpa061421/Br6522_mid_L002_ds.18cbee8fa2224f36b49a3c65ba6bb66d,/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-07-06_ASpa061421/Br6522_mid_L001_ds.7ba583fb13f94353be0a359a67df5887 \
--image=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/round3/1000145_dlpfc_mid_round3_B1.tif \
--slide=V10B01-052 \
--area=B1 \
--loupe-alignment=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/images_raw_align_json/round3/V10B01-052-B1_Br6522_mid.json \
--jobmode=local \
--localcores=8 \
--localmem=64

spaceranger count \
--id=DLPFC_Br8325_mid_manual_alignment_all \
--transcriptome=/dcl02/lieber/ajaffe/SpatialTranscriptomics/refdata-gex-GRCh38-2020-A \
--fastqs=/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-07-09_ASpa061421_01/7_L001_ds.6af5d4c7249342c58ef838a2bf6f785a,/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-07-09_ASpa061421_01/7_L002_ds.331ebbe8e28e41efba0e093568a70dda,/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-07-06_ASpa061421/Br8325_mid_L002_ds.dc9cc7652c7040469a73d37bbc2d8b5f,/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-07-06_ASpa061421/Br8325_mid_L001_ds.9258cd4b7691494fb5ce6ca781d0e386 \
--image=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/round3/1000145_dlpfc_mid_round3_C1.tif \
--slide=V10B01-052 \
--area=C1 \
--loupe-alignment=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/images_raw_align_json/round3/V10B01-052-C1_Br8325_mid.json \
--jobmode=local \
--localcores=8 \
--localmem=64

spaceranger count \
--id=DLPFC_Br8667_mid_manual_alignment_all \
--transcriptome=/dcl02/lieber/ajaffe/SpatialTranscriptomics/refdata-gex-GRCh38-2020-A \
--fastqs=/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-07-09_ASpa061421_01/8_L002_ds.ec1844ae28cb47d29095adccb5918ef6,/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-07-09_ASpa061421_01/8_L001_ds.c3c52ea52b654094a0f1551194150b71,/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-07-06_ASpa061421/Br8667_mid_L002_ds.a99814d7e0184b0984d00c593cf5e251,/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-07-06_ASpa061421/Br8667_mid_L001_ds.3427254002e14d3d8e2d957e11254ef3 \
--image=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/round3/1000145_dlpfc_mid_round3_D1.tif \
--slide=V10B01-052 \
--area=D1 \
--loupe-alignment=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/images_raw_align_json/round3/V10B01-052-D1_Br8667_mid.json \
--jobmode=local \
--localcores=8 \
--localmem=64


# restore working directory
cd $cwd