#!/bin/bash
#$ -pe local 8
#$ -l mem_free=10G,h_vmem=20G,h_fsize=100G
#$ -V
#$ -cwd

# run on JHPCE cluster
# qsub scripts/spaceranger/spaceranger_NextSeq_round1_resequence.sh

############################
# Script to run Space Ranger
############################

# locations of files:
# -------------------
# spaceranger reference: /dcl02/lieber/ajaffe/SpatialTranscriptomics/refdata-gex-GRCh38-2020-A
# fastq (NextSeq): /dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-06-17_ASpa041621/
# images (raw): /dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/

# summary spreadsheet: /dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/sample_info/Visium dlpfc 100520 Master.xlsx
# - contains sample ID, sample name, slide serial number, capture area ID


# load spaceranger module
module use /jhpce/shared/jhpce/modulefiles/libd
module load spaceranger

# run in outputs directory (spaceranger can only save outputs in current working directory)
cwd=$(pwd)
# mkdir -p outputs/NextSeq
cd processed-data/NextSeq


# run spaceranger count for each sample


spaceranger count \
--id=DLPFC_Br8325_mid_2 \
--transcriptome=/dcl02/lieber/ajaffe/SpatialTranscriptomics/refdata-gex-GRCh38-2020-A \
--fastqs=/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-09-22_ASpa082421/3_L001_ds.ba452df9e0b743009daca0811834d5de,/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-09-22_ASpa082421/3_L002_ds.35b89d8bb7214b40805f4799a9b4c943 \
--image=/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/Images/round4/V10B01-002_C1_Br8325_mid_DLPFC.tif \
--slide=V10B01-002 \
--area=C1 \
--loupe-alignment=/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/images_raw_align_json/round4/V10B01-002-C1.json \
--jobmode=local \
--localcores=8 \
--localmem=64

spaceranger count \
--id=DLPFC_Br2720_ant_2 \
--transcriptome=/dcl02/lieber/ajaffe/SpatialTranscriptomics/refdata-gex-GRCh38-2020-A \
--fastqs=/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-09-22_ASpa082421/4_L001_ds.7d02274e402d4635bc56094ffa6e0728,/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-09-22_ASpa082421/4_L002_ds.52b7909924c143458260ea881e9bad52 \
--image=/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/Images/round4/V10B01-002_D1_Br2720_ant_DLPFC.tif \
--slide=V10B01-002 \
--area=D1 \
--loupe-alignment=/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/images_raw_align_json/round4/V10B01-002-D1.json \
--jobmode=local \
--localcores=8 \
--localmem=64

spaceranger count \
--id=DLPFC_Br6432_ant_2 \
--transcriptome=/dcl02/lieber/ajaffe/SpatialTranscriptomics/refdata-gex-GRCh38-2020-A \
--fastqs=/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-09-22_ASpa082421/2_L001_ds.9bb0a982bc824035a270a0da8a43d250,/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-09-22_ASpa082421/2_L002_ds.58eddff4077348c58f82102f63d23000 \
--image=/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/Images/round4/V10B01-002_B1_Br6432_ant_DLPFC.tif \
--slide=V10B01-002 \
--area=B1 \
--loupe-alignment=/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/images_raw_align_json/round4/V10B01-002-B1.json \
--jobmode=local \
--localcores=8 \
--localmem=64


# restore working directory
cd $cwd