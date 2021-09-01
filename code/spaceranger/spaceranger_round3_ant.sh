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
# fastq (NextSeq): /dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-07-09_ASpa061421_01/  and /dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-07-06_ASpa061421/
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
--id=DLPFC_Br6471_ant_manual_alignment_all \
--transcriptome=/dcl02/lieber/ajaffe/SpatialTranscriptomics/refdata-gex-GRCh38-2020-A \
--fastqs=/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-07-09_ASpa061421_01/1_L001_ds.10ceb268148e48d39a388016db355a26,/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-07-09_ASpa061421_01/1_L002_ds.419c6990928344db96f69cd2ecf2a202,/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-07-06_ASpa061421/Br6471_ant_L001_ds.5fcf590a4b174c4693d4dad0eeaa5339,/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-07-06_ASpa061421/Br6471_ant_L002_ds.21e9938af86343c0912feb0191843e57 \
--image=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/round3/1000140_dlpfc_ant_round3_A1.tif \
--slide=V10U24-094 \
--area=A1 \
--loupe-alignment=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/images_raw_align_json/round3/V10U24-094-A1_Br6471_ant.json \
--jobmode=local \
--localcores=8 \
--localmem=64

spaceranger count \
--id=DLPFC_Br6522_ant_manual_alignment_all \
--transcriptome=/dcl02/lieber/ajaffe/SpatialTranscriptomics/refdata-gex-GRCh38-2020-A \
--fastqs=/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-07-09_ASpa061421_01/2_L002_ds.a90891a74f2d4f1794231f87f287d83a,/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-07-09_ASpa061421_01/2_L001_ds.2f7804989efd45488834df82a5aa5767,/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-07-06_ASpa061421/Br6522_ant_L002_ds.a2aa29ea8dfa497abb5f84a38c765ad4,/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-07-06_ASpa061421/Br6522_ant_L001_ds.ed89d24677ec48cb8aab1e40028add6c \
--image=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/round3/1000140_dlpfc_ant_round3_B1.tif \
--slide=V10U24-094 \
--area=B1 \
--loupe-alignment=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/images_raw_align_json/round3/V10U24-094-B1_Br6522_ant.json \
--jobmode=local \
--localcores=8 \
--localmem=64

spaceranger count \
--id=DLPFC_Br8325_ant_manual_alignment_all \
--transcriptome=/dcl02/lieber/ajaffe/SpatialTranscriptomics/refdata-gex-GRCh38-2020-A \
--fastqs=/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-07-09_ASpa061421_01/3_L002_ds.bbda1acf2d9e497db5d3a69be030c8b3,/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-07-09_ASpa061421_01/3_L001_ds.0a70da962c344a1da237b164df78a274,/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-07-06_ASpa061421/Br8325_ant_L002_ds.aae05687f92845dcb6ec71714bdc0b03,/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-07-06_ASpa061421/Br8325_ant_L001_ds.cece78473d63409ca31744302d3eb0a5 \
--image=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/round3/1000140_dlpfc_ant_round3_C1.tif \
--slide=V10U24-094 \
--area=C1 \
--loupe-alignment=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/images_raw_align_json/round3/V10U24-094-C1_Br8325_ant.json \
--jobmode=local \
--localcores=8 \
--localmem=64

spaceranger count \
--id=DLPFC_Br8667_ant_manual_alignment \
--transcriptome=/dcl02/lieber/ajaffe/SpatialTranscriptomics/refdata-gex-GRCh38-2020-A \
--fastqs=/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-07-09_ASpa061421_01/4_L002_ds.1f3c49a6bbae4232a591f2ca5dad96c9,/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-07-09_ASpa061421_01/4_L001_ds.92fff55d49904e7fa839a6096dbbba38,/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-07-06_ASpa061421/Br8667_ant_L002_ds.e24c104f32c143978d6f06d0bb69829a,/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-07-06_ASpa061421/Br8667_ant_L001_ds.3e939bf395104bc3973efe05adc26a9a \
--image=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/round3/1000140_dlpfc_ant_round3_D1.tif \
--slide=V10U24-094 \
--area=D1 \
--loupe-alignment=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/images_raw_align_json/round3/V10U24-094-D1_Br8667_ant.json \
--jobmode=local \
--localcores=8 \
--localmem=64


# restore working directory
cd $cwd