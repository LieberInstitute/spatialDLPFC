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
cd outputs/NextSeq


# run spaceranger count for each sample

# spaceranger count \
# --id=DLPFC_Br2720_post_extra_reads \
# --transcriptome=/dcl02/lieber/ajaffe/SpatialTranscriptomics/refdata-gex-GRCh38-2020-A \
# --fastqs=/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-05-21_ASpa050421/3v_a_L002_ds.fa84fd9e56f341508825c991393cdca3,/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-05-21_ASpa050421/3v_a_L001_ds.c7b5ef097ffa42ca99aa3c56ecc02654,/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-09-22_ASpa050421/r2_sample3v_a_L002_ds.05ba070833ee454bab4244cd7b0fdd35,/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-09-22_ASpa050421/r2_sample3v_a_L001_ds.8ed358ff823446d0815eb3af112cc4df \
# --image=/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/Images/round2/V10U24-091_C1_Br2720_DLPFC_post_3.tif \
# --slide=V10U24-091 \
# --area=C1 \
# --loupe-alignment=/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/images_raw_align_json/round2/V10U24-091-C1_Br2770_post.json \
# --jobmode=local \
# --localcores=8 \
# --localmem=64

spaceranger count \
--id=DLPFC_Br8667_ant_extra_reads \
--transcriptome=/dcl02/lieber/ajaffe/SpatialTranscriptomics/refdata-gex-GRCh38-2020-A \
--fastqs=/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-07-09_ASpa061421_01/4_L002_ds.1f3c49a6bbae4232a591f2ca5dad96c9,/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-07-09_ASpa061421_01/4_L001_ds.92fff55d49904e7fa839a6096dbbba38,/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-09-22_ASpa061421_01/r3_sample4_L001_ds.6350f2352ce242c8acfb9629a828f438,/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-09-22_ASpa061421_01/r3_sample4_L002_ds.578af6f8560b46f6814a6f24cdcddf35 \
--image=/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/Images/round3/1000140_dlpfc_ant_round3_D1.tif \
--slide=V10U24-094 \
--area=D1 \
--loupe-alignment=/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/images_raw_align_json/round3/V10U24-094-D1_Br8667_ant.json \
--jobmode=local \
--localcores=8 \
--localmem=64


# restore working directory
cd $cwd