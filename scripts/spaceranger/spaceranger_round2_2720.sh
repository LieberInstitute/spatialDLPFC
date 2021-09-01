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
# fastq (NextSeq): /dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-05-21_ASpa050421/
# images (raw): /dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/

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
--id=DLPFC_Br2720_ant_manual_alignment \
--transcriptome=/dcl02/lieber/ajaffe/SpatialTranscriptomics/refdata-gex-GRCh38-2020-A \
--fastqs=/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-05-21_ASpa050421/1v_a_L001_ds.49395d5e9db64678ba034e8bfa98bc42,/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-05-21_ASpa050421/1v_a_L002_ds.ca72d0eb7d304f8fa2db3a726586bada \
--image=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/round2/V10U24-091_A1_Br2720_DLPFC_ant_1.tif \
--slide=V10U24-091 \
--area=A1 \
--loupe-alignment=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/images_raw_align_json/round2/V10U24-091-A1_Br2770_ant.json \
--jobmode=local \
--localcores=8 \
--localmem=64

spaceranger count \
--id=DLPFC_Br2720_mid_manual_alignment \
--transcriptome=/dcl02/lieber/ajaffe/SpatialTranscriptomics/refdata-gex-GRCh38-2020-A \
--fastqs=/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-05-21_ASpa050421/2v_a_L002_ds.ba75a8fa3cf547bcb38099fb77458a04,/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-05-21_ASpa050421/2v_a_L001_ds.eca9eb391a484f7590a99fc76d930d61 \
--image=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/round2/V10U24-091_B1_Br2720_DLPFC_mid_2.tif \
--slide=V10U24-091 \
--area=B1 \
--loupe-alignment=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/images_raw_align_json/round2/V10U24-091-B1_Br2770_mid.json \
--jobmode=local \
--localcores=8 \
--localmem=64

spaceranger count \
--id=DLPFC_Br2720_post_manual_alignment \
--transcriptome=/dcl02/lieber/ajaffe/SpatialTranscriptomics/refdata-gex-GRCh38-2020-A \
--fastqs=/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-05-21_ASpa050421/3v_a_L002_ds.fa84fd9e56f341508825c991393cdca3,/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-05-21_ASpa050421/3v_a_L001_ds.c7b5ef097ffa42ca99aa3c56ecc02654 \
--image=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/round2/V10U24-091_C1_Br2720_DLPFC_post_3.tif \
--slide=V10U24-091 \
--area=C1 \
--loupe-alignment=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/images_raw_align_json/round2/V10U24-091-C1_Br2770_post.json \
--jobmode=local \
--localcores=8 \
--localmem=64


# restore working directory
cd $cwd