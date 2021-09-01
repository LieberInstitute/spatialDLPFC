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
--id=DLPFC_Br6432_ant_manual_alignment \
--transcriptome=/dcl02/lieber/ajaffe/SpatialTranscriptomics/refdata-gex-GRCh38-2020-A \
--fastqs=/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-05-21_ASpa050421/5v_a_L002_ds.9bbe35bb540e4087a0c7d59229f329f5,/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-05-21_ASpa050421/5v_a_L001_ds.b41990d1c68a498d8365087034b8d676 \
--image=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/round2/V10U24-092_A1_Br6432_DLPFC_ant_1.tif \
--slide=V10U24-092 \
--area=A1 \
--loupe-alignment=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/images_raw_align_json/round2/V10U24-092-A1_Br6432_ant.json \
--jobmode=local \
--localcores=8 \
--localmem=64

# spaceranger count \
# --id=DLPFC_Br6423_mid_manual_alignment \
# --transcriptome=/dcl02/lieber/ajaffe/SpatialTranscriptomics/refdata-gex-GRCh38-2020-A \
# --fastqs=/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-05-21_ASpa050421/6v_a_L002_ds.298bdefd8a1d4ab59c37c44bd34649b7,/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-05-21_ASpa050421/6v_a_L001_ds.d8a36b89304f426c9adbe8fcbb5ce06b \
# --image=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/round2/V10U24-092_B1_Br6432_DLPFC_mid_2.tif \
# --slide=V10U24-092 \
# --area=B1 \
# --loupe-alignment=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/images_raw_align_json/round2/V10U24-092-B1_Br6432_mid.json \
# --jobmode=local \
# --localcores=8 \
# --localmem=64

# spaceranger count \
# --id=DLPFC_Br6432_post_manual_alignment \
# --transcriptome=/dcl02/lieber/ajaffe/SpatialTranscriptomics/refdata-gex-GRCh38-2020-A \
# --fastqs=/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-05-21_ASpa050421/7v_a_L001_ds.8fd218130938421faecc3eba3bcbae83,/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-05-21_ASpa050421/7v_a_L002_ds.c2d6d4b32c234d9aa4ce24f4b609c242 \
# --image=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/round2/V10U24-092_C1_Br6432_DLPFC_post_3.tif \
# --slide=V10U24-092 \
# --area=C1 \
# --loupe-alignment=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/images_raw_align_json/round2/V10U24-092-C1_Br6432_post.json \
# --jobmode=local \
# --localcores=8 \
# --localmem=64


# restore working directory
cd $cwd