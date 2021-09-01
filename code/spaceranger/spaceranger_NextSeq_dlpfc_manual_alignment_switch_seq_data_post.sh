#!/bin/bash
#$ -pe local 8
#$ -l mem_free=10G,h_vmem=20G,h_fsize=100G
#$ -V
#$ -cwd

# run on JHPCE cluster
# qsub scripts/spaceranger/spaceranger_NextSeq_dlpfc_post.sh

############################
# Script to run Space Ranger
############################

# locations of files:
# -------------------
# spaceranger reference: /dcl02/lieber/ajaffe/SpatialTranscriptomics/refdata-gex-GRCh38-2020-A
# fastq (NextSeq): /dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/FASTQ/
# images (raw): /dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/

# summary spreadsheet: /dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/sample_info/Visium dlpfc 100520 Master.xlsx
# - contains sample ID, sample name, slide serial number, capture area ID


# load spaceranger module
module use /jhpce/shared/jhpce/modulefiles/libd
module load spaceranger

# run in outputs directory (spaceranger can only save outputs in current working directory)
cwd=$(pwd)
#mkdir -p outputs/NextSeq
cd outputs/NextSeq


# run spaceranger count for each sample

# spaceranger count \
# --id=DLPFC_Br2743_post_manual_alignment \
# --transcriptome=/dcl02/lieber/ajaffe/SpatialTranscriptomics/refdata-gex-GRCh38-2020-A \
# --fastqs=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/FASTQ/5_DLPFC_Br2743_mid \
# --image=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/Liebert_Institute_OTS-20-7748_rush_posterior_1.tif \
# --slide=V19B23-073 \
# --area=A1 \
# --loupe-alignment=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/images_raw_align_json/V19B23-073-A1_9_Br2743_post_manual_alignment.json \
# --jobmode=local \
# --localcores=8 \
# --localmem=64

# spaceranger count \
# --id=DLPFC_3942_post_manual_alignment \
# --transcriptome=/dcl02/lieber/ajaffe/SpatialTranscriptomics/refdata-gex-GRCh38-2020-A \
# --fastqs=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/FASTQ/6_DLPFC_Br3942_mid \
# --image=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/Liebert_Institute_OTS-20-7748_rush_posterior_2.tif \
# --slide=V19B23-073 \
# --area=B1 \
# --loupe-alignment=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/images_raw_align_json/V19B23-073-B1_10_Br3942_post_manual_alignment.json \
# --jobmode=local \
# --localcores=8 \
# --localmem=64

# spaceranger count \
# --id=DLPFC_Br6423_post_manual_alignment \
# --transcriptome=/dcl02/lieber/ajaffe/SpatialTranscriptomics/refdata-gex-GRCh38-2020-A \
# --fastqs=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/FASTQ/7_DLPFC_Br6423_mid \
# --image=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/Liebert_Institute_OTS-20-7748_rush_posterior_3.tif \
# --slide=V19B23-073 \
# --area=C1 \
# --loupe-alignment=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/images_raw_align_json/V19B23-073-C1_11_Br6423_post_manual_alignment.json \
# --jobmode=local \
# --localcores=8 \
# --localmem=64

spaceranger count \
--id=DLPFC_Br8492_post_manual_alignment \
--transcriptome=/dcl02/lieber/ajaffe/SpatialTranscriptomics/refdata-gex-GRCh38-2020-A \
--fastqs=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/FASTQ/12_DLPFC_Br8492_post \
--image=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/Liebert_Institute_OTS-20-7748_rush_posterior_4.tif \
--slide=V19B23-073 \
--area=D1 \
--loupe-alignment=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/images_raw_align_json/V19B23-073-D1_12_Br8492_post_manual_alignment.json \
--jobmode=local \
--localcores=8 \
--localmem=64

# restore working directory
cd $cwd