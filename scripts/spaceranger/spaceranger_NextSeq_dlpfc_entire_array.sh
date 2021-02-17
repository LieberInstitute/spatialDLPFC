#!/bin/bash
#$ -pe local 8
#$ -l mem_free=10G,h_vmem=20G,h_fsize=100G
#$ -V
#$ -cwd

# run on JHPCE cluster
# qsub scripts/spaceranger/spaceranger_NextSeq_dlpfc_tissue_alignment.sh

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

spaceranger count \
--id=DLPFC_Br2743_post_entire_array \
--transcriptome=/dcl02/lieber/ajaffe/SpatialTranscriptomics/refdata-gex-GRCh38-2020-A \
--fastqs=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/FASTQ/9_DLPFC_Br2743_post \
--image=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/Liebert_Institute_OTS-20-7748_rush_posterior_1.tif \
--slide=V19B23-075 \
--area=A1 \
--loupe-alignment=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/images_raw_align_json/V19B23-073-A1_9_Br2743_post.json \
--jobmode=local \
--localcores=8 \
--localmem=64

# spaceranger count \
# --id=DLPFC_Br2743_ant_entire_array \
# --transcriptome=/dcl02/lieber/ajaffe/SpatialTranscriptomics/refdata-gex-GRCh38-2020-A \
# --fastqs=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/FASTQ/1_DLPFC_Br2743_ant \
# --image=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/Liebert_Institute_OTS-20-7690_rush_anterior_1.tif \
# --slide=V19B23-075 \
# --area=A1 \
# --loupe-alignment=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/images_raw_align_json/V19B23-075-A1_1_Br2743_ant.json \
# --jobmode=local \
# --localcores=8 \
# --localmem=64

spaceranger count \
--id=DLPFC_3942_mid_entire_array \
--transcriptome=/dcl02/lieber/ajaffe/SpatialTranscriptomics/refdata-gex-GRCh38-2020-A \
--fastqs=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/FASTQ/6_DLPFC_Br3942_mid \
--image=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/Liebert_Institute_OTS-20-7748_rush-001_mid_2.tif \
--slide=V19B23-074 \
--area=B1 \
--loupe-alignment=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/images_raw_align_json/V19B23-074-B1_6_Br3942_mid.json \
--jobmode=local \
--localcores=8 \
--localmem=64

# restore working directory
cd $cwd