#!/bin/bash
#$ -pe local 8
#$ -l mem_free=10G,h_vmem=20G,h_fsize=100G
#$ -V
#$ -cwd

# run on JHPCE cluster
# qsub scripts/spaceranger/spaceranger_NextSeq_dlpfc_mid.sh

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
--id=DLPFC_Br2743_mid_manual_alignment_extra_reads \
--transcriptome=/dcl02/lieber/ajaffe/SpatialTranscriptomics/refdata-gex-GRCh38-2020-A \
--fastqs=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/FASTQ/8_DLPFC_Br8492_mid,/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-06-17_ASpa041621/Br8492_mid_L002_ds.deb2d4a1944a4a83a79f976bfb2c581a,/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-06-17_ASpa041621/Br8492_mid_L001_ds.b97357d3c73f4d89b975fcec46022c4d \
--image=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/round1/Liebert_Institute_OTS-20-7748_rush-001_mid_1.tif \
--slide=V19B23-074 \
--area=A1 \
--loupe-alignment=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/images_raw_align_json/V19B23-074-A1_5_Br2743_mid_manual_alignment.json \
--jobmode=local \
--localcores=8 \
--localmem=64

# spaceranger count \
# --id=DLPFC_3942_mid_manual_alignment \
# --transcriptome=/dcl02/lieber/ajaffe/SpatialTranscriptomics/refdata-gex-GRCh38-2020-A \
# --fastqs=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/FASTQ/9_DLPFC_Br2743_post \
# --image=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/Liebert_Institute_OTS-20-7748_rush-001_mid_2.tif \
# --slide=V19B23-074 \
# --area=B1 \
# --loupe-alignment=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/images_raw_align_json/V19B23-074-B1_6_Br3942_mid_manual_alignment.json \
# --jobmode=local \
# --localcores=8 \
# --localmem=64

# spaceranger count \
# --id=DLPFC_Br6423_mid_manual_alignment \
# --transcriptome=/dcl02/lieber/ajaffe/SpatialTranscriptomics/refdata-gex-GRCh38-2020-A \
# --fastqs=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/FASTQ/10_DLPFC_Br3942_post \
# --image=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/Liebert_Institute_OTS-20-7748_rush-001_mid_3.tif \
# --slide=V19B23-074 \
# --area=C1 \
# --loupe-alignment=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/images_raw_align_json/V19B23-074-C1_7_Br6432_mid_manual_alignment.json \
# --jobmode=local \
# --localcores=8 \
# --localmem=64

# spaceranger count \
# --id=DLPFC_Br8492_mid_manual_alignment_extra_reads \
# --transcriptome=/dcl02/lieber/ajaffe/SpatialTranscriptomics/refdata-gex-GRCh38-2020-A \
# --fastqs=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/FASTQ/11_DLPFC_Br6423_post,/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-06-17_ASpa041621/Br6423_post_L002_ds.9bb2573626b143a9803a060f1689b33b,/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-06-17_ASpa041621/Br6423_post_L001_ds.8d0297fe7adf4d1a9ee283ac92db27db \
# --image=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/round1/Liebert_Institute_OTS-20-7748_rush-001_mid_4.tif \
# --slide=V19B23-074 \
# --area=D1 \
# --loupe-alignment=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/images_raw_align_json/V19B23-074-D1_8_Br8492_mid_manual_alignment.json \
# --jobmode=local \
# --localcores=8 \
# --localmem=64

# restore working directory
cd $cwd