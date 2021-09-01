#!/bin/bash
#$ -pe local 8
#$ -l mem_free=10G,h_vmem=20G,h_fsize=100G
#$ -V
#$ -cwd

# run on JHPCE cluster
# qsub scripts/spaceranger/spaceranger_NextSeq_dlpfc_ant.sh

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
# mkdir -p outputs/NextSeq
cd outputs/NextSeq


# run spaceranger count for each sample

# spaceranger count \
# --id=DLPFC_Br2743_ant_manual_alignment \
# --transcriptome=/dcl02/lieber/ajaffe/SpatialTranscriptomics/refdata-gex-GRCh38-2020-A \
# --fastqs=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/FASTQ/2_DLPFC_Br3942_ant \
# --image=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/Liebert_Institute_OTS-20-7690_rush_anterior_1.tif \
# --slide=V19B23-075 \
# --area=A1 \
# --loupe-alignment=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/images_raw_align_json/V19B23-075-A1_1_Br2743_ant_manual_alignment.json \
# --jobmode=local \
# --localcores=8 \
# --localmem=64

# spaceranger count \
# --id=DLPFC_3942_ant_manual_alignment \
# --transcriptome=/dcl02/lieber/ajaffe/SpatialTranscriptomics/refdata-gex-GRCh38-2020-A \
# --fastqs=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/FASTQ/3_DLPFC_Br6423_ant \
# --image=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/Liebert_Institute_OTS-20-7690_rush_anterior_2.tif \
# --slide=V19B23-075 \
# --area=B1 \
# --loupe-alignment=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/images_raw_align_json/V19B23-075-B1_2_Br3942_ant_manual_alignment.json \
# --jobmode=local \
# --localcores=8 \
# --localmem=64

spaceranger count \
--id=DLPFC_Br6423_ant_manual_alignment_extra_reads \
--transcriptome=/dcl02/lieber/ajaffe/SpatialTranscriptomics/refdata-gex-GRCh38-2020-A \
--fastqs=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/FASTQ/1_DLPFC_Br2743_ant,/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-06-17_ASpa041621/Br2743_ant_L002_ds.fc72ee6430d840508be719765f7405d5,/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-06-17_ASpa041621/Br2743_ant_L001_ds.8f7de7ab26264e01b7511eca6ce97a81 \
--image=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/round1/Liebert_Institute_OTS-20-7690_rush_anterior_3.tif \
--slide=V19B23-075 \
--area=C1 \
--loupe-alignment=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/images_raw_align_json/V19B23-075-C1_3_Br6423_ant_manual_alignment.json \
--jobmode=local \
--localcores=8 \
--localmem=64

# spaceranger count \
# --id=DLPFC_Br8492_ant_manual_alignment \
# --transcriptome=/dcl02/lieber/ajaffe/SpatialTranscriptomics/refdata-gex-GRCh38-2020-A \
# --fastqs=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/FASTQ/4_DLPFC_Br8492_ant \
# --image=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/Liebert_Institute_OTS-20-7690_rush_anterior_4.tif \
# --slide=V19B23-075 \
# --area=D1 \
# --loupe-alignment=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/images_raw_align_json/V19B23-075-D1_4_Br8492_ant_manual_alignment.json \
# --jobmode=local \
# --localcores=8 \
# --localmem=64

# restore working directory
cd $cwd