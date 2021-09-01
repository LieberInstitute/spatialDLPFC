#!/bin/bash
#$ -pe local 8
#$ -l mem_free=10G,h_vmem=20G,h_fsize=100G
#$ -o /dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/Liebert_Institute_OTS-20-7690_rush_anterior_3.txt
#$ -e /dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/Liebert_Institute_OTS-20-7690_rush_anterior_3.txt
#$ -V
#$ -cwd


# Specify memory and other details. Note that 'mem_free' and 'h_vmem' specify
# per-core memory (12G * 5 cores = 60GB total, as we want), as indicated here:
# https://jhpce.jhu.edu/knowledge-base/how-to/#multicore

#  Make LIBD modules available, and load the "spaceranger" module
module use /jhpce/shared/jhpce/modulefiles/libd
module load spaceranger

#  The main Space Ranger command
spaceranger count \
--id=Liebert_Institute_OTS-20-7690_rush_anterior_3 \
--transcriptome=/dcl02/lieber/ajaffe/SpatialTranscriptomics/refdata-gex-GRCh38-2020-A \
--fastqs=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/FASTQ/3_DLPFC_Br6423_ant \
--image=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/Liebert_Institute_OTS-20-7690_rush_anterior_3.tif \
--slide=V19B23-075 \
--area=C1 \
--loupe-alignment=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/Lieber_Institute_OTS-20-7690_rush_anterior_3_V19B23-075-C1.json \
--jobmode=local \
--localcores=8 \
--localmem=64
