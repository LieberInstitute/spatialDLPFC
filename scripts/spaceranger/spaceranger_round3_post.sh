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
# fastq (NextSeq): /dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-07-09_ASpa061421_01/ and /dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-07-06_ASpa061421/
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
--id=DLPFC_Br6471_post_manual_alignment_all \
--transcriptome=/dcl02/lieber/ajaffe/SpatialTranscriptomics/refdata-gex-GRCh38-2020-A \
--fastqs=/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-07-09_ASpa061421_01/9_L002_ds.635a705b2c6e43abaa0e04b746c44c31,/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-07-09_ASpa061421_01/9_L001_ds.cbd0ef36129e468bafcfd20db0b79a20,/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-07-06_ASpa061421/Br6471_post_L002_ds.4d85a08656e24e89841d67e269935893,/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-07-06_ASpa061421/Br6471_post_L001_ds.e84086e206364d35990b5f46c48b9bb4 \
--image=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/round3/1000146_dlpfc_post_round3_A1.tif \
--slide=V10B01-053 \
--area=A1 \
--loupe-alignment=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/images_raw_align_json/round3/V10B01-053-A1_Br6471_post.json \
--jobmode=local \
--localcores=8 \
--localmem=64

spaceranger count \
--id=DLPFC_Br6522_post_manual_alignment_all \
--transcriptome=/dcl02/lieber/ajaffe/SpatialTranscriptomics/refdata-gex-GRCh38-2020-A \
--fastqs=/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-07-09_ASpa061421_01/10_L002_ds.c25895aff6a64d0aadec61b66376e7f5,/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-07-09_ASpa061421_01/10_L001_ds.ee8ee95852cb4f8cb47932641b76f45f,/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-07-06_ASpa061421/Br6522_post_L001_ds.83753540141649929e491149f6d21861,/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-07-06_ASpa061421/Br6522_post_L002_ds.82ad4d1a217f4e9580950d996238b314 \
--image=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/round3/1000146_dlpfc_post_round3_B1.tif \
--slide=V10B01-053 \
--area=B1 \
--loupe-alignment=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/images_raw_align_json/round3/V10B01-053-B1_Br6522_post.json \
--jobmode=local \
--localcores=8 \
--localmem=64

spaceranger count \
--id=DLPFC_Br8325_post_manual_alignment_all \
--transcriptome=/dcl02/lieber/ajaffe/SpatialTranscriptomics/refdata-gex-GRCh38-2020-A \
--fastqs=/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-07-09_ASpa061421_01/11_L002_ds.e53e6edd73244dbaac25291bb597338d,/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-07-09_ASpa061421_01/11_L001_ds.2afb581cacc143298347a45bb4ff8ec1,/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-07-06_ASpa061421/Br8325_post_L002_ds.fedaa7c0eef14127aedc3eb2afd82de4,/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-07-06_ASpa061421/Br8325_post_L001_ds.3007d815d09047afb4865f0d577fe7e5 \
--image=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/round3/1000146_dlpfc_post_round3_C1.tif \
--slide=V10B01-053 \
--area=C1 \
--loupe-alignment=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/images_raw_align_json/round3/V10B01-053-C1_Br8325_post.json \
--jobmode=local \
--localcores=8 \
--localmem=64

spaceranger count \
--id=DLPFC_Br8667_post_manual_alignment_all \
--transcriptome=/dcl02/lieber/ajaffe/SpatialTranscriptomics/refdata-gex-GRCh38-2020-A \
--fastqs=/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-07-09_ASpa061421_01/12_L002_ds.c7239ae756c44bd0b66474526d6a8c26,/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-07-09_ASpa061421_01/12_L001_ds.c300fd085f23466d8846aea932ec9823,/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-07-06_ASpa061421/Br8667_post_L002_ds.d81fb591b3f24525820be6f6c2d98857,/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-07-06_ASpa061421/Br8667_post_L001_ds.50875e92a52241f5b7c42a5a510f5ce0 \
--image=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/round3/1000146_dlpfc_post_round3_D1.tif \
--slide=V10B01-053 \
--area=D1 \
--loupe-alignment=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/images_raw_align_json/round3/V10B01-053-D1_Br8667_post.json \
--jobmode=local \
--localcores=8 \
--localmem=64


# restore working directory
cd $cwd