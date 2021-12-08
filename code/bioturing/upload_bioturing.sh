#!/bin/bash

## Load aws module
module load aws/2.2.35

## Upload main SpatialExperiment object
aws --profile bioturing-spatial s3 cp /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/spe_merged_final.Rdata s3://bioturing-spatial/spatialDLPFC/spe_merged_final.Rdata

## Locate raw images
grep "\-image" /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/code/spaceranger/*.sh
ls -lh /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/Images/round*/*.tif | grep -v cluster | grep -v nuclei

## Locate images using the SpaceRanger ouptut (ignore 3 of them: the 2 duplicated and the 1 that had to be fixed)
for i in /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/NextSeq/*/_invocation /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/NextSeq/*/*/_invocation
do
    grep "tissue_image_paths" $i
done

## Upload images
aws --profile bioturing-spatial s3 cp /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/Images/round1/Liebert_Institute_OTS-20-7690_rush_anterior_1.tif s3://bioturing-spatial/spatialDLPFC/Images/round1/Liebert_Institute_OTS-20-7690_rush_anterior_1.tif
aws --profile bioturing-spatial s3 cp /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/Images/round1/Liebert_Institute_OTS-20-7690_rush_anterior_2.tif s3://bioturing-spatial/spatialDLPFC/Images/round1/Liebert_Institute_OTS-20-7690_rush_anterior_2.tif
aws --profile bioturing-spatial s3 cp /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/Images/round1/Liebert_Institute_OTS-20-7690_rush_anterior_3.tif s3://bioturing-spatial/spatialDLPFC/Images/round1/Liebert_Institute_OTS-20-7690_rush_anterior_3.tif
aws --profile bioturing-spatial s3 cp /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/Images/round1/Liebert_Institute_OTS-20-7690_rush_anterior_4.tif s3://bioturing-spatial/spatialDLPFC/Images/round1/Liebert_Institute_OTS-20-7690_rush_anterior_4.tif
aws --profile bioturing-spatial s3 cp /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/Images/round1/Liebert_Institute_OTS-20-7748_rush-001_mid_1.tif s3://bioturing-spatial/spatialDLPFC/Images/round1/Liebert_Institute_OTS-20-7748_rush-001_mid_1.tif
aws --profile bioturing-spatial s3 cp /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/Images/round1/Liebert_Institute_OTS-20-7748_rush-001_mid_2.tif s3://bioturing-spatial/spatialDLPFC/Images/round1/Liebert_Institute_OTS-20-7748_rush-001_mid_2.tif
aws --profile bioturing-spatial s3 cp /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/Images/round1/Liebert_Institute_OTS-20-7748_rush-001_mid_3.tif s3://bioturing-spatial/spatialDLPFC/Images/round1/Liebert_Institute_OTS-20-7748_rush-001_mid_3.tif
aws --profile bioturing-spatial s3 cp /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/Images/round1/Liebert_Institute_OTS-20-7748_rush-001_mid_4.tif s3://bioturing-spatial/spatialDLPFC/Images/round1/Liebert_Institute_OTS-20-7748_rush-001_mid_4.tif
aws --profile bioturing-spatial s3 cp /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/Images/round1/Liebert_Institute_OTS-20-7748_rush_posterior_1.tif s3://bioturing-spatial/spatialDLPFC/Images/round1/Liebert_Institute_OTS-20-7748_rush_posterior_1.tif
aws --profile bioturing-spatial s3 cp /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/Images/round1/Liebert_Institute_OTS-20-7748_rush_posterior_2.tif s3://bioturing-spatial/spatialDLPFC/Images/round1/Liebert_Institute_OTS-20-7748_rush_posterior_2.tif
aws --profile bioturing-spatial s3 cp /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/Images/round1/Liebert_Institute_OTS-20-7748_rush_posterior_3.tif s3://bioturing-spatial/spatialDLPFC/Images/round1/Liebert_Institute_OTS-20-7748_rush_posterior_3.tif
aws --profile bioturing-spatial s3 cp /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/Images/round1/Liebert_Institute_OTS-20-7748_rush_posterior_4.tif s3://bioturing-spatial/spatialDLPFC/Images/round1/Liebert_Institute_OTS-20-7748_rush_posterior_4.tif
aws --profile bioturing-spatial s3 cp /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/Images/round2/V10U24-091_A1_Br2720_DLPFC_ant_1.tif s3://bioturing-spatial/spatialDLPFC/Images/round2/V10U24-091_A1_Br2720_DLPFC_ant_1.tif
aws --profile bioturing-spatial s3 cp /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/Images/round2/V10U24-091_B1_Br2720_DLPFC_mid_2.tif s3://bioturing-spatial/spatialDLPFC/Images/round2/V10U24-091_B1_Br2720_DLPFC_mid_2.tif
aws --profile bioturing-spatial s3 cp /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/Images/round2/V10U24-091_C1_Br2720_DLPFC_post_3.tif s3://bioturing-spatial/spatialDLPFC/Images/round2/V10U24-091_C1_Br2720_DLPFC_post_3.tif
aws --profile bioturing-spatial s3 cp /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/Images/round2/V10U24-092_A1_Br6432_DLPFC_ant_1.tif s3://bioturing-spatial/spatialDLPFC/Images/round2/V10U24-092_A1_Br6432_DLPFC_ant_1.tif
aws --profile bioturing-spatial s3 cp /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/Images/round2/V10U24-092_B1_Br6432_DLPFC_mid_2.tif s3://bioturing-spatial/spatialDLPFC/Images/round2/V10U24-092_B1_Br6432_DLPFC_mid_2.tif
aws --profile bioturing-spatial s3 cp /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/Images/round2/V10U24-092_C1_Br6432_DLPFC_post_3.tif s3://bioturing-spatial/spatialDLPFC/Images/round2/V10U24-092_C1_Br6432_DLPFC_post_3.tif
aws --profile bioturing-spatial s3 cp /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/Images/round3/1000140_dlpfc_ant_round3_A1.tif s3://bioturing-spatial/spatialDLPFC/Images/round3/1000140_dlpfc_ant_round3_A1.tif
aws --profile bioturing-spatial s3 cp /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/Images/round3/1000140_dlpfc_ant_round3_B1.tif s3://bioturing-spatial/spatialDLPFC/Images/round3/1000140_dlpfc_ant_round3_B1.tif
aws --profile bioturing-spatial s3 cp /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/Images/round3/1000140_dlpfc_ant_round3_C1.tif s3://bioturing-spatial/spatialDLPFC/Images/round3/1000140_dlpfc_ant_round3_C1.tif
aws --profile bioturing-spatial s3 cp /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/Images/round3/1000140_dlpfc_ant_round3_D1.tif s3://bioturing-spatial/spatialDLPFC/Images/round3/1000140_dlpfc_ant_round3_D1.tif
aws --profile bioturing-spatial s3 cp /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/Images/round3/1000145_dlpfc_mid_round3_A1.tif s3://bioturing-spatial/spatialDLPFC/Images/round3/1000145_dlpfc_mid_round3_A1.tif
aws --profile bioturing-spatial s3 cp /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/Images/round3/1000145_dlpfc_mid_round3_B1.tif s3://bioturing-spatial/spatialDLPFC/Images/round3/1000145_dlpfc_mid_round3_B1.tif
aws --profile bioturing-spatial s3 cp /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/Images/round3/1000145_dlpfc_mid_round3_C1.tif s3://bioturing-spatial/spatialDLPFC/Images/round3/1000145_dlpfc_mid_round3_C1.tif
aws --profile bioturing-spatial s3 cp /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/Images/round3/1000145_dlpfc_mid_round3_D1.tif s3://bioturing-spatial/spatialDLPFC/Images/round3/1000145_dlpfc_mid_round3_D1.tif
aws --profile bioturing-spatial s3 cp /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/Images/round3/1000146_dlpfc_post_round3_A1.tif s3://bioturing-spatial/spatialDLPFC/Images/round3/1000146_dlpfc_post_round3_A1.tif
aws --profile bioturing-spatial s3 cp /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/Images/round3/1000146_dlpfc_post_round3_B1.tif s3://bioturing-spatial/spatialDLPFC/Images/round3/1000146_dlpfc_post_round3_B1.tif
aws --profile bioturing-spatial s3 cp /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/Images/round3/1000146_dlpfc_post_round3_C1.tif s3://bioturing-spatial/spatialDLPFC/Images/round3/1000146_dlpfc_post_round3_C1.tif
aws --profile bioturing-spatial s3 cp /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/Images/round3/1000146_dlpfc_post_round3_D1.tif s3://bioturing-spatial/spatialDLPFC/Images/round3/1000146_dlpfc_post_round3_D1.tif
## Don't upload this image since it's not used in the SPE object
#aws --profile bioturing-spatial s3 cp /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/Images/round4/V10B01-002_A1_Br2720_mid_DLPFC.tif s3://bioturing-spatial/spatialDLPFC/Images/round4/V10B01-002_A1_Br2720_mid_DLPFC.tif
# aws --profile bioturing-spatial s3 rm s3://bioturing-spatial/spatialDLPFC/Images/round4/V10B01-002_A1_Br2720_mid_DLPFC.tif
aws --profile bioturing-spatial s3 cp /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/Images/round4/V10B01-002_B1_Br6432_ant_DLPFC.tif s3://bioturing-spatial/spatialDLPFC/Images/round4/V10B01-002_B1_Br6432_ant_DLPFC.tif
aws --profile bioturing-spatial s3 cp /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/Images/round4/V10B01-002_C1_Br8325_mid_DLPFC.tif s3://bioturing-spatial/spatialDLPFC/Images/round4/V10B01-002_C1_Br8325_mid_DLPFC.tif
aws --profile bioturing-spatial s3 cp /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/Images/round4/V10B01-002_D1_Br2720_ant_DLPFC.tif s3://bioturing-spatial/spatialDLPFC/Images/round4/V10B01-002_D1_Br2720_ant_DLPFC.tif

## Locate SpaceRanger _invocation files and upload them
## These files will enable matching with the raw images as shown below:
# $ grep "tissue_image_paths" /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/NextSeq/Round3/DLPFC_Br8667_ant_manual_alignment_all/_invocation
#     tissue_image_paths        = ["/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/round3/1000140_dlpfc_ant_round3_D1.tif"],
for i in /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/NextSeq/*/_invocation /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/NextSeq/*/*/_invocation
    do echo $i
    dir1=$(dirname $i)
    dirname=$(basename $dir1)
    filename=$(basename $i)
    #echo "$i s3://bioturing-spatial/spatialDLPFC/SpaceRanger/${dirname}/${filename}"
    aws --profile bioturing-spatial s3 cp $i s3://bioturing-spatial/spatialDLPFC/SpaceRanger/${dirname}/${filename}
done

## Upload scalefactors_json.json files
for i in /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/NextSeq/*/outs/spatial/scalefactors_json.json /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/NextSeq/*/*/outs/spatial/scalefactors_json.json
    do echo $i
    dir1=$(dirname $i)
    dir2=$(dirname $dir1)
    dir3=$(dirname $dir2)
    dirname=$(basename $dir3)
    filename=$(basename $i)
    #echo "$i s3://bioturing-spatial/spatialDLPFC/SpaceRanger/${dirname}/outs/spatial/${filename}"
    aws --profile bioturing-spatial s3 cp $i s3://bioturing-spatial/spatialDLPFC/SpaceRanger/${dirname}/outs/spatial/${filename}
done

## Ignored SpaceRanger output dirs:
# /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/NextSeq/Round4/DLPFC_Br2720_ant_2_old/
# /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/NextSeq/Round3/DLPFC_Br8667_ant_manual_alignment_all/
# /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/NextSeq/Round2/DLPFC_Br2720_post_manual_alignment/

## Remove directory we are not using. This one had an incorrect input image:
## it wasn't aligned well and was breaking VistoSeg as detailed at
## https://github.com/LieberInstitute/spatialDLPFC/issues/39
aws --profile bioturing-spatial s3 rm --recursive s3://bioturing-spatial/spatialDLPFC/SpaceRanger/DLPFC_Br2720_ant_2_old

## Remove duplicated ones: we sequenced more reads for these ones
aws --profile bioturing-spatial s3 rm --recursive s3://bioturing-spatial/spatialDLPFC/SpaceRanger/DLPFC_Br8667_ant_manual_alignment_all
aws --profile bioturing-spatial s3 rm --recursive s3://bioturing-spatial/spatialDLPFC/SpaceRanger/DLPFC_Br2720_post_manual_alignment


## Locate loupe alignment json files
for i in /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/images_raw_align_json/*.json /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/images_raw_align_json/*/*.json
    do echo $i
    filename=$(basename $i)
    #echo "$i s3://bioturing-spatial/spatialDLPFC/images_raw_align_json/${filename}"
    aws --profile bioturing-spatial s3 cp $i s3://bioturing-spatial/spatialDLPFC/images_raw_align_json/${filename}
done

## Re-upload on 2021-12-07 new tif image that Heena made on 2021-12-01 which
## was 2 days after the initial image upload
aws --profile bioturing-spatial s3 cp /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/Images/round4/V10B01-002_D1_Br2720_ant_DLPFC.tif s3://bioturing-spatial/spatialDLPFC/Images/round4/V10B01-002_D1_Br2720_ant_DLPFC.tif
## Companion new SpaceRanger output
# /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/NextSeq/Round4/DLPFC_Br2720_ant_2

## TODO (waiting for the new spe from Abby)
## Upload main SpatialExperiment object
# aws --profile bioturing-spatial s3 cp /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/spe_merged_final.Rdata s3://bioturing-spatial/spatialDLPFC/spe_merged_final.Rdata

## List files on the S3 bucket
aws --profile bioturing-spatial s3 ls --recursive --human-readable s3://bioturing-spatial/spatialDLPFC/
# 2021-11-29 10:37:52    2.5 GiB spatialDLPFC/Images/round1/Liebert_Institute_OTS-20-7690_rush_anterior_1.tif
# 2021-11-29 10:40:15    2.3 GiB spatialDLPFC/Images/round1/Liebert_Institute_OTS-20-7690_rush_anterior_2.tif
# 2021-11-29 10:41:00    2.4 GiB spatialDLPFC/Images/round1/Liebert_Institute_OTS-20-7690_rush_anterior_3.tif
# 2021-11-29 10:41:47    2.4 GiB spatialDLPFC/Images/round1/Liebert_Institute_OTS-20-7690_rush_anterior_4.tif
# 2021-11-29 10:42:31    2.2 GiB spatialDLPFC/Images/round1/Liebert_Institute_OTS-20-7748_rush-001_mid_1.tif
# 2021-11-29 10:43:14    2.2 GiB spatialDLPFC/Images/round1/Liebert_Institute_OTS-20-7748_rush-001_mid_2.tif
# 2021-11-29 10:43:57    2.2 GiB spatialDLPFC/Images/round1/Liebert_Institute_OTS-20-7748_rush-001_mid_3.tif
# 2021-11-29 10:44:43    2.2 GiB spatialDLPFC/Images/round1/Liebert_Institute_OTS-20-7748_rush-001_mid_4.tif
# 2021-11-29 10:45:26    3.0 GiB spatialDLPFC/Images/round1/Liebert_Institute_OTS-20-7748_rush_posterior_1.tif
# 2021-11-29 10:46:17    3.0 GiB spatialDLPFC/Images/round1/Liebert_Institute_OTS-20-7748_rush_posterior_2.tif
# 2021-11-29 10:50:01    2.6 GiB spatialDLPFC/Images/round1/Liebert_Institute_OTS-20-7748_rush_posterior_3.tif
# 2021-11-29 10:50:44    3.6 GiB spatialDLPFC/Images/round1/Liebert_Institute_OTS-20-7748_rush_posterior_4.tif
# 2021-11-29 10:51:42    2.1 GiB spatialDLPFC/Images/round2/V10U24-091_A1_Br2720_DLPFC_ant_1.tif
# 2021-11-29 10:52:21    2.0 GiB spatialDLPFC/Images/round2/V10U24-091_B1_Br2720_DLPFC_mid_2.tif
# 2021-11-29 10:52:57    2.1 GiB spatialDLPFC/Images/round2/V10U24-091_C1_Br2720_DLPFC_post_3.tif
# 2021-11-29 10:53:33    1.8 GiB spatialDLPFC/Images/round2/V10U24-092_A1_Br6432_DLPFC_ant_1.tif
# 2021-11-29 10:54:09    1.8 GiB spatialDLPFC/Images/round2/V10U24-092_B1_Br6432_DLPFC_mid_2.tif
# 2021-11-29 10:54:42    1.7 GiB spatialDLPFC/Images/round2/V10U24-092_C1_Br6432_DLPFC_post_3.tif
# 2021-11-29 10:55:39    2.2 GiB spatialDLPFC/Images/round3/1000140_dlpfc_ant_round3_A1.tif
# 2021-11-29 10:56:19    2.1 GiB spatialDLPFC/Images/round3/1000140_dlpfc_ant_round3_B1.tif
# 2021-11-29 10:56:57    2.0 GiB spatialDLPFC/Images/round3/1000140_dlpfc_ant_round3_C1.tif
# 2021-11-29 10:57:33    2.0 GiB spatialDLPFC/Images/round3/1000140_dlpfc_ant_round3_D1.tif
# 2021-11-29 10:58:07    2.9 GiB spatialDLPFC/Images/round3/1000145_dlpfc_mid_round3_A1.tif
# 2021-11-29 10:58:57    2.6 GiB spatialDLPFC/Images/round3/1000145_dlpfc_mid_round3_B1.tif
# 2021-11-29 10:59:40    2.6 GiB spatialDLPFC/Images/round3/1000145_dlpfc_mid_round3_C1.tif
# 2021-11-29 11:00:24    2.6 GiB spatialDLPFC/Images/round3/1000145_dlpfc_mid_round3_D1.tif
# 2021-11-29 11:01:23    2.5 GiB spatialDLPFC/Images/round3/1000146_dlpfc_post_round3_A1.tif
# 2021-11-29 11:02:04    2.4 GiB spatialDLPFC/Images/round3/1000146_dlpfc_post_round3_B1.tif
# 2021-11-29 11:02:46    2.4 GiB spatialDLPFC/Images/round3/1000146_dlpfc_post_round3_C1.tif
# 2021-11-29 11:03:28    2.4 GiB spatialDLPFC/Images/round3/1000146_dlpfc_post_round3_D1.tif
# 2021-11-29 11:04:57    3.1 GiB spatialDLPFC/Images/round4/V10B01-002_B1_Br6432_ant_DLPFC.tif
# 2021-11-29 11:05:46    3.1 GiB spatialDLPFC/Images/round4/V10B01-002_C1_Br8325_mid_DLPFC.tif
# 2021-12-07 10:18:41    3.5 GiB spatialDLPFC/Images/round4/V10B01-002_D1_Br2720_ant_DLPFC.tif
# 2021-12-07 12:22:10    2.2 KiB spatialDLPFC/SpaceRanger/DLPFC_Br2720_ant_2/_invocation
# 2021-12-07 12:30:07  166 Bytes spatialDLPFC/SpaceRanger/DLPFC_Br2720_ant_2/outs/spatial/scalefactors_json.json
# 2021-12-07 12:19:36    2.2 KiB spatialDLPFC/SpaceRanger/DLPFC_Br2720_ant_manual_alignment/_invocation
# 2021-12-07 12:27:28  166 Bytes spatialDLPFC/SpaceRanger/DLPFC_Br2720_ant_manual_alignment/outs/spatial/scalefactors_json.json
# 2021-12-07 12:19:44    2.2 KiB spatialDLPFC/SpaceRanger/DLPFC_Br2720_mid_manual_alignment/_invocation
# 2021-12-07 12:27:36  148 Bytes spatialDLPFC/SpaceRanger/DLPFC_Br2720_mid_manual_alignment/outs/spatial/scalefactors_json.json
# 2021-12-07 12:19:53    3.3 KiB spatialDLPFC/SpaceRanger/DLPFC_Br2720_post_extra_reads/_invocation
# 2021-12-07 12:27:44  157 Bytes spatialDLPFC/SpaceRanger/DLPFC_Br2720_post_extra_reads/outs/spatial/scalefactors_json.json
# 2021-12-07 12:17:47    1.5 KiB spatialDLPFC/SpaceRanger/DLPFC_Br2743_ant_manual_alignment/_invocation
# 2021-12-07 12:25:53  165 Bytes spatialDLPFC/SpaceRanger/DLPFC_Br2743_ant_manual_alignment/outs/spatial/scalefactors_json.json
# 2021-12-07 12:18:00    2.7 KiB spatialDLPFC/SpaceRanger/DLPFC_Br2743_mid_manual_alignment_extra_reads/_invocation
# 2021-12-07 12:26:01  146 Bytes spatialDLPFC/SpaceRanger/DLPFC_Br2743_mid_manual_alignment_extra_reads/outs/spatial/scalefactors_json.json
# 2021-12-07 12:18:08    1.5 KiB spatialDLPFC/SpaceRanger/DLPFC_Br2743_post_manual_alignment/_invocation
# 2021-12-07 12:26:09  167 Bytes spatialDLPFC/SpaceRanger/DLPFC_Br2743_post_manual_alignment/outs/spatial/scalefactors_json.json
# 2021-12-07 12:18:21    1.5 KiB spatialDLPFC/SpaceRanger/DLPFC_Br3942_ant_manual_alignment/_invocation
# 2021-12-07 12:26:17  157 Bytes spatialDLPFC/SpaceRanger/DLPFC_Br3942_ant_manual_alignment/outs/spatial/scalefactors_json.json
# 2021-12-07 12:18:30    1.5 KiB spatialDLPFC/SpaceRanger/DLPFC_Br3942_mid_manual_alignment/_invocation
# 2021-12-07 12:26:25  155 Bytes spatialDLPFC/SpaceRanger/DLPFC_Br3942_mid_manual_alignment/outs/spatial/scalefactors_json.json
# 2021-12-07 12:18:38    1.5 KiB spatialDLPFC/SpaceRanger/DLPFC_Br3942_post_manual_alignment/_invocation
# 2021-12-07 12:26:33  157 Bytes spatialDLPFC/SpaceRanger/DLPFC_Br3942_post_manual_alignment/outs/spatial/scalefactors_json.json
# 2021-12-07 12:18:46    2.7 KiB spatialDLPFC/SpaceRanger/DLPFC_Br6423_ant_manual_alignment_extra_reads/_invocation
# 2021-12-07 12:26:41  157 Bytes spatialDLPFC/SpaceRanger/DLPFC_Br6423_ant_manual_alignment_extra_reads/outs/spatial/scalefactors_json.json
# 2021-12-07 12:18:54    1.5 KiB spatialDLPFC/SpaceRanger/DLPFC_Br6423_mid_manual_alignment/_invocation
# 2021-12-07 12:26:48  163 Bytes spatialDLPFC/SpaceRanger/DLPFC_Br6423_mid_manual_alignment/outs/spatial/scalefactors_json.json
# 2021-12-07 12:19:02    2.7 KiB spatialDLPFC/SpaceRanger/DLPFC_Br6423_post_extra_reads/_invocation
# 2021-12-07 12:26:57  167 Bytes spatialDLPFC/SpaceRanger/DLPFC_Br6423_post_extra_reads/outs/spatial/scalefactors_json.json
# 2021-12-07 12:22:26    2.2 KiB spatialDLPFC/SpaceRanger/DLPFC_Br6432_ant_2/_invocation
# 2021-12-07 12:30:15  157 Bytes spatialDLPFC/SpaceRanger/DLPFC_Br6432_ant_2/outs/spatial/scalefactors_json.json
# 2021-12-07 12:20:01    2.2 KiB spatialDLPFC/SpaceRanger/DLPFC_Br6432_ant_manual_alignment/_invocation
# 2021-12-07 12:27:52  157 Bytes spatialDLPFC/SpaceRanger/DLPFC_Br6432_ant_manual_alignment/outs/spatial/scalefactors_json.json
# 2021-12-07 12:20:09    2.2 KiB spatialDLPFC/SpaceRanger/DLPFC_Br6432_mid_manual_alignment/_invocation
# 2021-12-07 12:28:00  166 Bytes spatialDLPFC/SpaceRanger/DLPFC_Br6432_mid_manual_alignment/outs/spatial/scalefactors_json.json
# 2021-12-07 12:20:17    2.2 KiB spatialDLPFC/SpaceRanger/DLPFC_Br6432_post_manual_alignment/_invocation
# 2021-12-07 12:28:08  157 Bytes spatialDLPFC/SpaceRanger/DLPFC_Br6432_post_manual_alignment/outs/spatial/scalefactors_json.json
# 2021-12-07 12:20:25    3.3 KiB spatialDLPFC/SpaceRanger/DLPFC_Br6471_ant_manual_alignment_all/_invocation
# 2021-12-07 12:28:16  164 Bytes spatialDLPFC/SpaceRanger/DLPFC_Br6471_ant_manual_alignment_all/outs/spatial/scalefactors_json.json
# 2021-12-07 12:20:33    3.3 KiB spatialDLPFC/SpaceRanger/DLPFC_Br6471_mid_manual_alignment_all/_invocation
# 2021-12-07 12:28:24  166 Bytes spatialDLPFC/SpaceRanger/DLPFC_Br6471_mid_manual_alignment_all/outs/spatial/scalefactors_json.json
# 2021-12-07 12:20:41    3.3 KiB spatialDLPFC/SpaceRanger/DLPFC_Br6471_post_manual_alignment_all/_invocation
# 2021-12-07 12:28:32  166 Bytes spatialDLPFC/SpaceRanger/DLPFC_Br6471_post_manual_alignment_all/outs/spatial/scalefactors_json.json
# 2021-12-07 12:20:49    3.3 KiB spatialDLPFC/SpaceRanger/DLPFC_Br6522_ant_manual_alignment_all/_invocation
# 2021-12-07 12:28:40  156 Bytes spatialDLPFC/SpaceRanger/DLPFC_Br6522_ant_manual_alignment_all/outs/spatial/scalefactors_json.json
# 2021-12-07 12:20:58    3.3 KiB spatialDLPFC/SpaceRanger/DLPFC_Br6522_mid_manual_alignment_all/_invocation
# 2021-12-07 12:28:48  149 Bytes spatialDLPFC/SpaceRanger/DLPFC_Br6522_mid_manual_alignment_all/outs/spatial/scalefactors_json.json
# 2021-12-07 12:21:06    3.3 KiB spatialDLPFC/SpaceRanger/DLPFC_Br6522_post_manual_alignment_all/_invocation
# 2021-12-07 12:28:56  156 Bytes spatialDLPFC/SpaceRanger/DLPFC_Br6522_post_manual_alignment_all/outs/spatial/scalefactors_json.json
# 2021-12-07 12:21:14    3.3 KiB spatialDLPFC/SpaceRanger/DLPFC_Br8325_ant_manual_alignment_all/_invocation
# 2021-12-07 12:29:04  156 Bytes spatialDLPFC/SpaceRanger/DLPFC_Br8325_ant_manual_alignment_all/outs/spatial/scalefactors_json.json
# 2021-12-07 12:22:34    2.2 KiB spatialDLPFC/SpaceRanger/DLPFC_Br8325_mid_2/_invocation
# 2021-12-07 12:30:23  157 Bytes spatialDLPFC/SpaceRanger/DLPFC_Br8325_mid_2/outs/spatial/scalefactors_json.json
# 2021-12-07 12:21:22    3.3 KiB spatialDLPFC/SpaceRanger/DLPFC_Br8325_mid_manual_alignment_all/_invocation
# 2021-12-07 12:29:12  149 Bytes spatialDLPFC/SpaceRanger/DLPFC_Br8325_mid_manual_alignment_all/outs/spatial/scalefactors_json.json
# 2021-12-07 12:21:30    3.3 KiB spatialDLPFC/SpaceRanger/DLPFC_Br8325_post_manual_alignment_all/_invocation
# 2021-12-07 12:29:20  148 Bytes spatialDLPFC/SpaceRanger/DLPFC_Br8325_post_manual_alignment_all/outs/spatial/scalefactors_json.json
# 2021-12-07 12:19:10    1.5 KiB spatialDLPFC/SpaceRanger/DLPFC_Br8492_ant_manual_alignment/_invocation
# 2021-12-07 12:27:04  157 Bytes spatialDLPFC/SpaceRanger/DLPFC_Br8492_ant_manual_alignment/outs/spatial/scalefactors_json.json
# 2021-12-07 12:19:19    2.7 KiB spatialDLPFC/SpaceRanger/DLPFC_Br8492_mid_manual_alignment_extra_reads/_invocation
# 2021-12-07 12:27:12  146 Bytes spatialDLPFC/SpaceRanger/DLPFC_Br8492_mid_manual_alignment_extra_reads/outs/spatial/scalefactors_json.json
# 2021-12-07 12:19:27    1.5 KiB spatialDLPFC/SpaceRanger/DLPFC_Br8492_post_manual_alignment/_invocation
# 2021-12-07 12:27:20  158 Bytes spatialDLPFC/SpaceRanger/DLPFC_Br8492_post_manual_alignment/outs/spatial/scalefactors_json.json
# 2021-12-07 12:21:38    3.3 KiB spatialDLPFC/SpaceRanger/DLPFC_Br8667_ant_extra_reads/_invocation
# 2021-12-07 12:29:27  156 Bytes spatialDLPFC/SpaceRanger/DLPFC_Br8667_ant_extra_reads/outs/spatial/scalefactors_json.json
# 2021-12-07 12:21:54    3.3 KiB spatialDLPFC/SpaceRanger/DLPFC_Br8667_mid_manual_alignment_all/_invocation
# 2021-12-07 12:29:43  167 Bytes spatialDLPFC/SpaceRanger/DLPFC_Br8667_mid_manual_alignment_all/outs/spatial/scalefactors_json.json
# 2021-12-07 12:22:02    3.3 KiB spatialDLPFC/SpaceRanger/DLPFC_Br8667_post_manual_alignment_all/_invocation
# 2021-12-07 12:29:51  157 Bytes spatialDLPFC/SpaceRanger/DLPFC_Br8667_post_manual_alignment_all/outs/spatial/scalefactors_json.json
# 2021-12-07 19:24:56   12.6 KiB spatialDLPFC/code/check_dates.txt
# 2021-12-07 19:22:15   27.9 KiB spatialDLPFC/code/upload_bioturing.sh
# 2021-12-07 12:51:32  556.4 KiB spatialDLPFC/images_raw_align_json/V10B01-002-B1.json
# 2021-12-07 12:51:40  558.8 KiB spatialDLPFC/images_raw_align_json/V10B01-002-C1.json
# 2021-12-07 12:51:48  545.7 KiB spatialDLPFC/images_raw_align_json/V10B01-002-D1.json
# 2021-12-07 12:49:51  564.7 KiB spatialDLPFC/images_raw_align_json/V10B01-052-A1_Br6471_mid.json
# 2021-12-07 12:49:59  554.0 KiB spatialDLPFC/images_raw_align_json/V10B01-052-B1_Br6522_mid.json
# 2021-12-07 12:50:08  553.4 KiB spatialDLPFC/images_raw_align_json/V10B01-052-C1_Br8325_mid.json
# 2021-12-07 12:50:16  560.1 KiB spatialDLPFC/images_raw_align_json/V10B01-052-D1_Br8667_mid.json
# 2021-12-07 12:50:24  557.2 KiB spatialDLPFC/images_raw_align_json/V10B01-053-A1_Br6471_post.json
# 2021-12-07 12:50:33  555.1 KiB spatialDLPFC/images_raw_align_json/V10B01-053-B1_Br6522_post.json
# 2021-12-07 12:50:41  562.3 KiB spatialDLPFC/images_raw_align_json/V10B01-053-C1_Br8325_post.json
# 2021-12-07 12:50:50  562.2 KiB spatialDLPFC/images_raw_align_json/V10B01-053-D1_Br8667_post.json
# 2021-12-07 12:49:00  546.9 KiB spatialDLPFC/images_raw_align_json/V10U24-091-A1_Br2770_ant.json
# 2021-12-07 12:49:08  527.3 KiB spatialDLPFC/images_raw_align_json/V10U24-091-B1_Br2770_mid.json
# 2021-12-07 12:49:17  566.4 KiB spatialDLPFC/images_raw_align_json/V10U24-091-C1_Br2770_post.json
# 2021-12-07 12:49:25  548.3 KiB spatialDLPFC/images_raw_align_json/V10U24-092-A1_Br6432_ant.json
# 2021-12-07 12:49:34  550.1 KiB spatialDLPFC/images_raw_align_json/V10U24-092-B1_Br6432_mid.json
# 2021-12-07 12:49:42  536.6 KiB spatialDLPFC/images_raw_align_json/V10U24-092-C1_Br6432_post.json
# 2021-12-07 12:50:58  540.8 KiB spatialDLPFC/images_raw_align_json/V10U24-094-A1_Br6471_ant.json
# 2021-12-07 12:51:06  560.8 KiB spatialDLPFC/images_raw_align_json/V10U24-094-B1_Br6522_ant.json
# 2021-12-07 12:51:15  550.3 KiB spatialDLPFC/images_raw_align_json/V10U24-094-C1_Br8325_ant.json
# 2021-12-07 12:51:23  552.0 KiB spatialDLPFC/images_raw_align_json/V10U24-094-D1_Br8667_ant.json
# 2021-12-07 12:47:18  555.6 KiB spatialDLPFC/images_raw_align_json/V19B23-073-A1_9_Br2743_post_manual_alignment.json
# 2021-12-07 12:47:26  562.2 KiB spatialDLPFC/images_raw_align_json/V19B23-073-B1_10_Br3942_post_manual_alignment.json
# 2021-12-07 12:47:35  550.4 KiB spatialDLPFC/images_raw_align_json/V19B23-073-C1_11_Br6423_post_manual_alignment.json
# 2021-12-07 12:47:43  560.8 KiB spatialDLPFC/images_raw_align_json/V19B23-073-D1_12_Br8492_post_manual_alignment.json
# 2021-12-07 12:47:52  560.3 KiB spatialDLPFC/images_raw_align_json/V19B23-074-A1_5_Br2743_mid_manual_alignment.json
# 2021-12-07 12:48:00  555.9 KiB spatialDLPFC/images_raw_align_json/V19B23-074-B1_6_Br3942_mid_manual_alignment.json
# 2021-12-07 12:48:09  556.1 KiB spatialDLPFC/images_raw_align_json/V19B23-074-C1_7_Br6432_mid_manual_alignment.json
# 2021-12-07 12:48:17  562.8 KiB spatialDLPFC/images_raw_align_json/V19B23-074-D1_8_Br8492_mid_manual_alignment.json
# 2021-12-07 12:48:25  554.1 KiB spatialDLPFC/images_raw_align_json/V19B23-075-A1_1_Br2743_ant_manual_alignment.json
# 2021-12-07 12:48:35  555.6 KiB spatialDLPFC/images_raw_align_json/V19B23-075-B1_2_Br3942_ant_manual_alignment.json
# 2021-12-07 12:48:43  555.4 KiB spatialDLPFC/images_raw_align_json/V19B23-075-C1_3_Br6423_ant_manual_alignment.json
# 2021-12-07 12:48:51  567.5 KiB spatialDLPFC/images_raw_align_json/V19B23-075-D1_4_Br8492_ant_manual_alignment.json
# 2021-11-29 10:21:12    1.5 GiB spatialDLPFC/spe_merged_final.Rdata


## Upload this script
aws --profile bioturing-spatial s3 cp /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/code/bioturing/check_dates.txt s3://bioturing-spatial/spatialDLPFC/code/check_dates.txt
aws --profile bioturing-spatial s3 cp /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/code/bioturing/upload_bioturing.sh s3://bioturing-spatial/spatialDLPFC/code/upload_bioturing.sh

