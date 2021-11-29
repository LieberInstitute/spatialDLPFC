#!/bin/bash

## Load aws module
module load aws/2.2.35

## Upload main SpatialExperiment object
aws --profile bioturing-spatial s3 cp /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/spe_merged_final.Rdata s3://bioturing-spatial/spatialDLPFC/spe_merged_final.Rdata

## Locate raw images
grep "\-image" /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/code/spaceranger/*.sh
ls -lh /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/Images/round*/*.tif | grep -v cluster | grep -v nuclei

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
aws --profile bioturing-spatial s3 cp /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/Images/round4/V10B01-002_A1_Br2720_mid_DLPFC.tif s3://bioturing-spatial/spatialDLPFC/Images/round4/V10B01-002_A1_Br2720_mid_DLPFC.tif
aws --profile bioturing-spatial s3 cp /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/Images/round4/V10B01-002_B1_Br6432_ant_DLPFC.tif s3://bioturing-spatial/spatialDLPFC/Images/round4/V10B01-002_B1_Br6432_ant_DLPFC.tif
aws --profile bioturing-spatial s3 cp /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/Images/round4/V10B01-002_C1_Br8325_mid_DLPFC.tif s3://bioturing-spatial/spatialDLPFC/Images/round4/V10B01-002_C1_Br8325_mid_DLPFC.tif
aws --profile bioturing-spatial s3 cp /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/Images/round4/V10B01-002_D1_Br2720_ant_DLPFC.tif s3://bioturing-spatial/spatialDLPFC/Images/round4/V10B01-002_D1_Br2720_ant_DLPFC.tif


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
# 2021-11-29 11:04:09    3.1 GiB spatialDLPFC/Images/round4/V10B01-002_A1_Br2720_mid_DLPFC.tif
# 2021-11-29 11:04:57    3.1 GiB spatialDLPFC/Images/round4/V10B01-002_B1_Br6432_ant_DLPFC.tif
# 2021-11-29 11:05:46    3.1 GiB spatialDLPFC/Images/round4/V10B01-002_C1_Br8325_mid_DLPFC.tif
# 2021-11-29 11:06:33    3.1 GiB spatialDLPFC/Images/round4/V10B01-002_D1_Br2720_ant_DLPFC.tif
# 2021-11-29 10:21:12    1.5 GiB spatialDLPFC/spe_merged_final.Rdata
