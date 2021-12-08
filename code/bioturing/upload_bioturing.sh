#!/bin/bash

## Load aws module
module load aws/2.2.35


## TODO (waiting for the new spe from Abby)
## Upload main SpatialExperiment object
aws --profile bioturing-spatial s3 cp /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/spe_merged_final.Rdata s3://bioturing-spatial/spatialDLPFC/spe_merged_final.Rdata

## Locate images using the SpaceRanger output
for i in /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/NextSeq/*/_invocation /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/NextSeq/*/*/_invocation
do
    grep "tissue_image_paths" $i
done

## Ignored SpaceRanger output dirs:
# $ ls -lh /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/NextSeq_not_used
# total 98K
# drwxrws--- 4 aspangle lieber_lcolladotor 19 Oct 21 02:51 DLPFC_Br2720_ant_2_old
# drwxrws--- 4 lcollado lieber_lcolladotor 19 Jul  8 13:48 DLPFC_Br2720_ant_manual_alignment
# drwxrws--- 4 lcollado lieber_lcolladotor 19 Jul 15 19:57 DLPFC_Br6432_ant_manual_alignment
# drwxrws--- 4 lcollado lieber_lcolladotor 19 Jul 27 23:57 DLPFC_Br8325_mid_manual_alignment_all
# drwxrws--- 4 lcollado lieber_lcolladotor 19 Jul 28 00:57 DLPFC_Br8667_ant_manual_alignment_all
## 3 samples were dropped due to being re-done fully in round4
## One was re-done due to
## https://github.com/LieberInstitute/spatialDLPFC/issues/39

## Upload images
sh upload_images.sh

## Locate SpaceRanger _invocation files and upload them. These files include
## * tissue_image_paths: useful for locating the matching input image
## * loupe_alignment_file: useful for locating the matching Loupe json file
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


## Locate Loupe alignment json files
for i in /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/NextSeq/*/_invocation /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/NextSeq/*/*/_invocation
do
    grep "loupe_alignment_file" $i
done

## Upload Loupe alignment json files
sh upload_raw_images_align_json.sh

## Upload this script and companion ones
for i in /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/code/bioturing/*
    do echo $i
    filename=$(basename $i)
    aws --profile bioturing-spatial s3 cp $i s3://bioturing-spatial/spatialDLPFC/code/bioturing/${filename}
done

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
# 2021-11-29 10:52:21    2.0 GiB spatialDLPFC/Images/round2/V10U24-091_B1_Br2720_DLPFC_mid_2.tif
# 2021-11-29 10:52:57    2.1 GiB spatialDLPFC/Images/round2/V10U24-091_C1_Br2720_DLPFC_post_3.tif
# 2021-11-29 10:54:09    1.8 GiB spatialDLPFC/Images/round2/V10U24-092_B1_Br6432_DLPFC_mid_2.tif
# 2021-11-29 10:54:42    1.7 GiB spatialDLPFC/Images/round2/V10U24-092_C1_Br6432_DLPFC_post_3.tif
# 2021-11-29 10:55:39    2.2 GiB spatialDLPFC/Images/round3/1000140_dlpfc_ant_round3_A1.tif
# 2021-11-29 10:56:19    2.1 GiB spatialDLPFC/Images/round3/1000140_dlpfc_ant_round3_B1.tif
# 2021-11-29 10:56:57    2.0 GiB spatialDLPFC/Images/round3/1000140_dlpfc_ant_round3_C1.tif
# 2021-11-29 10:57:33    2.0 GiB spatialDLPFC/Images/round3/1000140_dlpfc_ant_round3_D1.tif
# 2021-11-29 10:58:07    2.9 GiB spatialDLPFC/Images/round3/1000145_dlpfc_mid_round3_A1.tif
# 2021-11-29 10:58:57    2.6 GiB spatialDLPFC/Images/round3/1000145_dlpfc_mid_round3_B1.tif
# 2021-11-29 11:00:24    2.6 GiB spatialDLPFC/Images/round3/1000145_dlpfc_mid_round3_D1.tif
# 2021-11-29 11:01:23    2.5 GiB spatialDLPFC/Images/round3/1000146_dlpfc_post_round3_A1.tif
# 2021-11-29 11:02:04    2.4 GiB spatialDLPFC/Images/round3/1000146_dlpfc_post_round3_B1.tif
# 2021-11-29 11:02:46    2.4 GiB spatialDLPFC/Images/round3/1000146_dlpfc_post_round3_C1.tif
# 2021-11-29 11:03:28    2.4 GiB spatialDLPFC/Images/round3/1000146_dlpfc_post_round3_D1.tif
# 2021-11-29 11:04:57    3.1 GiB spatialDLPFC/Images/round4/V10B01-002_B1_Br6432_ant_DLPFC.tif
# 2021-11-29 11:05:46    3.1 GiB spatialDLPFC/Images/round4/V10B01-002_C1_Br8325_mid_DLPFC.tif
# 2021-12-07 10:18:41    3.5 GiB spatialDLPFC/Images/round4/V10B01-002_D1_Br2720_ant_DLPFC.tif
# 2021-12-07 20:46:26    2.2 KiB spatialDLPFC/SpaceRanger/DLPFC_Br2720_ant_2/_invocation
# 2021-12-07 20:51:55  166 Bytes spatialDLPFC/SpaceRanger/DLPFC_Br2720_ant_2/outs/spatial/scalefactors_json.json
# 2021-12-07 20:44:22    2.2 KiB spatialDLPFC/SpaceRanger/DLPFC_Br2720_mid_manual_alignment/_invocation
# 2021-12-07 20:49:50  148 Bytes spatialDLPFC/SpaceRanger/DLPFC_Br2720_mid_manual_alignment/outs/spatial/scalefactors_json.json
# 2021-12-07 20:44:30    3.3 KiB spatialDLPFC/SpaceRanger/DLPFC_Br2720_post_extra_reads/_invocation
# 2021-12-07 20:49:59  157 Bytes spatialDLPFC/SpaceRanger/DLPFC_Br2720_post_extra_reads/outs/spatial/scalefactors_json.json
# 2021-12-07 20:42:41    1.5 KiB spatialDLPFC/SpaceRanger/DLPFC_Br2743_ant_manual_alignment/_invocation
# 2021-12-07 20:48:09  165 Bytes spatialDLPFC/SpaceRanger/DLPFC_Br2743_ant_manual_alignment/outs/spatial/scalefactors_json.json
# 2021-12-07 20:42:49    2.7 KiB spatialDLPFC/SpaceRanger/DLPFC_Br2743_mid_manual_alignment_extra_reads/_invocation
# 2021-12-07 20:48:17  146 Bytes spatialDLPFC/SpaceRanger/DLPFC_Br2743_mid_manual_alignment_extra_reads/outs/spatial/scalefactors_json.json
# 2021-12-07 20:42:58    1.5 KiB spatialDLPFC/SpaceRanger/DLPFC_Br2743_post_manual_alignment/_invocation
# 2021-12-07 20:48:26  167 Bytes spatialDLPFC/SpaceRanger/DLPFC_Br2743_post_manual_alignment/outs/spatial/scalefactors_json.json
# 2021-12-07 20:43:06    1.5 KiB spatialDLPFC/SpaceRanger/DLPFC_Br3942_ant_manual_alignment/_invocation
# 2021-12-07 20:48:34  157 Bytes spatialDLPFC/SpaceRanger/DLPFC_Br3942_ant_manual_alignment/outs/spatial/scalefactors_json.json
# 2021-12-07 20:43:15    1.5 KiB spatialDLPFC/SpaceRanger/DLPFC_Br3942_mid_manual_alignment/_invocation
# 2021-12-07 20:48:43  155 Bytes spatialDLPFC/SpaceRanger/DLPFC_Br3942_mid_manual_alignment/outs/spatial/scalefactors_json.json
# 2021-12-07 20:43:23    1.5 KiB spatialDLPFC/SpaceRanger/DLPFC_Br3942_post_manual_alignment/_invocation
# 2021-12-07 20:48:51  157 Bytes spatialDLPFC/SpaceRanger/DLPFC_Br3942_post_manual_alignment/outs/spatial/scalefactors_json.json
# 2021-12-07 20:43:31    2.7 KiB spatialDLPFC/SpaceRanger/DLPFC_Br6423_ant_manual_alignment_extra_reads/_invocation
# 2021-12-07 20:48:59  157 Bytes spatialDLPFC/SpaceRanger/DLPFC_Br6423_ant_manual_alignment_extra_reads/outs/spatial/scalefactors_json.json
# 2021-12-07 20:43:40    1.5 KiB spatialDLPFC/SpaceRanger/DLPFC_Br6423_mid_manual_alignment/_invocation
# 2021-12-07 20:49:08  163 Bytes spatialDLPFC/SpaceRanger/DLPFC_Br6423_mid_manual_alignment/outs/spatial/scalefactors_json.json
# 2021-12-07 20:43:48    2.7 KiB spatialDLPFC/SpaceRanger/DLPFC_Br6423_post_extra_reads/_invocation
# 2021-12-07 20:49:17  167 Bytes spatialDLPFC/SpaceRanger/DLPFC_Br6423_post_extra_reads/outs/spatial/scalefactors_json.json
# 2021-12-07 20:46:34    2.2 KiB spatialDLPFC/SpaceRanger/DLPFC_Br6432_ant_2/_invocation
# 2021-12-07 20:52:04  157 Bytes spatialDLPFC/SpaceRanger/DLPFC_Br6432_ant_2/outs/spatial/scalefactors_json.json
# 2021-12-07 20:44:38    2.2 KiB spatialDLPFC/SpaceRanger/DLPFC_Br6432_mid_manual_alignment/_invocation
# 2021-12-07 20:50:07  166 Bytes spatialDLPFC/SpaceRanger/DLPFC_Br6432_mid_manual_alignment/outs/spatial/scalefactors_json.json
# 2021-12-07 20:44:47    2.2 KiB spatialDLPFC/SpaceRanger/DLPFC_Br6432_post_manual_alignment/_invocation
# 2021-12-07 20:50:15  157 Bytes spatialDLPFC/SpaceRanger/DLPFC_Br6432_post_manual_alignment/outs/spatial/scalefactors_json.json
# 2021-12-07 20:44:55    3.3 KiB spatialDLPFC/SpaceRanger/DLPFC_Br6471_ant_manual_alignment_all/_invocation
# 2021-12-07 20:50:24  164 Bytes spatialDLPFC/SpaceRanger/DLPFC_Br6471_ant_manual_alignment_all/outs/spatial/scalefactors_json.json
# 2021-12-07 20:45:03    3.3 KiB spatialDLPFC/SpaceRanger/DLPFC_Br6471_mid_manual_alignment_all/_invocation
# 2021-12-07 20:50:32  166 Bytes spatialDLPFC/SpaceRanger/DLPFC_Br6471_mid_manual_alignment_all/outs/spatial/scalefactors_json.json
# 2021-12-07 20:45:11    3.3 KiB spatialDLPFC/SpaceRanger/DLPFC_Br6471_post_manual_alignment_all/_invocation
# 2021-12-07 20:50:40  166 Bytes spatialDLPFC/SpaceRanger/DLPFC_Br6471_post_manual_alignment_all/outs/spatial/scalefactors_json.json
# 2021-12-07 20:45:20    3.3 KiB spatialDLPFC/SpaceRanger/DLPFC_Br6522_ant_manual_alignment_all/_invocation
# 2021-12-07 20:50:49  156 Bytes spatialDLPFC/SpaceRanger/DLPFC_Br6522_ant_manual_alignment_all/outs/spatial/scalefactors_json.json
# 2021-12-07 20:45:28    3.3 KiB spatialDLPFC/SpaceRanger/DLPFC_Br6522_mid_manual_alignment_all/_invocation
# 2021-12-07 20:50:57  149 Bytes spatialDLPFC/SpaceRanger/DLPFC_Br6522_mid_manual_alignment_all/outs/spatial/scalefactors_json.json
# 2021-12-07 20:45:36    3.3 KiB spatialDLPFC/SpaceRanger/DLPFC_Br6522_post_manual_alignment_all/_invocation
# 2021-12-07 20:51:05  156 Bytes spatialDLPFC/SpaceRanger/DLPFC_Br6522_post_manual_alignment_all/outs/spatial/scalefactors_json.json
# 2021-12-07 20:45:45    3.3 KiB spatialDLPFC/SpaceRanger/DLPFC_Br8325_ant_manual_alignment_all/_invocation
# 2021-12-07 20:51:14  156 Bytes spatialDLPFC/SpaceRanger/DLPFC_Br8325_ant_manual_alignment_all/outs/spatial/scalefactors_json.json
# 2021-12-07 20:46:42    2.2 KiB spatialDLPFC/SpaceRanger/DLPFC_Br8325_mid_2/_invocation
# 2021-12-07 20:52:12  157 Bytes spatialDLPFC/SpaceRanger/DLPFC_Br8325_mid_2/outs/spatial/scalefactors_json.json
# 2021-12-07 20:45:53    3.3 KiB spatialDLPFC/SpaceRanger/DLPFC_Br8325_post_manual_alignment_all/_invocation
# 2021-12-07 20:51:22  148 Bytes spatialDLPFC/SpaceRanger/DLPFC_Br8325_post_manual_alignment_all/outs/spatial/scalefactors_json.json
# 2021-12-07 20:43:56    1.5 KiB spatialDLPFC/SpaceRanger/DLPFC_Br8492_ant_manual_alignment/_invocation
# 2021-12-07 20:49:25  157 Bytes spatialDLPFC/SpaceRanger/DLPFC_Br8492_ant_manual_alignment/outs/spatial/scalefactors_json.json
# 2021-12-07 20:44:05    2.7 KiB spatialDLPFC/SpaceRanger/DLPFC_Br8492_mid_manual_alignment_extra_reads/_invocation
# 2021-12-07 20:49:34  146 Bytes spatialDLPFC/SpaceRanger/DLPFC_Br8492_mid_manual_alignment_extra_reads/outs/spatial/scalefactors_json.json
# 2021-12-07 20:44:13    1.5 KiB spatialDLPFC/SpaceRanger/DLPFC_Br8492_post_manual_alignment/_invocation
# 2021-12-07 20:49:42  158 Bytes spatialDLPFC/SpaceRanger/DLPFC_Br8492_post_manual_alignment/outs/spatial/scalefactors_json.json
# 2021-12-07 20:46:01    3.3 KiB spatialDLPFC/SpaceRanger/DLPFC_Br8667_ant_extra_reads/_invocation
# 2021-12-07 20:51:30  156 Bytes spatialDLPFC/SpaceRanger/DLPFC_Br8667_ant_extra_reads/outs/spatial/scalefactors_json.json
# 2021-12-07 20:46:09    3.3 KiB spatialDLPFC/SpaceRanger/DLPFC_Br8667_mid_manual_alignment_all/_invocation
# 2021-12-07 20:51:39  167 Bytes spatialDLPFC/SpaceRanger/DLPFC_Br8667_mid_manual_alignment_all/outs/spatial/scalefactors_json.json
# 2021-12-07 20:46:18    3.3 KiB spatialDLPFC/SpaceRanger/DLPFC_Br8667_post_manual_alignment_all/_invocation
# 2021-12-07 20:51:47  157 Bytes spatialDLPFC/SpaceRanger/DLPFC_Br8667_post_manual_alignment_all/outs/spatial/scalefactors_json.json
# 2021-12-07 21:43:34   12.6 KiB spatialDLPFC/code/bioturing/check_dates.txt
# 2021-12-07 21:43:42    3.9 KiB spatialDLPFC/code/bioturing/image_files.sh
# 2021-12-07 21:43:51   17.4 KiB spatialDLPFC/code/bioturing/upload_bioturing.sh
# 2021-12-07 21:44:00    7.4 KiB spatialDLPFC/code/bioturing/upload_images.sh
# 2021-12-07 21:44:08    7.8 KiB spatialDLPFC/code/bioturing/upload_raw_images_align_json.sh
# 2021-12-07 21:16:43  556.4 KiB spatialDLPFC/images_raw_align_json/V10B01-002-B1.json
# 2021-12-07 21:16:52  558.8 KiB spatialDLPFC/images_raw_align_json/V10B01-002-C1.json
# 2021-12-07 21:16:34  545.7 KiB spatialDLPFC/images_raw_align_json/V10B01-002-D1.json
# 2021-12-07 21:15:06  564.7 KiB spatialDLPFC/images_raw_align_json/V10B01-052-A1_Br6471_mid.json
# 2021-12-07 21:15:33  554.0 KiB spatialDLPFC/images_raw_align_json/V10B01-052-B1_Br6522_mid.json
# 2021-12-07 21:16:17  560.1 KiB spatialDLPFC/images_raw_align_json/V10B01-052-D1_Br8667_mid.json
# 2021-12-07 21:15:15  557.2 KiB spatialDLPFC/images_raw_align_json/V10B01-053-A1_Br6471_post.json
# 2021-12-07 21:15:41  555.1 KiB spatialDLPFC/images_raw_align_json/V10B01-053-B1_Br6522_post.json
# 2021-12-07 21:15:59  562.3 KiB spatialDLPFC/images_raw_align_json/V10B01-053-C1_Br8325_post.json
# 2021-12-07 21:16:25  562.2 KiB spatialDLPFC/images_raw_align_json/V10B01-053-D1_Br8667_post.json
# 2021-12-07 21:14:11  527.3 KiB spatialDLPFC/images_raw_align_json/V10U24-091-B1_Br2770_mid.json
# 2021-12-07 21:14:23  566.4 KiB spatialDLPFC/images_raw_align_json/V10U24-091-C1_Br2770_post.json
# 2021-12-07 21:14:36  550.1 KiB spatialDLPFC/images_raw_align_json/V10U24-092-B1_Br6432_mid.json
# 2021-12-07 21:14:45  536.6 KiB spatialDLPFC/images_raw_align_json/V10U24-092-C1_Br6432_post.json
# 2021-12-07 21:14:55  540.8 KiB spatialDLPFC/images_raw_align_json/V10U24-094-A1_Br6471_ant.json
# 2021-12-07 21:15:24  560.8 KiB spatialDLPFC/images_raw_align_json/V10U24-094-B1_Br6522_ant.json
# 2021-12-07 21:15:50  550.3 KiB spatialDLPFC/images_raw_align_json/V10U24-094-C1_Br8325_ant.json
# 2021-12-07 21:16:08  552.0 KiB spatialDLPFC/images_raw_align_json/V10U24-094-D1_Br8667_ant.json
# 2021-12-07 21:12:35  555.6 KiB spatialDLPFC/images_raw_align_json/V19B23-073-A1_9_Br2743_post_manual_alignment.json
# 2021-12-07 21:13:02  562.2 KiB spatialDLPFC/images_raw_align_json/V19B23-073-B1_10_Br3942_post_manual_alignment.json
# 2021-12-07 21:13:31  550.4 KiB spatialDLPFC/images_raw_align_json/V19B23-073-C1_11_Br6423_post_manual_alignment.json
# 2021-12-07 21:13:57  560.8 KiB spatialDLPFC/images_raw_align_json/V19B23-073-D1_12_Br8492_post_manual_alignment.json
# 2021-12-07 21:12:27  560.3 KiB spatialDLPFC/images_raw_align_json/V19B23-074-A1_5_Br2743_mid_manual_alignment.json
# 2021-12-07 21:12:53  555.9 KiB spatialDLPFC/images_raw_align_json/V19B23-074-B1_6_Br3942_mid_manual_alignment.json
# 2021-12-07 21:13:23  556.1 KiB spatialDLPFC/images_raw_align_json/V19B23-074-C1_7_Br6432_mid_manual_alignment.json
# 2021-12-07 21:13:49  562.8 KiB spatialDLPFC/images_raw_align_json/V19B23-074-D1_8_Br8492_mid_manual_alignment.json
# 2021-12-07 21:12:13  554.1 KiB spatialDLPFC/images_raw_align_json/V19B23-075-A1_1_Br2743_ant_manual_alignment.json
# 2021-12-07 21:12:44  555.6 KiB spatialDLPFC/images_raw_align_json/V19B23-075-B1_2_Br3942_ant_manual_alignment.json
# 2021-12-07 21:13:14  555.4 KiB spatialDLPFC/images_raw_align_json/V19B23-075-C1_3_Br6423_ant_manual_alignment.json
# 2021-12-07 21:13:40  567.5 KiB spatialDLPFC/images_raw_align_json/V19B23-075-D1_4_Br8492_ant_manual_alignment.json
# 2021-11-29 10:21:12    1.5 GiB spatialDLPFC/spe_merged_final.Rdata

## Re-upload this script after updating the ls output and re-saving this
## script
aws --profile bioturing-spatial s3 cp /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/code/bioturing/upload_bioturing.sh s3://bioturing-spatial/spatialDLPFC/code/bioturing/upload_bioturing.sh

