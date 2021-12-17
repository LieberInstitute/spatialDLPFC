#!/bin/bash

## Load aws module
module load aws/2.2.35

## Upload main SpatialExperiment object
aws --profile bioturing-spatial s3 cp /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/spe_final.Rdata s3://bioturing-spatial/spatialDLPFC/spe_2021-12-17.Rdata

## Locate images using the SpaceRanger output
for i in /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rerun_spaceranger/*/_invocation
do
    grep "tissue_image_paths" $i
done

## Upload images
sh upload_images.sh

## Locate SpaceRanger _invocation files and upload them. These files include
## * tissue_image_paths: useful for locating the matching input image
## * loupe_alignment_file: useful for locating the matching Loupe json file
for i in /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rerun_spaceranger/*/_invocation
    do echo $i
    dir1=$(dirname $i)
    dirname=$(basename $dir1)
    filename=$(basename $i)
    #echo "$i s3://bioturing-spatial/spatialDLPFC/SpaceRanger/${dirname}/${filename}"
    aws --profile bioturing-spatial s3 cp $i s3://bioturing-spatial/spatialDLPFC/SpaceRanger/${dirname}/${filename}
done

## Upload scalefactors_json.json files
for i in /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rerun_spaceranger/*/outs/spatial/scalefactors_json.json
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
for i in /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rerun_spaceranger/*/_invocation
do
    grep "loupe_alignment_file" $i
done

## Upload Loupe alignment json files
sh upload_raw_images_align_json.sh

## Upload this script and companion ones
for i in /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/code/bioturing/*.sh
    do echo $i
    filename=$(basename $i)
    aws --profile bioturing-spatial s3 cp $i s3://bioturing-spatial/spatialDLPFC/code/bioturing/${filename}
done

## Upload Visium_IF mock images
aws --profile bioturing-spatial s3 cp /dcs04/lieber/lcolladotor/deconvolution_LIBD4030/DLPFC_Visium_IF/raw-data/Images/20211124_VIF_DLPFC_Mock_Polaris.tar.gz s3://bioturing-spatial/20211124_VIF_DLPFC_Mock/

## List files on the S3 bucket
aws --profile bioturing-spatial s3 ls --recursive --human-readable s3://bioturing-spatial/
