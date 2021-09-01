#!/bin/bash
#$ -pe local 10
#$ -l mem_free=20G,h_vmem=20G,h_stack=256M,h_fsize=100G
#$ -o /dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/spotspotcheck-master/main2.txt
#$ -e /dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/spotspotcheck-master/main2.txt


echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

export TZ=America/New_York
matlab -nodesktop -nosplash -r "addpath(genpath('/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/spotspotcheck-master')); main2; exit;" 

echo "**** Job ends ****"

date
