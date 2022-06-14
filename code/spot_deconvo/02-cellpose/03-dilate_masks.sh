#!/bin/bash
#$ -cwd
#$ -N "dilate_masks"
#$ -o ../../../processed-data/spot_deconvo/02-cellpose/03-dilate_masks.log
#$ -e ../../../processed-data/spot_deconvo/02-cellpose/03-dilate_masks.log
#$ -l bluejay,mf=20G,h_vmem=20G

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

module load cellpose/2.0
python 03-dilate_masks.py

echo "**** Job ends ****"
date
