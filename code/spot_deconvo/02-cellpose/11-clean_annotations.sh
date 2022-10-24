#!/bin/bash
#$ -cwd
#$ -N "clean_annotations"
#$ -o ../../../processed-data/spot_deconvo/02-cellpose/11-clean_annotations.log
#$ -e ../../../processed-data/spot_deconvo/02-cellpose/11-clean_annotations.log
#$ -l mf=10G,h_vmem=10G

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

module load cellpose/2.0
python 11-clean_annotations.py

echo "**** Job ends ****"
date
