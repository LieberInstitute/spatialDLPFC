#!/bin/bash
#$ -cwd
#$ -N "check_confidence"
#$ -o ../../../processed-data/spot_deconvo/02-cellpose/12-check_confidence.log
#$ -e ../../../processed-data/spot_deconvo/02-cellpose/12-check_confidence.log
#$ -l mf=10G,h_vmem=10G

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

module load cellpose/2.0
python 12-check_confidence.py

echo "**** Job ends ****"
date
