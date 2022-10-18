#!/bin/bash
#$ -cwd
#$ -N "explore_models"
#$ -o ../../../processed-data/spot_deconvo/02-cellpose/10-explore_models.log
#$ -e ../../../processed-data/spot_deconvo/02-cellpose/10-explore_models.log
#$ -l mf=10G,h_vmem=10G

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

module load cellpose/2.0
python 10-explore_models.py

echo "**** Job ends ****"
date
