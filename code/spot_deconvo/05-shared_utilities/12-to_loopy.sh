#!/bin/bash
#$ -cwd
#$ -N "to_loopy"
#$ -o ../../../processed-data/spot_deconvo/05-shared_utilities/logs/12-to_loopy.log
#$ -e ../../../processed-data/spot_deconvo/05-shared_utilities/logs/12-to_loopy.log
#$ -l mf=10G,h_vmem=10G

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

module load conda_R/4.2.x
Rscript 12-to_loopy.R

echo "**** Job ends ****"
date
