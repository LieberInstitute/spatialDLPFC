#!/bin/bash
#$ -cwd
#$ -N "export_results"
#$ -o ../../../processed-data/spot_deconvo/05-shared_utilities/logs/11-export_results.log
#$ -e ../../../processed-data/spot_deconvo/05-shared_utilities/logs/11-export_results.log
#$ -l mf=10G,h_vmem=10G

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

module load conda_R/4.2.x
Rscript 11-export_results.R

echo "**** Job ends ****"
date
