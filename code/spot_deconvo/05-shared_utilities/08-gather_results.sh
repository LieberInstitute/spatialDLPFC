#!/bin/bash
#$ -cwd
#$ -N "gather_results"
#$ -o ../../../processed-data/spot_deconvo/05-shared_utilities/logs/08-gather_results.log
#$ -e ../../../processed-data/spot_deconvo/05-shared_utilities/logs/08-gather_results.log
#$ -l mf=5G,h_vmem=5G

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

module load conda_R/4.2.x
Rscript 08-gather_results.R

echo "**** Job ends ****"
date
