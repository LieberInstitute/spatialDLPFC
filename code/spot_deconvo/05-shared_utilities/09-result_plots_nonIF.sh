#!/bin/bash
#$ -cwd
#$ -N "result_plots_nonIF"
#$ -o ../../../processed-data/spot_deconvo/05-shared_utilities/logs/09-result_plots_nonIF_broad.log
#$ -e ../../../processed-data/spot_deconvo/05-shared_utilities/logs/09-result_plots_nonIF_broad.log
#$ -l mf=20G,h_vmem=20G

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

module load conda_R/4.2.x
Rscript 09-result_plots_nonIF.R

echo "**** Job ends ****"
date
