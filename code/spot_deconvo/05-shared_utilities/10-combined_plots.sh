#!/bin/bash
#$ -cwd
#$ -N "combined_plots"
#$ -o ../../../processed-data/spot_deconvo/05-shared_utilities/logs/10-combined_plots_nonIF.log
#$ -e ../../../processed-data/spot_deconvo/05-shared_utilities/logs/10-combined_plots_nonIF.log
#$ -l bluejay,mf=20G,h_vmem=20G

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

module load conda_R/4.2.x
Rscript 10-combined_plots.R

echo "**** Job ends ****"
date
