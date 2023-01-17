#!/bin/bash
#$ -cwd
#$ -N "compare_settings"
#$ -o ../../../processed-data/spot_deconvo/04-spotlight/04-compare_settings_layer.log
#$ -e ../../../processed-data/spot_deconvo/04-spotlight/04-compare_settings_layer.log
#$ -l mf=10G,h_vmem=10G

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

module load conda_R/4.2.x
Rscript 04-compare_settings.R

echo "**** Job ends ****"
date
