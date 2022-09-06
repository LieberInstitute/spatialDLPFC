#!/bin/bash
#$ -cwd
#$ -N "IF"
#$ -o ../../../processed-data/spot_deconvo/04-spotlight/01-IF_layer.log
#$ -e ../../../processed-data/spot_deconvo/04-spotlight/01-IF_layer.log
#$ -l mf=20G,h_vmem=20G,h_fsize=50G

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

module load conda_R/devel
Rscript 01-IF.R

echo "**** Job ends ****"
date
