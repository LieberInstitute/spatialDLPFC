#!/bin/bash
#$ -cwd
#$ -N "nonIF"
#$ -o ../../../processed-data/spot_deconvo/04-spotlight/02-nonIF_broad.log
#$ -e ../../../processed-data/spot_deconvo/04-spotlight/02-nonIF_broad.log
#$ -l mf=20G,h_vmem=20G,h_fsize=50G

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

module load conda_R/devel
Rscript 02-nonIF.R

echo "**** Job ends ****"
date
