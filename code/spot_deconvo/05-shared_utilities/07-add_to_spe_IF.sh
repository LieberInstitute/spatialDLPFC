#!/bin/bash
#$ -cwd
#$ -N "add_to_spe_IF"
#$ -o ../../../processed-data/spot_deconvo/05-shared_utilities/logs/07-add_to_spe_IF.log
#$ -e ../../../processed-data/spot_deconvo/05-shared_utilities/logs/07-add_to_spe_IF.log
#$ -l mf=10G,h_vmem=10G

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

module load conda_R/4.2.x
Rscript 07-add_to_spe_IF.R

echo "**** Job ends ****"
date
