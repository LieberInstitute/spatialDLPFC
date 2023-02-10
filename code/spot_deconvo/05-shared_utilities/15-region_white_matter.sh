#!/bin/bash
#$ -cwd
#$ -N "region_white_matter"
#$ -o ../../../processed-data/spot_deconvo/05-shared_utilities/logs/15-region_white_matter.log
#$ -e ../../../processed-data/spot_deconvo/05-shared_utilities/logs/15-region_white_matter.log
#$ -l mf=15G,h_vmem=15G

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

module load conda_R/4.2.x
Rscript 15-region_white_matter.R

echo "**** Job ends ****"
date
