#!/bin/bash
#$ -cwd
#$ -N "sample_table"
#$ -o ../../../processed-data/spot_deconvo/01-tangram/logs/08-sample_table.log
#$ -e ../../../processed-data/spot_deconvo/01-tangram/logs/08-sample_table.log
#$ -l mf=10G,h_vmem=10G

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

module load conda_R/devel
Rscript 08-sample_table.R

echo "**** Job ends ****"
date
