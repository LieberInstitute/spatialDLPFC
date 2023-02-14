#!/bin/bash
#$ -cwd
#$ -l mem_free=10G,h_vmem=10G
#$ -N "write_individual"
#$ -o ../../../processed-data/synapse_upload/02-metadata_files/03-write_individual.log
#$ -e ../../../processed-data/synapse_upload/02-metadata_files/03-write_individual.log

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

module load conda_R/4.2.x
Rscript 03-write_individual.R

echo "**** Job ends ****"
date
