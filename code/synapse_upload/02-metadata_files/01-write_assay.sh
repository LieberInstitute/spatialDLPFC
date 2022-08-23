#!/bin/bash
#$ -cwd
#$ -l mem_free=3G,h_vmem=3G
#$ -N "write_assay"
#$ -o ../../../processed-data/synapse_upload/02-metadata_files/01-write_assay.log
#$ -e ../../../processed-data/synapse_upload/02-metadata_files/01-write_assay.log

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

module load conda_R/devel
Rscript 01-write_assay.R

echo "**** Job ends ****"
date
