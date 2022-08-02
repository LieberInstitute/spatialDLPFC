#!/bin/bash
#$ -cwd
#$ -l mem_free=10G,h_vmem=10G
#$ -N "write_manifest"
#$ -o ../../../processed-data/synapse_upload/02-metadata_files/04-write_manifest.log
#$ -e ../../../processed-data/synapse_upload/02-metadata_files/04-write_manifest.log

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

module load conda_R/devel
Rscript 04-write_manifest.R

echo "**** Job ends ****"
date
