#!/bin/bash
#$ -cwd
#$ -l mem_free=5G,h_vmem=5G
#$ -N "download"
#$ -o ../../../processed-data/synapse_upload/03-download/01-download.log
#$ -e ../../../processed-data/synapse_upload/03-download/01-download.log

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

module load synapse/2.6.0
python 01-download.py

echo "**** Job ends ****"
date
