#!/bin/bash
#$ -cwd
#$ -l mem_free=3G,h_vmem=3G
#$ -N "prepare_fastq"
#$ -o ../../../processed-data/synapse_upload/01-prepare_fastq/01-prepare_fastq.log
#$ -e ../../../processed-data/synapse_upload/01-prepare_fastq/01-prepare_fastq.log

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

module load conda_R/4.3
Rscript 01-prepare_fastq.R

echo "**** Job ends ****"
date
