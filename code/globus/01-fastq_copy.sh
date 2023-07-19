#!/bin/bash
#$ -cwd
#$ -l rnet,mem_free=3G,h_vmem=3G,h_fsize=300G
#$ -N "fastq_copy"
#$ -o ../../processed-data/globus/01-fastq_copy.log
#$ -e ../../processed-data/globus/01-fastq_copy.log

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

module load conda_R/4.3
Rscript 01-fastq_copy.R

echo "**** Job ends ****"
date
