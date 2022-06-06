#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=70G,h_vmem=70G,h_fsize=100G
#$ -N "r_to_python"
#$ -o ../../processed-data/01-nn_run/00-r_to_python.log
#$ -e ../../processed-data/01-nn_run/00-r_to_python.log

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

module load conda_R/devel
Rscript 00-r_to_python.R

echo "**** Job ends ****"
date
