#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=70G,h_vmem=70G,h_fsize=100G
#$ -N "r_to_python_nonIF"
#$ -o ../../../processed-data/spot_deconvo/01-tangram/logs/01-r_to_python_nonIF.log
#$ -e ../../../processed-data/spot_deconvo/01-tangram/logs/01-r_to_python_nonIF.log

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

module load conda_R/devel
Rscript 01-r_to_python_nonIF.R

echo "**** Job ends ****"
date
