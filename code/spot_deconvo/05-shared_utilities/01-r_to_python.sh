#!/bin/bash
#$ -cwd
#$ -l mem_free=270G,h_vmem=270G,h_fsize=100G
#$ -N "r_to_python"
#$ -o ../../../processed-data/spot_deconvo/05-shared_utilities/logs/01-r_to_python.log
#$ -e ../../../processed-data/spot_deconvo/05-shared_utilities/logs/01-r_to_python.log

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

module load conda_R/devel
Rscript 01-r_to_python.R

echo "**** Job ends ****"
date
