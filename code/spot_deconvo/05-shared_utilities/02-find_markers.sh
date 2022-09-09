#!/bin/bash
#$ -cwd
#$ -l mem_free=20G,h_vmem=20G
#$ -N "find_markers"
#$ -o ../../../processed-data/spot_deconvo/05-shared_utilities/logs/02-find_markers_broad.log
#$ -e ../../../processed-data/spot_deconvo/05-shared_utilities/logs/02-find_markers_broad.log

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

module load conda_R/4.1.x
Rscript 02-find_markers.R

echo "**** Job ends ****"
date
