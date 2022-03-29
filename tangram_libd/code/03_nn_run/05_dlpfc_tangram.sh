#!/bin/bash
#$ -cwd
#$ -N "dlpfc_tangram_05"
#$ -o ../../processed-data/03_nn_run/logs/05_dlpfc_tangram_$TASK_ID.log
#$ -e ../../processed-data/03_nn_run/logs/05_dlpfc_tangram_$TASK_ID.log
#$ -l caracol,mf=64G,h_vmem=64G
#$ -t 1
#$ -tc 1

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

module load tangram/1.0.2

python 05_dlpfc_tangram.py

echo "**** Job ends ****"
date
