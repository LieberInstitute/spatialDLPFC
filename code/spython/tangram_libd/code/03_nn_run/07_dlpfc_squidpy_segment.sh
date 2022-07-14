#!/bin/bash
#$ -cwd
#$ -N "dlpfc_squidpy_segment"
#$ -o ../../processed-data/03_nn_run/logs/07_dlpfc_squidpy_segment_$TASK_ID.log
#$ -e ../../processed-data/03_nn_run/logs/07_dlpfc_squidpy_segment_$TASK_ID.log
#$ -l bluejay,mf=20G,h_vmem=20G
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
python 07_dlpfc_squidpy_segment.py

echo "**** Job ends ****"
date
