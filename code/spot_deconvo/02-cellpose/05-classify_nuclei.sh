#!/bin/bash
#$ -cwd
#$ -N "classify_nuclei"
#$ -o ../../../processed-data/spot_deconvo/02-cellpose/05-classify_nuclei_$TASK_ID.log
#$ -e ../../../processed-data/spot_deconvo/02-cellpose/05-classify_nuclei_$TASK_ID.log
#$ -l mf=20G,h_vmem=20G
#$ -t 1-4
#$ -tc 4

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

module load cellpose/2.0
python 05-classify_nuclei.py

echo "**** Job ends ****"
date
