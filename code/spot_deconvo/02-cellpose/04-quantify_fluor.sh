#!/bin/bash
#$ -cwd
#$ -N "quantify_fluor"
#$ -o ../../../processed-data/spot_deconvo/02-cellpose/04-quantify_fluor_$TASK_ID.log
#$ -e ../../../processed-data/spot_deconvo/02-cellpose/04-quantify_fluor_$TASK_ID.log
#$ -l mf=20G,h_vmem=20G
#$ -t 2-4
#$ -tc 2

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

module load cellpose/2.0
python 04-quantify_fluor.py

echo "**** Job ends ****"
date
