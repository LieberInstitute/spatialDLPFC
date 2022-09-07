#!/bin/bash
#$ -cwd
#$ -N "prepare_loopy"
#$ -o ../../../processed-data/spot_deconvo/02-cellpose/09-prepare_loopy_$TASK_ID.log
#$ -e ../../../processed-data/spot_deconvo/02-cellpose/09-prepare_loopy_$TASK_ID.log
#$ -l mf=5G,h_vmem=5G
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
python 09-prepare_loopy.py

echo "**** Job ends ****"
date
