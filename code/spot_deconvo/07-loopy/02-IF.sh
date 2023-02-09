#!/bin/bash
#$ -cwd
#$ -N "IF"
#$ -o ../../../processed-data/spot_deconvo/07-loopy/02-IF_$TASK_ID.log
#$ -e ../../../processed-data/spot_deconvo/07-loopy/02-IF_$TASK_ID.log
#$ -l mf=10G,h_vmem=10G
#$ -t 2-4
#$ -tc 3

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

module load loopy/1.0.0-next.8
python 02-IF.py

echo "**** Job ends ****"
date
