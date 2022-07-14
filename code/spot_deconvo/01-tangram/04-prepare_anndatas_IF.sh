#!/bin/bash
#$ -cwd
#$ -N "prepare_anndatas_IF"
#$ -o ../../../processed-data/spot_deconvo/01-tangram/logs/04-prepare_anndatas_IF_$TASK_ID.log
#$ -e ../../../processed-data/spot_deconvo/01-tangram/logs/04-prepare_anndatas_IF_$TASK_ID.log
#$ -l bluejay,mf=80G,h_vmem=80G,h_fsize=50G
#$ -t 2-4
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
python 04-prepare_anndatas_IF.py

echo "**** Job ends ****"
date
