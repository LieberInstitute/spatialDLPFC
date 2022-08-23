#!/bin/bash
#$ -cwd
#$ -N "prepare_anndatas_nonIF"
#$ -o ../../../processed-data/spot_deconvo/01-tangram/logs/02-prepare_anndatas_nonIF_$TASK_ID.log
#$ -e ../../../processed-data/spot_deconvo/01-tangram/logs/02-prepare_anndatas_nonIF_$TASK_ID.log
#$ -l bluejay,mf=120G,h_vmem=120G,h_fsize=50G
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
python 02-prepare_anndatas_nonIF.py

echo "**** Job ends ****"
date
