#!/bin/bash
#$ -cwd
#$ -N "pan_tangram"
#$ -j y
#$ -o ../../processed-data/03_nn_run/logs/02_pan_tangram_$TASK_ID.log
#$ -e ../../processed-data/03_nn_run/logs/02_pan_tangram_$TASK_ID.log
#$ -l gpu,mf=64G,h_vmem=64G
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

module load tangram/1.0.0

python 02_pan_tangram.py -i $SGE_TASK_ID

echo "**** Job ends ****"
date
