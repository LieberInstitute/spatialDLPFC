#!/bin/bash
#$ -cwd
#$ -N "dlpfc_prepare_anndata"
#$ -o ../../processed-data/03_nn_run/logs/05_dlpfc_prepare_anndata_$TASK_ID.log
#$ -e ../../processed-data/03_nn_run/logs/05_dlpfc_prepare_anndata_$TASK_ID.log
#$ -l bluejay,mf=10G,h_vmem=10G
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
python 05_dlpfc_prepare_anndata.py

echo "**** Job ends ****"
date
