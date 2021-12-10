#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=20G,h_vmem=20G,h_fsize=100G
#$ -N "sce_to_anndata_dlpfc_03"
#$ -o ../../processed-data/03_nn_run/logs/03_sce_to_anndata_dlpfc.log
#$ -e ../../processed-data/03_nn_run/logs/03_sce_to_anndata_dlpfc.log

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

module load conda_R/4.1.x
Rscript 03_sce_to_anndata_dlpfc.R

echo "**** Job ends ****"
date
