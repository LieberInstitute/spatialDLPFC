#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=50G,h_vmem=50G,h_fsize=100G
#$ -N "sce-to-anndata"
#$ -j y
#$ -o ../../processed-data/03_nn_run/logs/sce-2-anndata.log
#$ -e ../../processed-data/03_nn_run/logs/sce-2-anndata.log

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

## Load dependencies
module load conda_R/4.0.x

## List current modules
module list

Rscript 01_sce_to_anndata.R

echo "**** Job ends ****"
date
