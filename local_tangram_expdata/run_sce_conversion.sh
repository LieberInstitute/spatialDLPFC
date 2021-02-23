#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=30G,h_vmem=30G,h_fsize=100G
#$ -pe local 5
#$ -N "sce-to-anndata"
#$ -j y
#$ -o logs/sce-2-anndata.log

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

Rscript sce_to_anndata.R "JHPCE"

echo "**** Job ends ****"
date
