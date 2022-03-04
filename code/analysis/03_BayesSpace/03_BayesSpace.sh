#!/bin/bash
#$ -cwd
#$ -l mem_free=80G,h_vmem=80G,h_fsize=100G
#$ -N bayesSpace_k_grid
#$ -o logs/bayesSpace_k_grid.$TASK_ID.txt
#$ -e logs/bayesSpace_k_grid.$TASK_ID.txt
#$ -m e
#$ -t 2-15
#$ -tc 20

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

## Load the R module (absent since the JHPCE upgrade to CentOS v7)
module load conda_R/devel

## List current modules for reproducibility
module list

## Edit with your job command
Rscript 03_BayesSpace.R

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/
