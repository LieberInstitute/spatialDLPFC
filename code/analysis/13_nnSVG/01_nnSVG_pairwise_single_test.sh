#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=25G,h_vmem=25G,h_fsize=100G
#$ -N nnSVG_pairwise_single_test
#$ -o logs/01_nnSVG_pairwise_single_test.txt
#$ -e logs/01_nnSVG_pairwise_single_test.txt
#$ -m e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

## Load the R module (absent since the JHPCE upgrade to CentOS v7)
module load conda_R/4.2

## List current modules for reproducibility
module list

## Edit with your job command
Rscript 01_nnSVG_pairwise.R 4v16 1

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/
