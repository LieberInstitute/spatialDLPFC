#!/bin/bash
#$ -cwd
#$ -l mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -N nnSVG_pairwise_loop_7v12
#$ -o logs/01_nnSVG_pairwise_loop_7v12.$TASK_ID.txt
#$ -e logs/01_nnSVG_pairwise_loop_7v12.$TASK_ID.txt
#$ -m e
#$ -t 1-30
#$ -tc 10

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
Rscript 01_nnSVG_pairwise.R 7v12 $SGE_TASK_ID

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/


