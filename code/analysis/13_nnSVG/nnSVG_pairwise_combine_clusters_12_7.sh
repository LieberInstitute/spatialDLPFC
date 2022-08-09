#!/bin/bash
#$ -cwd
#$ -l mem_free=15G,h_vmem=15G,h_fsize=200G
#$ -N nnSVG_pairwise_combine_12_7
#$ -o logs/nnSVG_pairwise_combine_clusters_12_7.$TASK_ID.txt
#$ -e logs/nnSVG_pairwise_combine_clusters_12_7.$TASK_ID.txt
#$ -m e
#$ -pe local 8


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
Rscript nnSVG_pairwise_combine_clusters_12_7.R

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/
