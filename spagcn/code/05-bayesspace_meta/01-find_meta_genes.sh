#!/bin/bash
#$ -cwd
#$ -N "find_meta_genes_01"
#$ -o ../../processed-data/05-bayesspace_meta/01-find_meta_genes_$TASK_ID.log
#$ -e ../../processed-data/05-bayesspace_meta/01-find_meta_genes_$TASK_ID.log
#$ -l bluejay,mf=20G,h_vmem=20G
#$ -t 6-30
#$ -tc 10

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

#   The script internally reads the 'SGE_TASK_ID' environment variable
module load spagcn/1.2.0
python 01-find_meta_genes.py

echo "**** Job ends ****"
date
