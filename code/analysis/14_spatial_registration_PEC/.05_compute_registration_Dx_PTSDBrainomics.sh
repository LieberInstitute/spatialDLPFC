#!/bin/bash
#$ -cwd
#$ -l mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -N compute_registration_Dx_PTSDBrainomics
#$ -o logs/05_compute_registration_Dx_PTSDBrainomics.txt
#$ -e logs/05_compute_registration_Dx_PTSDBrainomics.txt
#$ -hold_jid pseudobulk_data_PTSDBrainomics
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
module load conda_R/4.2.x

## List current modules for reproducibility
module list

## Sumbit Rscript
Rscript 05_compute_registration_Dx.R PTSDBrainomics

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/


