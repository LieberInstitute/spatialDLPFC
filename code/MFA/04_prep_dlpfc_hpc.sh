#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=10G
#SBATCH --job-name=04_prep_dlpfc_hpc
#SBATCH -c 1
#SBATCH -t 1-0:00:00
#SBATCH -o ../../processed-data/MFA/logs/04_prep_dlpfc_hpc.txt
#SBATCH -e ../../processed-data/MFA/logs/04_prep_dlpfc_hpc.txt

set -e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

## Load the R module
module load conda_R/4.5

## List current modules for reproducibility
module list

Rscript 04_prep_dlpfc_hpc.R

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.2.4
## available from http://research.libd.org/slurmjobs/
