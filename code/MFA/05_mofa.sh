#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=10G
#SBATCH --job-name=05_mofa
#SBATCH -c 1
#SBATCH -t 1-00:00:00
#SBATCH -o /dev/null
#SBATCH -e /dev/null
#SBATCH --array=9-17%20

## Define loops and appropriately subset each variable for the array task ID
all_dataset=(DLPFC combined)
dataset=${all_dataset[$(( $SLURM_ARRAY_TASK_ID / 9 % 2 ))]}

all_num_factors=(2 3 4 5 6 7 8 9 10)
num_factors=${all_num_factors[$(( $SLURM_ARRAY_TASK_ID / 1 % 9 ))]}

## Explicitly pipe script output to a log
log_path=../../processed-data/MFA/logs/05_mofa_${dataset}_${num_factors}_${SLURM_ARRAY_TASK_ID}.txt

{
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

## Edit with your job command
Rscript 05_mofa.R --dataset ${dataset} --num_factors ${num_factors}

echo "**** Job ends ****"
date

} > $log_path 2>&1

## This script was made using slurmjobs version 1.3.0
## available from http://research.libd.org/slurmjobs/
