#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=20G
#SBATCH --job-name=07_coherence
#SBATCH -c 1
#SBATCH -t 1-00:00:00
#SBATCH -o /dev/null
#SBATCH -e /dev/null
#SBATCH --array=1-576%20

## Define loops and appropriately subset each variable for the array task ID
all_DLPFC_clus=(DLPFC_Sp09D01 DLPFC_Sp09D02 DLPFC_Sp09D03 DLPFC_Sp09D04 DLPFC_Sp09D05 DLPFC_Sp09D06 DLPFC_Sp09D07 DLPFC_Sp09D08 DLPFC_Sp09D09)
DLPFC_clus=${all_DLPFC_clus[$(( $SLURM_ARRAY_TASK_ID / 64 % 9 ))]}

all_HPC_clus=(HPC_visium_CA1 HPC_visium_CA2.4 HPC_visium_Choroid HPC_visium_GABA HPC_visium_GCL HPC_visium_ML HPC_visium_RHP HPC_visium_SL.SR HPC_visium_SLM.SGZ HPC_visium_SR.SLM HPC_visium_SUB HPC_visium_SUB.RHP HPC_visium_Vascular HPC_visium_WM.1 HPC_visium_WM.2 HPC_visium_WM.3)
HPC_clus=${all_HPC_clus[$(( $SLURM_ARRAY_TASK_ID / 4 % 16 ))]}

all_filter_by_region=(FALSE TRUE)
filter_by_region=${all_filter_by_region[$(( $SLURM_ARRAY_TASK_ID / 2 % 2 ))]}

all_clean_expression=(FALSE TRUE)
clean_expression=${all_clean_expression[$(( $SLURM_ARRAY_TASK_ID / 1 % 2 ))]}

## Explicitly pipe script output to a log
log_path=../../processed-data/MFA/logs/07_coherence_${DLPFC_clus}_${HPC_clus}_${filter_by_region}_${clean_expression}_${SLURM_ARRAY_TASK_ID}.txt

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
Rscript 07_coherence.R --DLPFC_clus ${DLPFC_clus} --HPC_clus ${HPC_clus} --filter_by_region ${filter_by_region} --clean_expression ${clean_expression}

echo "**** Job ends ****"
date

} > $log_path 2>&1

## This script was made using slurmjobs version 1.3.0
## available from http://research.libd.org/slurmjobs/

