#!/bin/bash
#$ -cwd
#$ -l mem_free=25G,h_vmem=25G,h_fsize=100G
#$ -N compute_registration_stats_SZBDMulti-Seq
#$ -o logs/02_compute_registration_stats_SZBDMulti-Seq.txt
#$ -e logs/02_compute_registration_stats_SZBDMulti-Seq.txt
#$ -hold_jid pseudobulk_data_SZBDMulti-Seq
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
Rscript 02_compute_registration_stats.R SZBDMulti-Seq

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/


