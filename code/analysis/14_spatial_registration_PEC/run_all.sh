#!/bin/bash
#$ -cwd
#$ -l mem_free=2G,h_vmem=2G,h_fsize=100G
#$ -N run_all
#$ -o logs/run_all.txt
#$ -e logs/run_all.txt
#$ -m e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

## List current modules for reproducibility
module list

# CMC
rm logs/01_pseudobulk_data_CMC.txt
qsub 01_pseudobulk_data_CMC.sh

rm logs/02_compute_registration_stats_CMC.txt
qsub 02_compute_registration_stats_CMC.sh

# DevBrain
rm logs/01_pseudobulk_data_DevBrain.txt
qsub 01_pseudobulk_data_DevBrain.sh

rm logs/02_compute_registration_stats_DevBrain.txt
qsub 02_compute_registration_stats_DevBrain.sh

# IsoHuB
rm logs/01_pseudobulk_data_IsoHuB.txt
qsub 01_pseudobulk_data_IsoHuB.sh

rm logs/02_compute_registration_stats_IsoHuB.txt
qsub 02_compute_registration_stats_IsoHuB.sh

# SZBD
rm logs/01_pseudobulk_data_SZBD.txt
qsub 01_pseudobulk_data_SZBD.sh

rm logs/02_compute_registration_stats_SZBD.txt
qsub 02_compute_registration_stats_SZBD.sh

# UCLA
rm logs/01_pseudobulk_data_UCLA.txt
qsub 01_pseudobulk_data_UCLA.sh

rm logs/02_compute_registration_stats_UCLA.txt
qsub 02_compute_registration_stats_UCLA.sh

## Compile and Summarize
rm logs/03_correlate_spatial.txt
qsub 03_correlate_spatial

echo "**** Job ends ****"
date