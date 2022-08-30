#!/bin/bash
#$ -cwd
#$ -N "deconvo_figure_IF"
#$ -o ../../../processed-data/spot_deconvo/05-shared_utilities/logs/03-deconvo_figure_IF_$TASK_ID.log
#$ -e ../../../processed-data/spot_deconvo/05-shared_utilities/logs/03-deconvo_figure_IF_$TASK_ID.log
#$ -l mf=10G,h_vmem=10G
#$ -t 1-4
#$ -tc 4

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

module load conda_R/devel
Rscript 03-deconvo_figure_IF.R

echo "**** Job ends ****"
date
