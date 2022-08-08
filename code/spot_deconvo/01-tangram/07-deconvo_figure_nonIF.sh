#!/bin/bash
#$ -cwd
#$ -N "deconvo_figure_nonIF"
#$ -o ../../../processed-data/spot_deconvo/01-tangram/logs/07-deconvo_figure_nonIF_$TASK_ID.log
#$ -e ../../../processed-data/spot_deconvo/01-tangram/logs/07-deconvo_figure_nonIF_$TASK_ID.log
#$ -l mf=10G,h_vmem=10G
#$ -t 1
#$ -tc 1

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

module load conda_R/devel
Rscript 07-deconvo_figure_nonIF.R

echo "**** Job ends ****"
date
