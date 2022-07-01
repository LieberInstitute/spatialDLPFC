#!/bin/bash
#$ -cwd
#$ -N "count_cells_all_channels"
#$ -o ../../../processed-data/spot_deconvo/02-cellpose/04-count_cells_all_channels_gfap.log
#$ -e ../../../processed-data/spot_deconvo/02-cellpose/04-count_cells_all_channels_gfap.log
#$ -l bluejay,mf=20G,h_vmem=20G
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

module load cellpose/2.0
python 04-count_cells_all_channels.py

echo "**** Job ends ****"
date
