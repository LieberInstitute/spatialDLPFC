#!/bin/bash
#$ -cwd
#$ -N "classify_nuclei_cart"
#$ -o ../../../processed-data/spot_deconvo/02-cellpose/08-classify_nuclei_cart_$TASK_ID.log
#$ -e ../../../processed-data/spot_deconvo/02-cellpose/08-classify_nuclei_cart_$TASK_ID.log
#$ -l mf=20G,h_vmem=20G
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

module load cellpose/2.0
python 08-classify_nuclei_cart.py

echo "**** Job ends ****"
date
