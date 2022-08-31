#!/bin/bash
#$ -cwd
#$ -N "cart"
#$ -o ../../../processed-data/spot_deconvo/02-cellpose/07-cart.log
#$ -e ../../../processed-data/spot_deconvo/02-cellpose/07-cart.log
#$ -l mf=10G,h_vmem=10G

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

module load cellpose/2.0
python 07-cart.py

echo "**** Job ends ****"
date
