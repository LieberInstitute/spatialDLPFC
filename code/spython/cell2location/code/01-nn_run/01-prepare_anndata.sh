#!/bin/bash
#$ -cwd
#$ -N "prepare_anndata"
#$ -o ../../processed-data/01-nn_run/01-prepare_anndata.log
#$ -e ../../processed-data/01-nn_run/01-prepare_anndata.log
#$ -l bluejay,mf=30G,h_vmem=30G

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

module load cell2location/0.8a0
python 01-prepare_anndata.py

echo "**** Job ends ****"
date
