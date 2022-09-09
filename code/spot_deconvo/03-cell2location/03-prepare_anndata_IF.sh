#!/bin/bash
#$ -cwd
#$ -N "prepare_anndata_IF"
#$ -o ../../../processed-data/spot_deconvo/03-cell2location/03-prepare_anndata_IF_broad.log
#$ -e ../../../processed-data/spot_deconvo/03-cell2location/03-prepare_anndata_IF_broad.log
#$ -l mf=80G,h_vmem=80G,h_fsize=50G

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

module load cell2location/0.8a0
python 03-prepare_anndata_IF.py

echo "**** Job ends ****"
date
