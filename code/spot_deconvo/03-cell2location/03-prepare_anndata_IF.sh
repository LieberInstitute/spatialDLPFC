#!/bin/bash
#$ -cwd
#$ -N "prepare_anndata_IF"
#$ -o ../../../processed-data/spot_deconvo/03-cell2location/03-prepare_anndata_IF.log
#$ -e ../../../processed-data/spot_deconvo/03-cell2location/03-prepare_anndata_IF.log
#$ -l bluejay,mf=20G,h_vmem=20G,h_fsize=50G

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

module load tangram/1.0.2
python 03-prepare_anndata_IF.py

echo "**** Job ends ****"
date
