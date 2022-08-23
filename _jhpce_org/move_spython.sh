#!/bin/bash
#$ -cwd
#$ -l mem_free=10G,h_vmem=10G,h_fsize=300G
#$ -N move_spython
#$ -o logs/move_spython.txt
#$ -e logs/move_spython.txt
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

## Edit with your job command
rsync -avh /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spython /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/code

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/
