#!/bin/bash
#$ -cwd
#$ -o logs/subset_SZBDMulti.log
#$ -e logs/subset_SZBDMulti.log
#$ -l caracol,mf=150G,h_vmem=150G,h_fsize=50G

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

module load tangram/1.0.2
python subset_SZBDMulti.py

echo "**** Job ends ****"
date
