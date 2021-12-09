#!/bin/bash
#$ -cwd
#$ -l mem_free=20G,h_vmem=20G,h_fsize=100G
#$ -pe local 8
#$ -o /users/hdivecha/VistoSeg/code/Logs/V10B01-002_D1_Br2720_ant_DLPFC_VNS.txt
#$ -e /users/hdivecha/VistoSeg/code/Logs/V10B01-002_D1_Br2720_ant_DLPFC_VNS.txt
#$ -m e
#$ -M heenadivecha@gmail.com


echo "**** Job starts ****"
date


echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"
echo "****"
echo "Sample id: V10B01-002_D1_Br2720_ant_DLPFC"
echo "****"

module load matlab/R2019a

toolbox='/users/hdivecha/VistoSeg/code'
fname='/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/Images/round4/V10B01-002_D1_Br2720_ant_DLPFC.tif'

matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), VNS('$fname',5)"

echo "**** Job ends ****"
date
