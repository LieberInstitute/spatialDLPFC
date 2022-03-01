#!/bin/bash
#$ -pe local 5
#$ -l mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -o /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/code/VistoSeg/IF/splitSlide_IF.txt
#$ -e /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/code/VistoSeg/IF/splitSlide_IF.txt
#$ -m e
#$ -M madhavitippani28@gmail.com

 
echo "**** Job starts ****"
date
 
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"
echo "****"

toolbox='/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/code/VistoSeg/IF/'
fname='/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/Images/VisiumIF/InForm_0%/V10B01-087.mat'

matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), splitSlide_IF('$fname')" 
mv /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/Images/VisiumIF/InForm_0%/V10B01-087_*1* /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/Images/VisiumIF/VistoSeg/
echo "**** Job ends ****"
date

