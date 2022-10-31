#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -pe local 4
#$ -o /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/code/VistoSeg/code/Logs/refineSegmentation.$TASK_ID.txt
#$ -e /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/code/VistoSeg/code/Logs/refineSegmentation.$TASK_ID.txt
#$ -m e
#$ -M madhavitippani28@gmail.com
#$ -t 19-30
#$ -tc 6


echo "**** Job starts ****"
date


echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"
echo "****"
echo "Sample id: $(cat /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/code/VistoSeg/code/ALLsamples.txt | awk '{print $1}' | awk "NR==${SGE_TASK_ID}") "
echo "****"

module load matlab/R2019a

toolbox='/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/code/VistoSeg/code'
fname=$(cat /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/code/VistoSeg/code/ALLsamples.txt | awk '{print $1}' | awk "NR==${SGE_TASK_ID}")
M=$(cat /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/code/VistoSeg/code/ALLsamples.txt | awk '{print $2}' | awk "NR==${SGE_TASK_ID}")

matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), refineVNS('$fname',$M)"

echo "**** Job ends ****"
date