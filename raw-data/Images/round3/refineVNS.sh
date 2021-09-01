#!/bin/bash
#$ -cwd
#$ -l mem_free=100G,h_vmem=100G,h_fsize=100G
#$ -o /dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/round3/logs/$TASK_ID.txt
#$ -e /dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/round3/logs/$TASK_ID.txt
#$ -t 1-12
#$ -tc 12
#$ -m a
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"
echo "****"
echo "Sample id: $(cat /dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/round3/DLPFC_round3_manual.txt | awk '{print $(NF)}' | awk "NR==${SGE_TASK_ID}")"
echo "****"

## Locate file and ids
FILE1=$(cat /dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/round3/DLPFC_round3_manual.txt | awk '{print $(NF-1)}'| awk "NR==${SGE_TASK_ID}")
M=$(cat /dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/round3/DLPFC_round3_manual.txt | awk '{print $(NF)}' | awk "NR==${SGE_TASK_ID}")
Code='/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/VisiumLIBD/code'
matlab -nodesktop -nosplash -r "addpath(genpath('$Code')), refineVNS('$FILE1',$M)"


echo "**** Job ends ****"
date
