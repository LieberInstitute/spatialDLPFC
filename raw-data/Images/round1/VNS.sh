#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=200G,h_vmem=200G,h_fsize=200G
#$ -o /dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/VNSlogs/$TASK_ID.txt
#$ -e /dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/VNSlogs/$TASK_ID.txt
#$ -m e
#$ -M madhavitippani28@gmail.com
#$ -t 5
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
echo "Sample id: $(cat /dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/listOfFiles.txt | awk '{print $NF}' | awk "NR==${SGE_TASK_ID}")"
echo "****"

FILE1=$(cat /dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/listOfFiles.txt | awk '{print $NF}' | awk "NR==${SGE_TASK_ID}")
echo $FILE1

cd /dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images

matlab -nodesktop -nosplash -r "VNS('$FILE1')"

echo "**** Job ends ****"
date

