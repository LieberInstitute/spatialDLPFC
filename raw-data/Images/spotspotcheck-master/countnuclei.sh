#!/bin/bash
#$ -cwd
#$ -l mem_free=200G,h_vmem=50G,h_fsize=100G
#$ -o /dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/spotspotcheck-master/logs/$TASK_ID.txt
#$ -e /dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/spotspotcheck-master/logs/$TASK_ID.txt
#$ -m e
#$ -M heenadivecha@gmail.com
#$ -t 1


 
echo "**** Job starts ****"
date
 
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"
echo "****"
echo "Sample id: $(countnuclei)"
echo "****"



matlab -nodesktop -nosplash -r "addpath(genpath('$toolbox')); rnascope_mouse('$FILE1','$toolbox')"

echo "**** Job ends ****"
date
