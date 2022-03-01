#!/bin/bash
#$ -pe local 6
#$ -l mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -o /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/code/VistoSeg/IF/InFormStitch.log
#$ -e /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/code/VistoSeg/IF/InFormStitch.log
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

module load matlab/R2019a

toolbox='/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/code/VistoSeg/IF'

filename='/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/Images/VisiumIF/InForm_0%/20220126_VIF_DLPFC_Real_Scan1_*_component_data.tif'
fname='V10B01-087'
fname2='/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/Images/VisiumIF/InForm_0%/20220126_VIF_DLPFC_Real_Scan1_[14415,65835]_component_data.tif'
matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), O = extractMD('$fname2'); InFormStitch('$filename',O,6,'$fname')"

echo "**** Job ends ****"
date
