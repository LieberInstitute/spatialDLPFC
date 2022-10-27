#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -pe local 4
#$ -N spatialDLPFC_rerun_countNuclei
#$ -o Logs/countNuclei.$TASK_ID.txt
#$ -e Logs/countNuclei.$TASK_ID.txt
#$ -m e
#$ -t 1-12
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

module load matlab/R2019a

## Read parameters
SAMPLE=$(awk 'BEGIN {FS="\t"} {print $1}' ../../spaceranger/spaceranger_parameters.txt | awk "NR==${SGE_TASK_ID}")
IMAGEPATH=$(awk 'BEGIN {FS="\t"} {print $4}' ../../spaceranger/spaceranger_parameters.txt | awk "NR==${SGE_TASK_ID}")
MASKPATH=$(echo ${IMAGEPATH} | sed "s/.tif/_nuclei_WS_final.mat/g")

echo "Processing sample ${SAMPLE} with mask ${MASKPATH}"
date

toolbox="/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/code/VistoSeg/code"
jsonname="/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rerun_spaceranger/${SAMPLE}/outs/spatial/scalefactors_json.json"
posname="/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rerun_spaceranger/${SAMPLE}/outs/spatial/tissue_positions_list.csv"


matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), countNuclei('${MASKPATH}','$jsonname','$posname')"

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/
