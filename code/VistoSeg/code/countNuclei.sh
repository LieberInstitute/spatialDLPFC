#!/bin/bash
#$ -cwd
#$ -l mem_free=25G,h_vmem=25G,h_fsize=100G
#$ -pe local 9
#$ -o Logs/V10B01-002_D1_Br2720_ant_DLPFC.txt
#$ -e Logs/V10B01-002_D1_Br2720_ant_DLPFC.txt
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

toolbox='/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/code/VistoSeg/code'
mask='/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/Images/round4/V10B01-002_D1_Br2720_ant_DLPFC_nuclei.mat'
jsonname='/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/NextSeq/Round4/DLPFC_Br2720_ant_2/outs/spatial/scalefactors_json.json'
posname='/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/NextSeq/Round4/DLPFC_Br2720_ant_2/outs/spatial/tissue_positions_list.csv'


matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), countNuclei('$mask','$jsonname','$posname')"

echo "**** Job ends ****"
date
