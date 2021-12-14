#!/bin/bash
#$ -cwd
#$ -l mem_free=25G,h_vmem=25G,h_fsize=100G
#$ -pe local 4
#$ -o Logs/DLPFC_Br3942_ant_manual_alignment.txt
#$ -e Logs/DLPFC_Br3942_ant_manual_alignment.txt
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
echo "Sample id: DLPFC_Br3942_ant_manual_alignment"
echo "****"

module load matlab/R2019a

toolbox='/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/code/VistoSeg/code'
mask='/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/Images/round1/Liebert_Institute_OTS-20-7690_rush_anterior_2_nuclei.mat'
jsonname='/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rerun_spaceranger/DLPFC_Br3942_ant_manual_alignment/outs/spatial/scalefactors_json.json'
posname='/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rerun_spaceranger/DLPFC_Br3942_ant_manual_alignment/outs/spatial/tissue_positions_list.csv'


matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), countNuclei('$mask','$jsonname','$posname')"

echo "**** Job ends ****"
date
