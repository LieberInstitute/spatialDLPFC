#!/bin/bash
#$ -cwd
#$ -N "four_sample_demo"
#$ -j y
#$ -o /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/code/spython/tangram_libd/processed-data/02_four_sample_demo/four_sample_demo_$TASK_ID.log
#$ -e /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/code/spython/tangram_libd/processed-data/02_four_sample_demo/four_sample_demo_$TASK_ID.log
#$ -l gpu,mf=64G,h_vmem=64G
#$ -t 1-4
#$ -tc 1

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

conda activate tangram

conda env list

module list

python code/02_four_sample_demo/02_four_sample_demo.py -i $SGE_TASK_ID

echo "**** Job ends ****"
date
