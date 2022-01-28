#!/bin/bash
#$ -cwd
#$ -N "run_spagcn_02"
#$ -o ../../processed-data/04-spatial_dlpfc/02-run_spagcn_$TASK_ID.log
#$ -e ../../processed-data/04-spatial_dlpfc/02-run_spagcn_$TASK_ID.log
#$ -l caracol,mf=60G,h_vmem=60G
#$ -t 3-30
#$ -tc 2

image_info=../../raw-data/04-spatial_dlpfc/sample_info.txt

image_path=$(awk "NR == $SGE_TASK_ID" $image_info | cut -f 4)
sample_name=$(awk "NR == $SGE_TASK_ID" $image_info | cut -f 1 | cut -d "_" -f 2-4 | sed -r 's/_(manual|extra)//')

module load spagcn/1.2.0
python 02-run_spagcn.py -i $image_path -s $sample_name
